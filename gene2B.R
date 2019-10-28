#!/usr/bin/env Rscript
library("tools")

# a few functions copied from Eberhardt LLC paper code to get B from gene expression level
# assumes single knockouts with knockout in the diagonal of the read gene matrix

mpinv<-function(X) {
    # Moore-Penrose inverse (by Thorsten Hothorn, downloaded from
    # http://tolstoy.newcastle.edu.au/R/help/99a/1118.htm)
    #
    # Example usage:
    # X <- cbind(1, diag(3)); # singular matrix
    # y <- 1:3
    # b2 <- mpinv(X)%*%y # least square fit using Moore-Penrose
    # X%*%b2 # == y 	
    #
    
    Eps <- 100 * .Machine$double.eps
    
    # singular value decomposition
    s <- svd(X)
    d <- s$d
    m <- length(d)
    if (!(is.vector(d)))
        return(t(s$v%*%(1/d)%*%t(s$u)))
    # remove eigenvalues equal zero
    d <- d[d > Eps]
    notnull <- length(d) 
    if (notnull == 1)
    {
        inv <- 1/d
    } else {
        #browser()
        #inv <- solve(diag(d))
        inv <- diag(1/d) #correction made by A. Hyttinen
    }
    # add rows, columns of zeros if needed 
    if (notnull != m)
    {
        inv <- cbind(inv, matrix(0, nrow=notnull, ncol=(m - notnull)))
        inv <- rbind(inv, matrix(0, nrow=(m-notnull), ncol=m))
    } 
    
    # compute Moore-Penrose
    mp <- s$v%*%inv%*%t(s$u)
    
    # set very small values to zero
    mp[abs(mp) < Eps] <- 0
    return(mp)
}

estimate.L0<-function(eqs, maxpa=3, lambda=0.01) {
    #L0 regularized solution of a linear equations system.
    #INPUT:
    # eqs - System of equations eqs$K%*%b =eqs$k.
    # maxpa - Maximum number of parents considered.
    # lambda - Regularization penalizing the number of parents.
    #OUTPUT:
    # b - solution vector or B directs effects matrix if 
    #     indexing structure eqs$P given
    
    n<-eqs$n
    
    #insert the solution with maxpa=0
    best.SSE<-sum(eqs$k^2)
    b<-rep(0,ncol(eqs$K))
    
    for (npa  in 1:maxpa ) {
        pasets<-combinations(ncol(eqs$K),maxpa) 
        
        for ( ipa in 1:nrow(pasets) ) {
            #sk<-eqs$k[index.k]
            sK<-eqs$K[,pasets[ipa,],drop=FALSE]
            estimate<-mpinv(t(sK)%*%sK)%*%t(sK)%*%eqs$k
            
            #evaluation criterion
            SSE<-sum( ( sK%*%estimate - eqs$k )^2 ) + lambda*npa
            
            
            if ( SSE < best.SSE ) {
                print(SSE)
                
                #print(paset)
                best.SSE<-SSE
                b<-rep(0,ncol(eqs$K))
                b[pasets[ipa,]]<-estimate
            }
        }
    }
    if ( all( !is.na( eqs$P ) ) ) {
        b.to.B(b, eqs$P, eqs$n )
    } else {
        b
    }
}

estimate.L1<-function(eqs, lambda=0.01, prev=NA) {
    #L1 regularized solution of a linear equations system.
    #INPUT:
    # eqs - System of equations eqs$K%*%b =eqs$k.
    # lambda - Regularization penalizing the number of parents.
    # prev - 
    #OUTPUT:
    # b - solution vector or B directs effects matrix if 
    #     indexing structure eqs$P given
    
    
    gamma<-2
    
    #function to be minimized
    f<-function(b) {
        sum( (eqs$K%*%b-eqs$k)^2 ) + lambda*(1/gamma)*sum( log(cosh(gamma*b)) )
    }
    
    #gradient
    g<-function(b) {
        as.vector(2*t(eqs$K)%*%(eqs$K%*%b-eqs$k) + lambda*tanh(gamma*b)) #notice that the 2 disappears
    }
    
    #determining the starting value
    #this is a fairly good estimate better to use some regularization already here
    if ( !all(is.na(prev)) ) {
        b0<-prev[t(eqs$P)]
    } else {
        b0<-as.vector( mpinv(t(eqs$K)%*%eqs$K+lambda*diag(ncol(eqs$K)))%*%t(eqs$K)%*%eqs$k )
    }
    
    if ( is.na(f(b0)) || is.infinite(f(b0)) ) {
        b0<-rep(0,ncol(eqs$K))
    }
    
    opt<-optim(b0,f,gr=g,method ='BFGS',control=list(maxit=10000,reltol=1e-15))
    b<-opt$par
    
    if (all(!is.na(eqs$P))) {
        b.to.B(b, eqs$P, eqs$n)
    } else {
        b
    }
}

estimate.L2<-function(eqs, lambda=0.01) {
    #L2 regularized solution of a linear equations system.
    #INPUT:
    # eqs - System of equations eqs$K%*%b =eqs$k.
    # lambda - Regularization penalizing the number of parents.
    #OUTPUT:
    # b - solution vector or B directs effects matrix if 
    #     indexing structure eqs$P given
    
    b<-as.vector( mpinv(t(eqs$K)%*%eqs$K+lambda*diag(ncol(eqs$K)))%*%t(eqs$K)%*%eqs$k )
    
    if (all(!is.na(eqs$P))) {
        b.to.B(b, eqs$P, eqs$n )
    } else {
        b
    }
} 

gene2T<-function(gene) {
    # create T matrix of all total causal effects for a table of gene expression measurements
    n = ncol(gene)
    diagonal = matrix(rep(diag(as.matrix(gene)), n), nrow=n, byrow=TRUE)
    gene[1:n,] / diagonal
}

T2B<-function(ntotals, lambda=0.01) {
    n = nrow(ntotals)
    B<-array(0,c(n,n))
    # have to process here only arcs to one node at a time
    for (node in 1:n) {
        # Creating the equation structure here.
        # This is very easy since all observed quantaties are totals.
        # xi intervened, quite simple equations!
        # k becomes a row of "T" without the diagonal since it is the knockout
        k<- ntotals[node, -node]
        # every row in K is a column of "T" without the diagonal, so "T" without diag transposed
        K<-t(ntotals[-node,-node])
        ### Solving the system
        # cat("solving for node", node, '\n')
        B[node,-node]<-estimate.L1(list(n=n, P=NA, K=K, k=as.double(k)), lambda=lambda)
    }
    B
}


read = function(path) {
    if (file_ext(path) == "csv") {return(read.csv(path, header=T, row.names=1))}
    else {return(read.table(path))}
}


# arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
    stop("Supply two arguments: gene infile, B outfile.", call.=FALSE)
}
gene.fname = args[1]
outfname = args[2]

cat("Reading", gene.fname, "\n")
gene = read(gene.fname)
cat("gene -> T\n")
ntotals = gene2T(gene)
cat("T -> B\n")
B = T2B(ntotals, lambda=0.0001)
cat("Writing to", outfname, "\n")
write.table(B, outfname, sep=" ", row.names=F, col.names=F)



