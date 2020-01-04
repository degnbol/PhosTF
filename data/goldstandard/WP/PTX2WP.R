# packages
library(Matrix)
# functions
flatten = function(x) as.vector(as.matrix(x))

setwd("~/cwd/data/goldstandard/WP")

P_edges = read.table("../P/detectable/edges.tsv", sep="\t", quote="", header=T)
T_ = flatten(read.table("../T/T.txt", quote=""))
P  = flatten(read.table("../P/detectable/P.txt", quote=""))

row_names = c(P, T_)
col_names = P

i = match(P_edges$Target, row_names)
j = match(P_edges$P, col_names)
x = P_edges$Relationship
# make sure levels are in expected order
levels(x)

# has to be numeric, we convert back when dense. factors in x are reduced to 1, 2, 3
dimensions = c(length(row_names), length(col_names))
dimension_names = list(row_names, col_names)
W = as.matrix(sparseMatrix(as.numeric(i), as.numeric(j), x=as.numeric(x), dims=dimensions, dimnames=dimension_names))
W[W == 2] = "+"
W[W == 3] = "-"

# check things
sum(W == 1)
sum(W == "+")
sum(W == "-")
dim(W)

write.table(W, file="WP.mat", sep=" ", quote=F)

