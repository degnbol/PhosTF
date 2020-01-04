# packages
library(Matrix)
# functions
flatten = function(x) as.vector(as.matrix(x))

setwd("~/cwd/data/goldstandard/WT")

T_edges = read.table("../T/edges.tsv", sep="\t", quote="", header=T)
T_ = flatten(read.table("../T/T.txt", quote=""))
P  = flatten(read.table("../P/detectable/P.txt", quote=""))
X  = flatten(read.table("../X/detectable/X.txt", quote=""))

row_names = c(P, T_, X)
col_names = T_

i = match(T_edges$Target, row_names)
j = match(T_edges$TF, col_names)
x = T_edges$Mode
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

write.table(W, file="WT.mat", sep=" ", quote=F)

