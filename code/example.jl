include("mtx.jl") #import functions from mtx.jl
c=@__DIR__ #get current dir
c=c*"/../data/example/" #concatenate current dir and path to example folder

println("Loading example matrix, healthy sample, 2000 genes, 500 cells")
x = MBG(c*"healthy_2000_500.h5","hg19") # load example h5 file, human species sample (hg19)
println("Matrix loaded, matrix shape: "*x.M.m*", "*x.M.n)
println("Running svd")
s=svds(x.M,nsv=3)[1];
println("SVD done")
