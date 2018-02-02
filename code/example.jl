include("mtx.jl") #import functions from mtx.jl
c=@__DIR__ #get current dir
c=c*"/../data/example/" #

x = MBG(c*"healthy_2000_500.h5","hg19") # load example h5 file, human species sample (hg19)
println(x)