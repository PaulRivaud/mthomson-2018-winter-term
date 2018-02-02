using MatrixMarket
using JLD
using MultivariateStats
using HDF5
using MAT

"""
    MBG
Structure that contains:
`M`: Sparse Matrix
`B`: Column names (columns)
`G`: Gene names (rows)
`GE`: Gene Emsembl IDs
"""
mutable struct MBG{T<:SparseMatrixCSC}
    M::T #data,SparseMatrixCSC
    B::Vector{String} #barcode names
    G::Vector{String} #Gene symbols (names)
    GE::Vector{String} #Ensembl gene IDs
end

"""
    MBG(pathM,pathB,pathG)
Build an MBG struct from an mtx file `pathM`,
a tsv file `pathB` containing the barcodes,
a tsv file `pathG` containing the genes.
Return MBG struct x.
"""
MBG(pathM::String,pathB::String,pathG::String) = MBG(read_csc(pathM),read_barcodes(pathB),read_genes(pathG),read_genes_ensembl(pathG))

"""
    MBG(x)
Copy an existing MBG object `x`
"""
MBG(x::MBG) = deepcopy(x)

"""
    MBG(p,species)
Build an MBG struct from an H5 file `p`.
Objects are accessed internally using `species` in their respective path.
Return MBG struct x.
"""
function MBG(p::String,species::String="mm10")
    f=h5open(p)
    D=read(f[species*"/data"])#::Vector{Int32} #count values
    GI=read(f[species*"/indices"])::Vector{Int64} #row indices
    G=read(f[species*"/gene_names"])::Vector{String} #gene names
    GE=read(f[species*"/genes"])::Vector{String} #gene names
    B=read(f[species*"/barcodes"]) #barcode names
    IP=read(f[species*"/indptr"])::Vector{Int64} #col pointers
    S=read(f[species*"/shape"])::Vector{Int32} #shape tuple (m,n) m rows n cols
    GI+=1; #Julia indexing starts at 1
    IP+=1; #Julia indexing starts at 1
    G = [uppercase(x) for x in G] #meet uppercase standards for gene names
    #if the gene indices are sorted in descending order for each cell (true for h5 files output by the cellranger count pipeline)
    #this is necessary to properly build a SparseMatrixCSC obj in Julia
    if GI[1]>GI[2]
        for col in 1:S[2]
            a = IP[col] #a and b define a column (a block in the matrix), hence the start and end indices of the block
            b = IP[col+1]-1
            GI[a:b] = GI[b:-1:a] #reversing the gene indices for column col
            D[a:b] = D[b:-1:a] #reversing the matching data values
        end
    end
    #construct sparse matrix with shape, index pointers, row indices and non-zero values
    M=SparseMatrixCSC(S[1],S[2],IP,GI,D);
    #feeding the MBG struct
    x=MBG(Float64.(M),B,G,GE)
end

"""
    read_csc(pathM)
Read an mtx file pointed to by `pathM` and return a SparseMatrixCSC object.
"""
function read_csc(pathM::String)
    x=MatrixMarket.mmread(pathM);
    Float64.(x)
end

"""
    read_barcodes(tsvPath)
Read a tsv file and return its values in an array of strings.
"""
function read_barcodes(tsvPath::String)
    f=open(tsvPath)
    lines=readlines(f)
    a=String[]
    for l in lines
        push!(a,uppercase(l))
    end
    close(f)
    return a
end

"""
    read_barcodes(tsvPath, sample)
Read a tsv file and return its values concatenated with the `sample` name in an array of strings.
"""
function read_barcodes(tsvPath::String,sample::String)
    f=open(tsvPath)
    lines=readlines(f)
    a=String[]
    for l in lines
        push!(a,sample*"_"*uppercase(l))
    end
    close(f)
    return a
end

"""
    read_genes(tsvPath)
Read a tsv file and return its values in an array of strings.
"""
function read_genes(tsvPath::String)
    f=open(tsvPath)
    lines=readlines(f)
    a=String[] #Array{}
    for l in lines
        push!(a,uppercase(split(l,"\t")[2]))
    end
    close(f)
    return a
end

"""
    read_genes_ensembl(tsvPath)
Read a tsv file and return its values in an array of strings.
"""
function read_genes_ensembl(tsvPath::String)
    f=open(tsvPath)
    lines=readlines(f)
    a=String[] #Array{}
    for l in lines
        push!(a,uppercase(split(l,"\t")[1]))
    end
    close(f)
    return a
end

"""
    save_MBG(p,x,species)
Save the elements of an MBG object to an h5 file.
`p`: path of the H5 file to be written.
`x`: MBG object.
`species`: experiment species, used to build the H5 tree structure. Default is "mm10".
"""
function save_MBG(p::String,x::MBG,species="mm10"::String)
    #to hdf5
    if !endswith(p,".h5")
        p*=".h5";
    end
    h5open(p,"w") do file
        write(file,species*"/barcodes",x.B)
        write(file,species*"/gene_names",x.G)
        write(file,species*"/genes",x.GE)
        write(file,species*"/data",x.M.nzval)
        write(file,species*"/indices",x.M.rowval)
        write(file,species*"/indptr",x.M.colptr)
        write(file,species*"/shape",[x.M.m,x.M.n])
    end
end

"""
    save_stds(p,s,species)
Save the standard deviations array `s` to an h5 file.
`p`: path of the H5 file to be written.
`s`: Vector{Float64}
`species`: experiment species, used to build the H5 tree structure. Default is "mm10".
"""
function save_stds(p::String,s::Vector{Float64},species="mm10"::String)
    #to hdf5
    if !endswith(p,".h5")
        p*=".h5";
    end
    h5open(p,"w") do file
        write(file,species*"/stds",s)
    end
end

"""
    catSparseMatrices(V)
Concatenate the sparse matrices of vector `V`.
Return SparseMatrixCSC R.
"""
function catSparseMatrices(V::Vector{SparseMatrixCSC{Float64,Int64}})
    n=0
    for M in V
        n += M.n #update total number of columns
    end
    R = spzeros(Float64,Int64,V[1].m,n) #empty sparse matrix with m rows n cols
    s=1 #initial start
    e=0 #initial end (mock init though)
    for M in V
        e+=M.n #update end position
        R[:,s:e]=M #feed values
        s=e+1 #update start position
    end
    return R
end

"""
    catMBGs(V)
Concatenate the MBGs of vector `V`/
Return MBG struct.
"""
function catMBGs(V::Vector{MBG{SparseMatrixCSC{Float64,Int64}}})
    M=[x.M for x in V] #vector of SparseMatrixCSC
    B=[x.B for x in V] #vector of Vector{String}
    #MBG constructor with cat-ed sparse matrices, cat-ed barcode names, one gene list and one ensembl gene linst
    MBG(catSparseMatrices(M),vcat(B...),V[1].G,V[1].GE)
end

"""
    catH5s(V)
Concatenate the matrices from the H5 files of `V`.
`V` contains string paths only.
Return MBG struct.
"""
catH5s(V::Vector{String},species="mm10"::String) = catMBGs([MBG(P,species) for P in V])

"""
    add_sample_name!(x,name)
Add a sample `name` to each barcode of `x`.
Example: with `name` set to "N", ["BC1","BC2"] becomes ["N_BC1","N_BC2"].
Overwrites barcode names vector.
"""
function add_sample_name!(x::MBG,name::String)
    name=name*"_";
    x.B = name.*x.B;
end

"""
    normalizecols!(X, p)
Normalize the non-empty columns of `X` using a `p`-norm.
Overwrite the content of `X`.
"""
function normalize_cols!(X::SparseMatrixCSC,m=100::Int64)
    for col = 1:X.n
        rr = X.colptr[col]:X.colptr[col+1]-1 #pointer range for column col
        s = sum(@view(X.nzval[rr])) #grab non zeros values in that range and compute associated sum
        X.nzval[rr] .= X.nzval[rr] .* m ./ s #divide no zeros of X within rr by n

    end
    X #conventional return (not necessary)
end

"""
    normalizerows!(X)
Normalize the rows of `X` using their respective standard deviation values.
std=sqrt((1/(n-1)).sum((xi-mean)^2)).
Overwrite the content of `X`.
"""
function normalize_rows!(x::SparseMatrixCSC)
    m=zeros(x.m) #array of means (m rows)
    s=zeros(x.m) #array of stds
    #compute means
    for (i,v) in zip(x.rowval,x.nzval) #ith row index, vth value
        m[i]+=v
    end
    m./=x.n #divide by number of columns
    #compute stds
    for (i,v) in zip(x.rowval,x.nzval)
        s[i]+=(v-m[i])^2
    end
    s./=x.n-1
    s=sqrt.(s)
    #println(s)
    #divide by std
    for (k,i) in enumerate(x.rowval) #kth entry, ith row index
        x.nzval[k]/=s[i]
    end
    return s
end

"""
    round!(x,dec)
Round the non-zero values of a SparseMatrixCSC object `x`,
down to `dec` decimals.
Overwrite `x`.nzval.
"""
function round!(x::SparseMatrixCSC,dec=4::Int64)
    for i = 1:length(x.nzval)
        x.nzval[i]=round(x.nzval[i],dec)
    end
end

"""
    norm_sd!(x)
Apply column and row normalization to the data of MBG object `x`.
Overwrite the values of `x`.M.
`m` is the factor used for the column normalization after the sum division.
!! deal with saving stds
!! deal with trimming
"""
function norm_sd!(x::MBG,m=100::Int64)
    ii = cutoff(x.M)
    normalize_cols!(x.M,m)
    select_genes!(x,ii)
    s=normalize_rows!(x.M)
    round!(x.M)
    println(s)
end

"""
    cutoff(x,m)
Find rows (genes) of `x which max value is equal or greater than `m`.
Return an array of indices that satisfy this condition.
"""
function cutoff(x::SparseMatrixCSC,m=5::Int64)
    vmax=vec(maximum(x,2)) #get max of each row (axis 2)
    find(x->(x>=m),vmax) #return indices of rows where the condition is met
end

"""
    select_genes!(x,ii)
Splice sparse matrix `x` with row indices vector `ii`.
Ordering of `ii` is preserved.
Overwrite `x`.
"""
function select_genes!(x::MBG,ii::Vector{Int64})
    x.M=x.M[ii,:]
    x.G=x.G[ii]
    x.GE=x.GE[ii]
end

"""
    select_genes!(x,gn)
Splice sparse matrix `x` with gene names vector `gn`.
Ordering of `gn` is preserved.
Overwrite `x`.
"""
function select_genes!(x::MBG,gn::Vector{String})
    gn = [uppercase(x) for x in gn]
    ii=sfindin(x.G,gn) #find indices in x.G of elements in gn
    select_genes!(x,ii)
end

"""
    select_cells!(x,ii)
Splice sparse matrix `x` with column indices vector `ii`.
Ordering of `ii` is preserved.
Overwrite `x`.
"""
function select_cells!(x::MBG,ii::Vector{Int64})
    x.M=x.M[:,ii]
    x.B=x.B[ii]
end

"""
    select_cells!(x,cn)
Splice sparse matrix `x` with cell names vector `cn`.
Ordering of `cn` is preserved.
Overwrite `x`.
"""
function select_cells!(x::MBG,cn::Vector{String})
    ii=sfindin(x.B,cn)
    select_cells!(x,ii)
end

"""
    sfindin(a,b)
Sorted version of the Base.findin function.
Return the indices of elements in a that
"""
function sfindin(a, b)
    ind = Array{Int,1}(0)
    adict=Dict();
    for (i,ai) in enumerate(a)
        adict[ai]=i
    end
    @inbounds for g in b
        push!(ind, adict[g])
    end
    ind
end
