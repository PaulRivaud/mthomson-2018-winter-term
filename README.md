# Introduction

This Github repository contains data files and code examples to get you started in the Computational Biology class taught by Matt Thomson at Caltech. The goal of this course is to have the students analyze some single-cell sequencing data using programming languages (Julia, R, Python, Matlab, etc.). The content of this README is listed below:  
  
[Basic information regarding data files](#basic-information-regarding-data-files)  
__JULIA__  
[Installing Julia (and Juno)](#installing-julia-and-juno)  
[Installing packages in Julia](#installing-packages-in-julia)  
[Running the julia example](#running-the-julia-example)  
[Loading mtx files in Julia](#loading-mtx-files-in-julia)  
[Loading labels in Julia](#loading-labels-in-julia)  
__R__  
[Loading mtx files in R](#loading-mtx-files-in-r)  
[Loading labels in R](#loading-labels-in-r)  
__PYTHON__  
[Loading mtx files in Python](#loading-mtx-files-in-python)  
[Loading labels in Python](#loading-labels-in-python)  
__MATLAB__  
[Loading mtx files in Matlab](#loading-mtx-files-in-matlab)  
[Loading labels in Matlab](#loading-labels-in-matlab)  

# Basic information regarding data files  

Single-cell sequencing data is often stored as sparse matrix objects to cope with the data low density (~10-15% of the entries are non-zero entries). Working with sparse matrices is computationally more efficiency but requires more rigor when it comes down to keeping track of row and column labels, which are stored in separate arrays (**Reminder: Array indexing starts at 0 in Python, 1 in Julia, R and Matlab**). Dataframe objects ie dense matrices (R, Python) deal with that aspect but tend to perform slower and can be hard to load in memory when the datasets get larger.  
  
Once unziped, the data folders contain three files:  
-- `matrix.mtx`: the read values of the gene expression matrix and their respective row and column indices, in a matrix market format.  
-- `barcodes.tsv`: a file containing the column (cell) labels.  
-- `genes.tsv`: a file containing the row (gene) labels.  

* Sparse Matrices (MatrixMarket/.mtx format):  
[Info here](https://math.nist.gov/MatrixMarket/formats.html#MMformat)  
Mtx files store MatrixMarket format matrices. The base principle is to store the row and column indices of each non-zero entry in the matrix. The MatrixMarket format contains three distinct parts:  
    - Comment lines that start with `%`
    - Header line: total number of rows, total number of columns, total number of non-zeros entries (space separated)
    - Entries: row index, column index, entry value (space separated)  

  The gene and barcode labels are stored in separate files and need to be read separately.  

* Sparse matrices (MatrixMarket deconstructed into multiple vectors):  
[Info here (10X's HDF5 format)](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices)  
10X Genomics also deconstructs sparse matrices into multiple vectors: the sparse matrix is a Matrix Market object in RAM when used, but is not stored as one under their HDF5 format. Storing all these vectors under one file enables users to not lose track of matrices respective's labels. The different vectors are:
    - `data (or nzval)`: non-zero entry values (length: number of entries)
    - `indices (or rowval)`: row indices (length: number of entries)
    - `indptr (or colptr)`: column index pointers. Index of the start of each column (length: number of columns +1, the last value indicates the end index of the last column +1).
    - `barcodes`: barcode labels (length: number of columns)
    - `gene_names`: gene labels (length: number of rows)  
We can illustrate this concept with a Julia example. Lets' consider the following matrix `M_dense`:  
<pre><code>M_dense
5×3 Array{Int64,2}:  
1000  1200    0  
   0     0    2  
   2   400    0  
   0     0  500  
   0     0    0  
</pre></code>  
The matching sparse matrix `M` is:
<pre><code>M
5×3 SparseMatrixCSC{Int64,Int64} with 6 stored entries:
  [1, 1]  =  1000
  [3, 1]  =  2
  [1, 2]  =  1200
  [3, 2]  =  400
  [2, 3]  =  2
  [4, 3]  =  500
</pre></code> 
The deconstructed vectors of M are:
<pre><code>M.nzval
Int64[6]
1000
2
1200
400
2
500

M.rowval 
Int64[6]
1
3
1
3
2
4

M.colptr 
Int64[4]
1
3
5
7
</pre></code>


# Julia  

## Installing Julia (and Juno)
Visit the JuliaComputing web page https://juliacomputing.com/products/juliapro.html to navigate your way to the download section. You may be asked to enter an email address to have access to a free download.

* Linux:  
Select `JuliaPro-0.6.2.1 - MKL (for Linux)`, or the latest release. Once the download is complete, just uncompress the archive and add the bin folder to the PATH variable in the .bashrc file or from the terminal. The command should look similar to the following:

<pre><code>export PATH=$PATH:/home/user/Downloads/julia-235245113/bin</code></pre>

* Mac OSX:  
Select `JuliaPro-0.6.2.1-MKL (for Mac)`, or the latest release. Once the download is complete, uncompress the archive and go through the installation process.

Those installations will install Julia (command line version) and should install Juno as well, a Julia IDE based on Atom. If you have trouble installing Juno, we recommend to look at this page: https://github.com/JunoLab/uber-juno/blob/master/setup.md.

## Installing packages in Julia
Packages are installed thanks to Julia commands (https://docs.julialang.org/en/stable/stdlib/pkg/). The most common way to install and use packages is:
<pre><code>Pkg.add("packageName")  
using packageName</pre></code>  

## Running the Julia example
Prior to running the example featured in that repository, run the install_packages.jl script using the following:
<pre><code>julia install_packages.jl</pre></code>  
Look up the [Installing packages](#installing-packages) section to learn more about package management.

## Loading mtx files in Julia
Once the MatrixMarket package installed, you can use the following:
<pre><code>"""
    read_csc(pathM)
Read an mtx file pointed to by `pathM` and return a SparseMatrixCSC object.
"""
function read_csc(pathM::String)
     x=MatrixMarket.mmread(pathM);
     Float64.(x)
end
</code></pre>
Note:  
Strings in Julia must be delimited by double quotes (`"`), not single quotes (`'`).  
The dot used after `Float64` applies [broadcasting](https://docs.julialang.org/en/stable/manual/arrays/#Broadcasting-1). It enables an operation to be applied to every entry in an object (it is similar to mapping).

## Loading labels in Julia
<pre><code>"""
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
</code></pre>

<pre><code>"""
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
</code></pre>

# R  

## Loading mtx files in R
R base package [`Matrix`](https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/externalFormats.html) lets you load .mtx files directly. The returned object is a [dgTMatrix](https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/dgTMatrix-class.html).You can load a file as shown below:  
<pre><code>> readMM("path_to_matrix/matrix.mtx")  
32738 x 1985 sparse Matrix of class "dgTMatrix"</pre></code> 

Note: Even though Matrix is a base R package, it has to be loaded manually (through the `Packages` section in R studio or `library("Matrix")` in the R console.)

## Loading labels in R  
acrwea

# Python  

## Loading mtx files in Python
Mtx files can be loaded in Python using the `io` module of Scipy. There are multiple ways to install Scipy (Homebrew, pip, conda, etc.). For example, type in your terminal:
<pre><code>pip install scipy</pre></code>  
You can then use in Python:
<pre><code>import scipy.io
M = scipy.io.mmread('path_to_matrix/matrix.mtx')
</pre></code>  
M is a [coo_matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html) object, you can use M.[tocsc](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.tocsc.html#scipy.sparse.coo_matrix.tocsc)() to convert the coo_matrix to a [csc_matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html). The latter sorts the columns in the object which facilitates column-based operations.

## Loading labels in Python  
The [csv](https://docs.python.org/2/library/csv.html) module will help you read cell and gene labels from the barcodes.tsv and genes.tsv files.
<pre><code>barcodes = [row[0] for row in csv.reader(open('my_path/barcodes.tsv'), delimiter="\t")]  
genes =row[1].upper() for row in csv.reader(open('my_path/genes.tsv'), delimiter="\t")]
</pre></code>  

Note that you can also use very basic file processing if you find it easier:
<pre><code>barcodes = []
with open('mypath/barcodes.tsv') as f:
    for line in f:
        barcodes.append(line.strip('\n'))
</pre></code>  

# Matlab

## Loading mtx files in Matlab  
ff

## Loading labels in Matlab
wefw
