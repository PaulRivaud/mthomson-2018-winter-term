# Introduction

This Github repository contains data files and code examples to get you started in the Computational Biology class taught by Matt Thomson at Caltech. The goal of this course is to have the students analyze some single-cell sequencing data using programming languages (Julia, R, Python, Matlab, etc.). The content of this README is listed below:  
  
[Basic information regarding data files](#basic-information-regarding-data-files)  
__Julia__  
[Installing Julia (and Juno)](#installing-julia-and-juno)  
[Installing packages in Julia](#installing-packages-in-julia)  
[Running the julia example](#running-the-julia-example)  
__R__  
[Loading mtx files in R](#loading-mtx-files-in-r)  
__Python__  
[Loading mtx files in Python](#loading-mtx-files-in-python)  
[Loading H5 files in Python](#loading-h5-files-in-python)  

# Basic information regarding data files
Single-cell sequencing data is often stored as sparse matrix objects to cope with the data low density (~10-15% of the entries are non-zero entries). Working with sparse matrices is computationally more efficiency but requires more rigor when it comes down to keeping track of row and column labels, which are stored in separate arrays (**Reminder: Array indexing starts at 0 in Python, 1 in Julia, R and Matlab**). Dataframe objects (R, Python) deal with that aspect but tend to perform slower and can be hard to load in memory when the datasets get larger.  
  
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
<pre><code>#non-zero entry values
M.nzval
Int64[6]
1000
2
1200
400
2
500

#row index for each entry
M.rowval 
Int64[6]
1
3
1
3
2
4

#start of each column. Length is #col +1
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

# R  

## Loading mtx files in R
R base package [`Matrix`](https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/externalFormats.html) lets you load .mtx files directly. The returned object is a [dgTMatrix](https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/dgTMatrix-class.html).You can load a file as shown below:  
<pre><code>> readMM("path_to_matrix/matrix.mtx")  
32738 x 1985 sparse Matrix of class "dgTMatrix"</pre></code> 

Note: Even though Matrix is a base R package, it has to be loaded manually (`Packages` section in R studio or `library("Matrix")` in the R console.)

# Python  

## Loading mtx files in Python
rgewrerg  

## Loading H5 files in Python  
wegfw
