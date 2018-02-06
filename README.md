# Introduction

This Github repository contains data files and code examples to get you started in the Computational Biology class taught by Matt Thomson at Caltech. The goal of this course is to have the students analyze some single-cell sequencing data using programming languages (Julia, R, Python, Matlab, etc.). The content of this README is listed below:  
  
[Basic information regarding data files](#basic-information-regarding-data-files)  
__Julia__  
[Installing Julia (and Juno)](#installing-julia-and-juno)  
[Installing packages in Julia](#installing-packages-in-julia)  
[Running the julia example](#running-the-julia-example)  
__R__  
[Loading mtx files in R](#loading-mtx-files-in-r)

# Basic information regarding data files
Single-cell sequencing data is often stored as sparse matrix objects to cope with the data low density (~10-15% of the entries are non-zero entries). Working with sparse matrices is computationally more efficiency but requires more rigor when it comes down to keeping track of row and column labels, which are stored in separate arrays (**Reminder: Array indexing starts at 0 in Python, 1 in Julia, R and Matlab**). Dataframe objects (R, Python) deal with that aspect but tend to perform slower and can be hard to load in memory when the datasets get larger.  
* Sparse Matrices (MatrixMarket/.mtx format):  
[Info here](https://math.nist.gov/MatrixMarket/formats.html#MMformat)  
Mtx files store MatrixMarket format matrices. The base principle is to store the row and column indices of each non-zero entry in the matrix. The MatrixMarket format contains three distinct parts:  
    - Comment lines that start with `%`
    - Header line: total number of rows, total number of columns, total number of non-zeros entries (space separated)
    - Entries: row index, column index, entry value (space separated)  
The gene and barcode labels are stored in separate files and need to be read separately.  
* Sparse matrices (h5 format):  
[Info here (10X's HDF5 format)](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices)  
HDF5 can store MatrixMarket objects, but can also deconstruct the sparse matrix into multiple vectors: the sparse matrix is a Matrix Market object in RAM when used, but is not stored as one. The different vectors are:
    - data: non-zero entry values (length: number of entries)
    - indices: row indices (length: number of entries)
    - indptr: column index pointers. Index of the start of each column (length: number of columns +1, the last value indicates the end index of the last column +1).
    - barcodes: barcode labels (length: number of columns)
    - gene_names: gene labels (length: number of rows)  
Storing all these vectors under one file enables users to not lose track of matrices respective's labels.  

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

`Note:` Even though Matrix is a base R package, it has to be loaded manually (`Packages` section in R studio or `library("Matrix")` in the R console.)
