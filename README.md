# Introduction

This Github repository contains data files and code examples to get you started in the Computational Biology class taught by Matt Thomson at Caltech. The goal of this course is to have the students analyze some single-cell sequencing data using programming languages (Julia, R, Python, Matlab, etc.). The content of this README is listed below:  
[Basic information regarding data files](#basic-information-regarding-data-files)  
[Installing Julia (and Juno)](#installing-julia-and-juno)  
[Installing packages in Julia](#installing-packages-in-julia)  
[Running the julia example](#running-the-julia-example)

# Basic information regarding data files
Single-cell sequencing data is often stored as sparse matrix objects to cope with the data low density (~10-15% of the entries are non-zero entries). Working with sparse matrices is computationally more efficiency but requires more rigor when it comes down to keeping track of row and column labels, which are stored in separate arrays. Dataframe objects (R, Python) deal with that aspect but tend to perform slower and can be hard to load in memory when the datasets get larger.

# Installing Julia (and Juno)
Visit the JuliaComputing web page https://juliacomputing.com/products/juliapro.html to navigate your way to the download section. You may be asked to enter an email address to have access to a free download.

* Linux:  
Select "JuliaPro-0.6.2.1 - MKL (for Linux)". Once the download is complete, just uncompress the archive and add the bin folder to the PATH variable in the .bashrc file or from the terminal. The command should look similar to the following:

<pre><code>export PATH=$PATH:/home/user/Downloads/julia-235245113/bin</code></pre>

* Mac OSX:  
Select "JuliaPro-0.6.2.1-MKL (for Mac)". Once the download is complete, uncompress the archive and go through the installation process.

Those installations will install Julia (command line version) and should install Juno as well, a Julia IDE based on Atom. If you have trouble installing Juno, we recommend to look at this page: https://github.com/JunoLab/uber-juno/blob/master/setup.md.

# Installing packages in Julia
Packages are installed thanks to Julia commands (https://docs.julialang.org/en/stable/stdlib/pkg/). The most common way to install and use packages is:
<pre><code>Pkg.add("packageName")  
using packageName</pre></code>  

# Running the Julia example
Prior to running the example featured in that repository, run the install_packages.jl script using the following:
<pre><code>julia install_packages.jl</pre></code>  
Look up the [Installing packages](#installing-packages) section to learn more about package management.

# Loading mtx files in R

