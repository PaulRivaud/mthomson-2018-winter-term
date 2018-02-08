println("Installing packages...")
Pkg.update()
Pkg.add("MatrixMarket");
Pkg.add("JLD")
Pkg.add("MultivariateStats")
Pkg.add("HDF5")
Pkg.add("MAT")
