# URT
Fast Unit-Root Tests in C++ with wrappers for R and Python

# Description
URT was originally designed to procure speed while keeping a high level of flexibility for the user when running unit-root tests. The original code is in C++ and uses for all linear algebra computations three of the most widely used C++ libraries in this domain, Armadillo, Blaze and Eigen. The user can switch from one library to another as he wishes and compare performaces. While some are faster than other I wanted to give all of them a chance are they are under active development. They all offer great flexibility as they can be compiled using their own BLAS/LAPACK wrappers or by calling external libraries as Intel MKL and OpenBLAS for example. 
URT can also be used under R and Python. The R version is currenty using Armadillo and developped under Rcpp. The Python version is currently using Blaze and developped under Cython.

# What is inside this package ?
- OLS regression
- Augmented Dickey-Fuller test
- Dickey-Fuller Generalized Least Squares test
- Phillips-Perron test
- Kwiatkowski–Phillips–Schmidt–Shin test

# Design
The C++ template class OLS: to get fast unit-root tests, we need a fast and flexible OLS regression allowing to get the parameters (regressor coefficients) solution of the multiple linear equation y = X.b, as well as their variances to compute their t-statistics.
The C++ template class UniRoot: base class from which all tests will inherit  
