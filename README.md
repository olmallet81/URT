# URT
Fast Unit-Root Tests in C++ with wrappers for R and Python

# Description
URT is designed to procure speed while keeping a high level of flexibility for the user when running unit-root tests. 
The original code is in C++ and calls for all linear algebra computations three of the most widely used C++ libraries in this domain, Armadillo, Blaze and Eigen. The user can switch from one library to another as he wishes and compare performaces. While some are faster than other I wanted to give all of them a chance are they are under active development. They all offer great flexibility as they can be compiled using their own BLAS/LAPACK wrappers or by calling external libraries as Intel MKL and OpenBLAS for example. 
URT can also be used under R and Python. The R version is currenty using Armadillo and developped under Rcpp. The Python version is currently using Blaze and developped under Cython.
URT contains an OLS regression and four different unit-root tests among the most famous: ADF, DF-GLS, Phillips-Perron and KPSS. ADF and DF-GLS allow for lag length optimization through different methods as information criterion minimization and t-statistic method.

# Why such a project and how can you contribute ?
I have been for a while developping tools for algorithmic trading and it is no secret that unit-root tests are widely used in this domain to decide if a time serie is (weakly) stationary or not and construct on this idea a profitable mean-reversion strategy. Nowadays you often have to look at smaller and smaller time frames to find such trading opportunities and that means on the back-testing side using more and more historical data to test whether the strategy can be profitable or not. I found frustrating that the available libraries under R and Python, which are commonly used in the first steps when building a trading algorithm, were too slow or did not offer enough flexibility. My goal was then to develop a library that could be used under higher level languages to get a first idea on the profitability of a strategy and also when developping a more serious back-testing model on more historical data under a lower level language as C++.
In algorithmic trading we have to find the right amount of data to test for stationarity. If not enough data the back-testing will be faster but the test precision will be smaller and less reliable, if too many data the back-testing will be slower and the test precision will be greater and more reliable.
I have tried to find the correct set up for each linear algebra C++ library (Armadillo, Blaze and Eigen) in order to get the fastest results on a standard sample size of 1000. If anyone can find a faster configuration for one of them he is more than welcome to bring his contribution to this project.

# What is inside this package ?
- OLS regression
- Augmented Dickey-Fuller test
- Dickey-Fuller Generalized Least Squares test
- Phillips-Perron test
- Kwiatkowski–Phillips–Schmidt–Shin test
- Rcpp wrapper for R
- Cython wrapper for Python

# Design
The C++ template class OLS: to get fast unit-root tests, we need a fast and flexible OLS regression allowing to get the parameters (regressor coefficients) solution of the multiple linear equation y = X.b, as well as their variances to compute their t-statistics.
The C++ template class UniRoot: base class from which all tests will inherit  
