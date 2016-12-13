# URT
Fast Unit-Root Tests and OLS regression in C++ with wrappers for R and Python

# Description
URT is designed to procure speed while keeping a high level of flexibility for the user when running unit-root tests. 
The original code is in C++, its core is based on three of the most widely used linear algebra C++ libraries: Armadillo, Blaze and Eigen. The user can switch from one library to another as he wishes and compare performaces. While some are faster than other depending on the Vector and Matrix dimensions I wanted to give them all a chance as they are under active development and future updates could improve their respective performances. They all offer great flexibility as they can be compiled using their own BLAS/LAPACK wrappers or by calling external libraries as for instance Intel MKL and OpenBLAS. 
URT can also be used under R and Python. The R version is currenty using Armadillo and developped under Rcpp using the RcppArmadillo R package. The Python version is currently using Blaze and developped under Cython.
URT contains an OLS regression and four different unit-root tests among the most famous and most used: ADF, DF-GLS, Phillips-Perron and KPSS. ADF and DF-GLS allow for lag length optimization through different methods as information criterion minimization and t-statistic method.

# Why such a project and how can you contribute ?
I have been for a while developping tools for algorithmic trading and it is no secret that unit-root tests are widely used in this domain to decide if a time serie is (weakly) stationary or not and construct on this idea a profitable mean-reversion strategy. Nowadays you often have to look at smaller and smaller time frames to find such trading opportunities and that means on the back-testing side using more and more historical data to test whether the strategy can be profitable on the long term or not. I found frustrating that the available libraries under R and Python, which are commonly used in the first steps when building a trading algorithm, were too slow or did not offer enough flexibility. To that extent I wanted to develop a library that could be used under higher level languages to get a first idea on the profitability of a strategy and also when developping a more serious back-testing model on a larger amount of historical data under a lower level language as C++.
In algorithmic trading we have to find the right sample size to test for stationarity. If we use a too short sample the back-testing will be faster but the test precision will be smaller and the results will be less reliable, on the contrary if we use a too large sample the back-testing will be slower but the test precision will be greater and the results will be more reliable. Hence, when testing for stationarity we have to always keep this tradeoff in mind. Optimal sample size are in general between 100 to 2000, leading to relatively small size arrays. I have then decided not to use parallelism when doing Matrix/Vector operations as it could slow down the code on such small dimensions. Through my benchmarks between linear algebra libraries I am using the sequential versions of Intel MKL and OpenBLAS. Although Armadillo does not allow parallelism yet, Blaze and Eigen do but I made sure to turn off this ability. However, parallelism is used to speed up the lag length optimization in ADF and DF-GLS tests by using OpenMP. 
During my experimentations I have tried to find the correct set up for each linear algebra C++ library (Armadillo, Blaze and Eigen compiled with Intel MKL or OpenBLAS) in order to get the fastest results on a standard sample size of 1000. If anyone can find a faster configuration for one of them he is more than welcome to bring his contribution to this project.

# What is inside this package ?
- OLS regression
- Augmented Dickey-Fuller test
- Dickey-Fuller Generalized Least Squares test
- Phillips-Perron test
- Kwiatkowski–Phillips–Schmidt–Shin test
- Rcpp wrapper for R
- Cython wrapper for Python

# Innovation
Unit-root tests use lags in order to reduce as much as possible auto-correlation in the tested data serie. The test p-value is lag dependent as the critical values will be different depending on the number of lags, several studies have shown this dependency. However very few unit-root tests librairies take this phenomenom into account and return wrong p-values for large number of lags.  

# Design
- C++ template class OLS: 
To get fast unit-root tests, we need a fast and flexible OLS regression allowing to get the parameters (regressor coefficients) solution of the multiple linear equation y = X.b, as well as their variances to compute their t-statistics. This statistics will be used by the unit-root tests to decide whether the serie is (weakly) stationary or not.
The OLS regression is run by simply declaring an OLS object with the parameters:
    - a Vector y containing the dependent variable
    - a Matrix x containing the independent variables (it can include intercept, trend, etc...)
    - a control named stats will compute additional statistics if true as R2, adjusted R2, F statistic and Durbin-Watson statistic

- C++ template class UnitRoot: 
Abstract base class from which all unit-root tests will inherit, it contains all the variables and functions the derived classes ADF, DFGLS, PP and KPSS will need.

  This class has 3 pure virtual functions:
    - statistic() which computes the test statistic
    - pvalue() which calls statistic() and computes the p-value
    - show() with calls pvalue() and output the test results
    
- C++ template class ADF: 
Derived class from UnitRoot, 
