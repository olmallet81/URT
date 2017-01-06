# URT
Fast Unit-Root Tests and OLS regression in C++ with wrappers for R and Python

# Description
URT is a library designed to procure speed while keeping a high level of flexibility for the user.

The core code is in C++ and based on three of the most widely used C++ linear algebra libraries: Armadillo, Blaze and Eigen. The user can switch from one library to another and compare performaces. While some are faster than other depending on array dimensions all of them have been given a chance as they are under active development and future updates might improve their respective performances. They can all be compiled by calling external libraries as for instance Intel MKL and OpenBLAS or by using their own BLAS/LAPACK wrappers.

URT can also be used under R and Python. The R version is currenty using Armadillo and developped under Rcpp using the RcppArmadillo R package. The Python version is currently using Blaze and developped under Cython.
URT contains an OLS regression and four of the most famous and most used unit-root tests: ADF, DF-GLS, Phillips-Perron and KPSS. ADF and DF-GLS allow for lag length optimization through different methods as information criterion minimization and t-statistic method. Test p-values can be computed via an extension of the method proposed by Cheung and Lai back in 1995 or by bootstrap. 

# Why such a project and how can you contribute ?
I have been developping tools for algorithmic trading for a while and it is no secret that unit-root tests are widely used in this domain to decide whether a time serie is (weakly) stationary or not and construct on this idea a profitable mean-reversion strategy. Nowadays you often have to look at smaller and smaller time frames to find such trading opportunities and that means on the back-testing side using more and more historical data to test whether the strategy can be profitable on the long term or not. I found frustrating that the available libraries under R and Python, which are commonly used in the first steps when building a trading algorithm, were too slow or did not offer enough flexibility. To that extent I wanted to develop a library that could be used under higher level languages to get a first idea on the profitability of a strategy and also when developping a more serious back-tester on a larger amount of historical data under a lower level language as C++. 

In algorithmic trading we have to find the right sample size to test for stationarity. If we use a too short sample of historica data on a rolling window the back-testing will be faster but the test precision will be smaller and the results will be less reliable, on the contrary if we use a too large sample the back-testing will be slower but the test precision will be greater and the results will be more reliable. Hence, when testing for stationarity we have to always keep this tradeoff in mind. Optimal sample size are in general between 100 to 2000, leading to relatively small size arrays. I have then decided not to use parallelism when calling Matrix/Vector operations as it would not bring any speed improvement and on the contrary would slow down the code when applied on such small dimensions. Although Armadillo does not allow for parallelism yet, Blaze and Eigen do, I made sure to turn off this ability. However, parallelism is used to speed up the lag length optimization by information criterion minimization in ADF and DF-GLS tests by using OpenMP. All of these libraries are now using vectorization (from SSE to AVX), activating this feature can improve greatly the general performance.

During my experimentations I have tried to find the correct set up for each C++ linear algebra library (Armadillo, Blaze and Eigen compiled with Intel MKL or OpenBLAS) in order to get the fastest results on a sample size of 1000. If anyone can find a faster configuration for one of them, or more generally, if anyone has anything to propose that could make the C++ code or the Cython and Rcpp wrappers faster, he is more than welcome to bring his contribution to this project.

# What is inside this package ?
- Ordinary Least Squares regression
- Augmented Dickey-Fuller test
- Dickey-Fuller Generalized Least Squares test
- Phillips-Perron test
- Kwiatkowski–Phillips–Schmidt–Shin test
- Lag dependent unit-root test p-values

