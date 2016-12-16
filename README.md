# URT
Fast Unit-Root Tests and OLS regression in C++ with wrappers for R and Python

# Description
URT is a library designed to procure speed while keeping a high level of flexibility for the user when running unit-root tests. 
The original code is in C++, its core is based on three of the most widely used linear algebra C++ libraries: Armadillo, Blaze and Eigen. The user can switch from one library to another as he wishes and compare performaces. While some are faster than other depending on the Vector and Matrix dimensions I wanted to give them all a chance as they are under active development and future updates could improve their respective performances. They all offer great flexibility as they can be compiled using their own BLAS/LAPACK wrappers or by calling external libraries as for instance Intel MKL and OpenBLAS. 
URT can also be used under R and Python. The R version is currenty using Armadillo and developped under Rcpp using the RcppArmadillo R package. The Python version is currently using Blaze and developped under Cython.
URT contains an OLS regression and four different unit-root tests among the most famous and most used: ADF, DF-GLS, Phillips-Perron and KPSS. ADF and DF-GLS allow for lag length optimization through different methods as information criterion minimization and t-statistic method.

# Why such a project and how can you contribute ?
I have been for a while developping tools for algorithmic trading and it is no secret that unit-root tests are widely used in this domain to decide whether a time serie is (weakly) stationary or not and construct on this idea a profitable mean-reversion strategy. Nowadays you often have to look at smaller and smaller time frames to find such trading opportunities and that means on the back-testing side using more and more historical data to test whether the strategy can be profitable on the long term or not. I found frustrating that the available libraries under R and Python, which are commonly used in the first steps when building a trading algorithm, were too slow or did not offer enough flexibility. To that extent I wanted to develop a library that could be used under higher level languages to get a first idea on the profitability of a strategy and also when developping a more serious back-testing model on a larger amount of historical data under a lower level language as C++.
In algorithmic trading we have to find the right sample size to test for stationarity. If we use a too short sample the back-testing will be faster but the test precision will be smaller and the results will be less reliable, on the contrary if we use a too large sample the back-testing will be slower but the test precision will be greater and the results will be more reliable. Hence, when testing for stationarity we have to always keep this tradeoff in mind. Optimal sample size are in general between 100 to 2000, leading to relatively small size arrays. I have then decided not to use parallelism when doing Matrix/Vector operations as it could slow down the code on such small dimensions. Through my benchmarks between linear algebra libraries I am using the sequential versions of Intel MKL and OpenBLAS. Although Armadillo does not allow parallelism yet, Blaze and Eigen do but I made sure to turn off this ability. However, parallelism is used to speed up the lag length optimization in ADF and DF-GLS tests by using OpenMP. 
During my experimentations I have tried to find the correct set up for each linear algebra C++ library (Armadillo, Blaze and Eigen compiled with Intel MKL or OpenBLAS) in order to get the fastest results on a standard sample size of 1000. If anyone can find a faster configuration for one of them he is more than welcome to bring his contribution to this project.
More generally, if anyone has an idea about any kind of modifications that could allow the code to run significantly faster, feel free to propose. 

# What is inside this package ?
- Ordinary Least Squares regression
- Augmented Dickey-Fuller test
- Dickey-Fuller Generalized Least Squares test
- Phillips-Perron test
- Kwiatkowski–Phillips–Schmidt–Shin test
- Lag dependent unit-root test p-values
- Bootstrapped unit-root test p-values 
- Wrapper to expose URT to R 
- Wrapper to expose URT to Python

# Innovation
Unit-root tests use lags in order to reduce auto-correlation as much as possible in the serie being tested. The test p-value is lag dependent as the critical values will be different depending on the number of lags, several studies have shown this dependency and it can easily been proved by Monte-Carlo simulations. However, very few unit-root tests librairies take this phenomenom into account and return wrong p-values for a large number of lags. The method used in this project is the one explained by Cheung and Lai in "Lag Order and Critical Values of the Augmented Dickey-Fuller Test" (1995). This method has been pushed further and adapted to other unit-root tests. 

The methodology is simple, starting from a chosen set of sample sizes and a chosen set of number of lags, it consists in 3 steps:
- step 1: generate a non-stationary random sample of a given size (Wiener process) for the case ADF, DF-GLS and Phillips-Perron tests and a stationary random sample of a given size (Gaussian noise) for the case the KPSS test
- step 2: compute the corresponding test statistic for a given number of lags
- repeat step 1 and 2 many times to get a sample of test statistics from a given couple (sample size, number of lags)
- step 3: sort the statistics sample obtained to get their distribution and record the critical value for each confidence level of your choice
- repeat step 1 to 3 for all combinations (sample size, number of lags) and fit by OLS regression these critical values for each required significance level to the equation proposed by Cheung and Lai:

    ![screenshot from 2016-12-16 17-10-54](https://cloud.githubusercontent.com/assets/20603093/21269345/b6abd474-c3b2-11e6-8247-d43163a11b39.png)
    
    where CR(N,k) is the critical value estimate for a sample of size N and a number of lags k (and for a given significance level), T = N - k being the effective number of observations and e(N,k) the model residuals

In order to increase the precision of the method some terms have been added to each sum while trying to get significant heteroskedasticity consistent t-statistics for the regression coefficients obtained. Both sample sizes and number of lags sets proposed by Cheung and Kai have been expanded. For the most important critical values that is the ones at the confidence levels 1%, 5% and 10% for ADF, DF-GLS and Phillips-Perron tests and at 99%, 95% and 90% for the KPSS test, Monte-Carlo critical values have been computed using a high number of simulations and for reduced set of sizes and lags to compare and improve the estimated critical values precision by modifying the initial set of sizes and lags and by adding or removing some terms to the equation proposed by Cheung and Lai.

The coefficients obtained by OLS regression for each unit-root test and each significance level are reported in the header files in ./URT/include:
- Coeff_adf.hpp for ADF test
- Coeff_dfgls.hpp for DF-GLS test
- Coeff_pp.hpp for Phillips-Perron tests (t-statistic and normalized statistic)
- Coeff_kpss.hpp for KPSS test

NB: each index 0 array contains the asymptotic estimate of the critical value for the corresponding significance level *Tau(0)* and the coefficients of the first term of the equation *Tau(i)*, each index 1 array contains contains the coefficients of the second term of the equation *Phi(j)*.

# Requirement
To use this package you will need at least one of these 3 free linear algebra libraries:
- Armadillo version >= 7.600.1
- Blaze version >= 3.0
- Eigen version >= 3.3.1

It is not necessary to install Intel MKL and/or OpenBLAS (both are now free) for BLAS/LAPACK routines as these 3 libraries have their own wrapper, however I recommend to link to one of them as they will run a lot faster especially Armadillo and Blaze, for Eigen the difference is small. All these libraries will obviously need to be on your path.

NB: if you decide to link to Intel MKL or OpenBLAS, please use their sequential and not parallel version. Intel MKL has the two versions already by default, however OpenBLAS needs to be built from source as sequential with USE_THREAD=0. 
 
# Compilation
URT is not header only to provide an easy way to be exported as a shared library to Rcpp and Cython. Build the shared library using the provided makefile located in ./URT/build. The makefile has been written for Linux and GNU/gcc, it can be easily adapted to run under Windows/OSX and with another compiler, adapt this makefile to your own requirements. 

All linear algebra librairies now use vectorization, in general it should be enabled by default for all of them when compiling with -march=native.

In URT, the user can set the following variables to get the desired shared library:
- USE_OPENMP = 1 to activate parallelism in URT
- USE_BLAZE = 1 to use Blaze library
- USE_EIGEN = 1 to use Eigen library
- USE_MKL = 1 to use Intel MKL library
- USE_BLAS = 1 to use OpenBLAS library

The default configuration when running *make* is no parallelism and Armadillo without external BLAS/LAPACK libraries.

Example: *make USE_OPENMP=1 USE_BLAZE=1 USE_MKL=1* => the shared library libURT.so will be compiled using OpenMP and Blaze with Intel MK.

NB: when compiling with Intel MKL or OpenBLAS, static version of these libraries have been chosen, the shared library obtained will be larger in size but URT will run faster under C++ and the difference will be more important when wrapped for R and Python. You will need to adjust the path of these static libraries in the makefile to your own paths. You are free to rather link to their dynamic version.

# Design

## Introduction
All URT classes and functions are within the namespace *urt*. As URT allows the use of three different linear algebra libraries, convienent typedefs are defined for manipulating arrays, Vector and Matrix defined as below:

- ### with Armadillo
   ```c++
   namespace urt {
      template <typename T>
      using Matrix = arma::Mat<T>;
      template <typename T>
      using Vector = arma::Col<T>;
   }
   ```

- ### with Blaze
   ```c++
   namespace urt {
      template <typename T>
      using Matrix = blaze::DynamicMatrix<T, blaze::columnMajor>;
      template <typename T>
      using CMatrix = blaze::CustomMatrix<T, blaze::unaligned, blaze::unpadded, blaze::columnMajor>;
      template <typename T>
      using Vector = blaze::DynamicVector<T>;
      template <typename T>
      using CVector = blaze::CustomVector<T, blaze::unaligned, blaze::unpadded>;
   }
   ```

- ### with Eigen
   ```c++
   namespace urt {
      template <typename T>
      using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
      template <typename T>
      using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
   }
   ```

## C++ template class OLS: 
To get fast unit-root tests, we need a fast and flexible OLS regression allowing to get the parameters (regressor coefficients) solution of the multiple linear equation y = X.b, as well as their variances to compute the t-statistics. These statistics will be used by the unit-root tests to decide whether the serie is (weakly) stationary or not.

- ### Constructor
    The OLS regression is run by declaring an OLS object using the following constructor:
    ```c++
    OLS<T>::OLS(const Vector<T>& y, const Matrix<T>& x, bool stats = false)
    ```
    with:
    - *y* = vector of the dependent variable
    - *x* = matrix of the independent variables (it can include intercept, constant trend, etc...)
    - *stats* if turned to *true*, additional statistics will be computed as R-squared, adjusted R-squared, F statistic and Durbin-Watson statistic
     
- ### Member variables (public)
   - *param* = regressors coefficients
   - *t_stat* = regressors t-statistics
   - *resid* = regression residuals
   - *var* = regressors variances
   - *MSE* = mean of squares for error
   - *SSR* = sum of squares residuals
   - *R2* = R-squared
   - *adj_R2* = adjusted R-squared
   - *F_stat* = F-statistic
   - *DW_stat* = Durbin-Watson statistic
   - *IC* = information criterion
   - *nobs* = number of observations
   - *nreg* = number of regressors
   - *ndef* = number of degrees of freedom
   - *lags* = number of lags
    
NB: *IC* and *lags* are for the case when OLS is called by UnitRoot for lag length optimization by information criterion minimization.
    
- ### Member functions (public)
    - *get_stats()* taking same *x* and *y* as arguments as the constructor, computes the additional OLS regression statistics and detects the presence of an intercept term
    - *show()* outputs the results
    
- ### Additionnal tools (not members of OLS)
URT provides 3 functions allowing to add quickly constant terms to a Matrix:
    - *add_intercept()* inserts a column of ones into a Matrix as shown in the example above
    - *add_constant_trend()* inserts a column (1,2,3,...)
    - *add_quadratic_trend()* inserts a column (1,4,9,...)
    
Code example using Armadillo:

```c++
#include "./URT/include/URT.hpp"

int main()
{
    int nrows = 1000;
    int ncols = 10;

    // generating random arrays
    urt::Vector<double> y = arma::randn<urt::Vector<double>>(nrows);
    urt::Matrix<double> x = arma::randn<urt::Matrix<double>>(nrows, ncols);

    urt::add_intercept(x);

    urt::OLS<double> fit(y, x, true);

    fit.show();
}
    
```
NB: I made the choice not to copy the Vector and Matrix arguments when declaring an object OLS for performance reasons, when the Matrix becomes large it can quickly lead to a significative difference in term of performance. Also, if *stats* has not been set to "true" the function "get_stats()" will not be called and the intercept will not be detected in the output.

## C++ template class UnitRoot: 
Abstract base class from which all unit-root tests will inherit, it contains all the variables and functions the derived classes ADF, DFGLS, PP and KPSS will need.

  This class has 3 pure virtual functions:
    - "statistic()" computes the test statistic
    - "pvalue()" calls "statistic()" and computes the p-value
    - "show()" calls "pvalue()" and outputs the test results
    
- C++ template class ADF: 
Derived class from UnitRoot, this class has 2 constructors:
    - 1st constructor to compute the test for a given lag
    - 2nd constructor to compute the test with lag length optimization
The constructors accept the following arguments:
    - Vector "data" containing the data on which the test will be processed
    - the number of "lags" (1st constructor only)
    - the optimization "method" (2nd constructor only)
    - the type of "trend" (optional, default to constant or intercept term)
    - control "regression" (optional, default to false) to indicate if the additional statistics of the OLS regression must be computed
Once an object ADF declare, one of the overriden pure virtual functions from the base class must be called to get the test results. The user can switch from a test to another by simply modifying the arguments above.
The user can switch from a test for a given lag to a test with lag length optimization by using "method"
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
