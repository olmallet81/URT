# URT
Fast Unit-Root Tests and OLS regression in C++ with wrappers for R and Python

# Description
URT is a library designed to procure speed while keeping a high level of flexibility for the user.

The core code is in C++ and based on three of the most widely used C++ linear algebra libraries: Armadillo, Blaze and Eigen. The user can switch from one library to another as he wishes and compare performaces. While some are faster than other depending on array dimensions all of them have been given a chance as they are under active development and future updates could improve their respective performances. They all offer great flexibility as they can be compiled using their own BLAS/LAPACK wrappers or by calling external libraries as for instance Intel MKL and OpenBLAS. 

URT can also be used under R and Python. The R version is currenty using Armadillo and developped under Rcpp using the RcppArmadillo R package. The Python version is currently using Blaze and developped under Cython.
URT contains an OLS regression and four of the most famous and most used unit-root tests: ADF, DF-GLS, Phillips-Perron and KPSS. ADF and DF-GLS allow for lag length optimization through different methods as information criterion minimization and t-statistic method. Test p-values can be computed via an extension of the method proposed by Cheung and Lai back in 1995 or by bootstrap. 

# Why such a project and how can you contribute ?
I have been developping tools for algorithmic trading for a while and it is no secret that unit-root tests are widely used in this domain to decide whether a time serie is (weakly) stationary or not and construct on this idea a profitable mean-reversion strategy. Nowadays you often have to look at smaller and smaller time frames to find such trading opportunities and that means on the back-testing side using more and more historical data to test whether the strategy can be profitable on the long term or not. I found frustrating that the available libraries under R and Python, which are commonly used in the first steps when building a trading algorithm, were too slow or did not offer enough flexibility. To that extent I wanted to develop a library that could be used under higher level languages to get a first idea on the profitability of a strategy and also when developping a more serious back-tester on a larger amount of historical data under a lower level language as C++. 

In algorithmic trading we have to find the right sample size to test for stationarity. If we use a too short sample of historica data on a rolling window the back-testing will be faster but the test precision will be smaller and the results will be less reliable, on the contrary if we use a too large sample the back-testing will be slower but the test precision will be greater and the results will be more reliable. Hence, when testing for stationarity we have to always keep this tradeoff in mind. Optimal sample size are in general between 100 to 2000, leading to relatively small size arrays. I have then decided not to use parallelism when calling Matrix/Vector operations as it would not bring any speed improvement and on the contrary would slow down the code when applied on such small dimensions. Although Armadillo does not allow for parallelism yet, Blaze and Eigen do, I made sure to turn off this ability. However, parallelism is used to speed up the lag length optimization by information criterion minimization in ADF and DF-GLS tests by using OpenMP. All of these libraries are now using vectorization (from SSE to AVX), activating this feature can improve greatly the general performance.

During my experimentations I have tried to find the correct set up for each C++ linear algebra library (Armadillo, Blaze and Eigen compiled with Intel MKL or OpenBLAS) in order to get the fastest results on a sample size of 1000. If anyone can find a faster configuration for one of them or more generally, if anyone has something to propose that could make the C++ code faster or even the Cython and Rcpp wrappers, he is more than welcome to bring his contribution to this project.

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
- repeat step 1 to 3 for all combinations of sample size and number of lags and fit by OLS regression these critical values for each required significance level to the equation proposed by Cheung and Lai:

    ![screenshot from 2016-12-16 17-10-54](https://cloud.githubusercontent.com/assets/20603093/21269345/b6abd474-c3b2-11e6-8247-d43163a11b39.png)
    
    where CR(N,k) is the critical value estimate for a sample of size N and a number of lags k (and for a given significance level), T = N - k being the effective number of observations and Epsilon(N,k) the model residuals

In order to increase the precision of the method some terms have been added by increasing degree of the first sum and/or the second sum, while trying to get significant heteroskedasticity consistent t-statistics for the regression coefficients obtained. Both sample sizes and number of lags sets proposed by Cheung and Kai have been expanded. For the most important critical values that is the ones at the confidence levels 1%, 5% and 10% for ADF, DF-GLS and Phillips-Perron tests and at 99%, 95% and 90% for the KPSS test, Monte-Carlo critical values have been computed using a high number of simulations and for reduced set of sizes and lags to compare and improve the estimated critical values precision by modifying the initial set of sizes and lags and by adding or removing some terms to the equation proposed by Cheung and Lai.

The coefficients obtained by OLS regression for each unit-root test and each significance level are reported in the header files in ./URT/include:
- Coeff_adf.hpp for ADF test
- Coeff_dfgls.hpp for DF-GLS test
- Coeff_pp.hpp for Phillips-Perron tests (t-statistic and normalized statistic)
- Coeff_kpss.hpp for KPSS test

NB: each index 0 array contains the asymptotic estimate of the critical value for the corresponding significance level *Tau(0)* and the coefficients of the first term of the equation *Tau(i)*, each index 1 array contains contains the coefficients of the second term of the equation *Phi(j)*.

# Requirements
To use this package you will need at least one of these 3 free linear algebra libraries:
- Armadillo version >= 7.600.1
- Blaze version >= 3.0
- Eigen version >= 3.3.1

It is not necessary to install Intel MKL and/or OpenBLAS (both are now free) for BLAS/LAPACK routines as these 3 libraries have their own wrapper, however I recommend to link to one of them as they will run a lot faster especially Armadillo and Blaze, for Eigen the difference is small. All these libraries will obviously need to be on your path.

NB: if you decide to link to Intel MKL or OpenBLAS, please use their sequential and not parallel version. Intel MKL has the two versions already by default, however OpenBLAS needs to be built from source as sequential with USE_THREAD=0. 
 
# Compilation
URT is not header only to provide an easy way to be exported as a shared library to Rcpp and Cython. Build the shared library using the provided makefile located in ./URT/build. The makefile has been written for Linux and GNU/gcc, it can be easily adapted to run under Windows/OSX and with another compiler, adapt this makefile to your own requirements. 

All linear algebra librairies now use vectorization, in general it should be enabled by default for all of them when compiling with -march=native flag.

In URT, the user can set the following variables to get the desired shared library:
- USE_OPENMP = 1 to activate parallelism in URT
- USE_BLAZE = 1 to use Blaze library
- USE_EIGEN = 1 to use Eigen library
- USE_MKL = 1 to use Intel MKL library
- USE_BLAS = 1 to use OpenBLAS library

The default configuration when running *make* is no parallelism and Armadillo without external BLAS/LAPACK libraries.

Example: *make USE_OPENMP=1 USE_BLAZE=1 USE_MKL=1* => the shared library libURT.so will be compiled using OpenMP and Blaze with Intel MKL.

NB: when compiling with Intel MKL or OpenBLAS, static version of these libraries have been chosen, the shared library obtained will be larger in size but URT will run faster under C++ and the difference will be more important when wrapped for R and Python. You will need to adjust the path of these static libraries in the makefile to your own paths. You are free to rather link to their dynamic version.

# Design

## Introduction
All URT classes and functions are within the namespace *urt*. As URT allows the use of three different linear algebra libraries, convienent typedefs Vector and Matrix have been defined for manipulating arrays:

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
Declared in ./URT/include/OLS.hpp, defined in ./URT/src/OLS.cpp.

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
   - *DW_stat* = Durbin-Watson test statistic
   - *IC* = information criterion
   - *nobs* = number of observations
   - *nreg* = number of regressors
   - *ndef* = number of degrees of freedom
   - *lags* = number of lags
    
NB: *IC* and *lags* are for the case when OLS is called by UnitRoot for lag length optimization by information criterion minimization.
    
- ### Member functions (public)
    - *get_stats()* takes same *x* and *y* as arguments as the constructor, computes the additional OLS regression statistics and detects the presence of an intercept term
    - *show()* outputs the results
    
- ### Additionnal tools (not members of OLS)
URT provides 3 functions located in ./URT/include/Tools.hpp, to add quickly constant terms to a Matrix object:
    - *add_intercept()* inserts a column of ones into a Matrix as shown in the example above
    - *add_constant_trend()* inserts a column (1,2,3,...)
    - *add_quadratic_trend()* inserts a column (1,4,9,...)
    
- ### Code example using Armadillo:

    ```c++
    #include "./URT/include/URT.hpp"

    int main()
    {
       int nrows = 1000;
       int ncols = 10;

       // generating random arrays
       urt::Vector<double> y = arma::randn<urt::Vector<double>>(nrows);
       urt::Matrix<double> x = arma::randn<urt::Matrix<double>>(nrows, ncols);

       // adding intercept to matrix of independent variables
       urt::add_intercept(x);

       // running OLS regression
       urt::OLS<double> fit(y, x, true);

       // outputting regression results
       fit.show();
    }   
    ```
NB: the choice has been made not to copy Vector and Matrix, arguments of OLS constructor for performance reasons. Indeed, when the Matrix becomes large it can quickly lead to a significative difference in term of performance. Also, if *stats* has not been set to "true" the function "get_stats()" will not be called and the intercept will not be detected in the output.

## C++ template class UnitRoot
Declared in ./URT/include/UnitRoot.hpp, defined in ./URT/src/UnitRoot.cpp.

Abstract base class from which all unit-root tests will inherit, it contains all the variables and functions the derived classes ADF, DFGLS, PP and KPSS will need.

- ### Member variables (public)
    - *lags_type* = default number of lags, long (Schwert l12-rule) or short (Schwert l4-rule)
    - *method* = method for lag length optimization, for ADF and DF-GLS tests only, possible choices are:
        - *AIC* = Aikaike Information Criterion
        - *BIC* = Bayesian Information Criterion
        - *HQC* = Hannah-Quinn Criterion
        - *MAIC* = Modified Aikaike Information Criterion
        - *MBIC* = Modified Bayesian Information Criterion
        - *MHQC* = Modified Hannah-Quinn Criterion
        - *T-STAT* = optimal lag length selection using (absolute) threshold
    - *test_type* for Phillips-Perron test only, possibles choices are:
        - *tau* = t-statistic test (default value)
        - *rho* = normalized statistic test
    - *trend* = regression trend, possible choices are:
        - *nc* = no constant (for ADF and PP tests only)
        - *c* = constant (for all tests, default value)
        - *ct* = constant trend (for all tests)
        - *ctt* = quadratic trend (for ADF test only)
    - *level* = statistic (absolute) threshold for optimal lag length selection, for T-STAT method only, set to *1.64* by default
    - *lags* = model number of lags
    - *max_lags* = maximum number of lags, for models with lag length optimization, initialized to Schwert l12-rule by default
    - *niter* = number of iterations when computing test p-value by bootstrap, initialized to 1000 by default
    - *bootstrap* = if set to *true*, test p-value will be computed by bootstrap, set to *false* by default
    - *regression* = if set to *true*, OLS regression results will be outputted when running *show()* method, set to *false* by default
    
- ### Member functions (public)
    - *get_stat()* to return the test statistic
    - *get_pval()* to return the test p-value
    - *get_ols()* to return the OLS regression results (returns object of type urt::OLS)
    - *get_trends()* to return possible trends for the current test (returns std::vector)
    
    This class has also 3 pure virtual functions:
    - *statistic()* computes the test statistic
    - *pvalue()* calls *statistic()* and computes the p-value
    - *show()* calls *pvalue()* and outputs the test results
    
 NB: UnitRoot template class is designed in a way that test statistic and test p-value will not be re-calculated if all parameters remain identical. For example, if user calls *pvalue()* method and after *show()*, *pvalue()* will not be run again unless the user has modified at least one parameter of the current test.  

## C++ template class ADF
Declared in ./URT/include/ADF.hpp, defined in ./URT/src/ADF.cpp.

Derived template class from UnitRoot, this class has 2 constructors:

- ### Constructor for computing ADF test for a given number of lags

    ```c++
    ADF(const Vector<T>& data, int lags, const std::string& trend = "c", bool regression = false)
    ```

- ### Constructor for computing ADF test with lag length optimization

    ```c++
    ADF(const Vector<T>& data, const std::string& method, const std::string& trend = "c", bool regression = false)
    ```
    
- ### Code example using Armadillo:

    ```c++
    #include "./URT/include/URT.hpp"

    int main()
    {
       int nobs = 1000;

       // generating random non-stationary data
       urt::Vector<double> data = arma::cumsum(arma::randn<urt::Vector<double>>(nobs));

       // initializing ADF test with 10 lags and constant trend
       urt::ADF<double> test(data, 10, "ct");

       // outputting test results
       test.show();
    
       // switching to test with lag length optimization and p-value computation by bootstrap with 10000 iterations
       test.method = "AIC";
       test.bootstrap = true;
       test.niter = 10000;
    
       // outputting test results
       test.show();  
    }
    ```
   
## C++ template class DFGLS
Declared in ./URT/include/DFGLS.hpp, defined in ./URT/src/DFGLS.cpp.

Derived template class from UnitRoot, this class has 2 constructors:

- ### Constructor for computing DF-GLS test for a given number of lags

    ```c++
    DFGLS(const Vector<T>& data, int lags, const std::string& trend = "c", bool regression = false)
    ```

- ### Constructor for computing DF-GLS test with lag length optimization

    ```c++
    DFGLS(const Vector<T>& data, const std::string& method, const std::string& trend = "c", bool regression = false)
    ```   
    
- ### Code example using Armadillo:

    ```c++
    #include "./URT/include/URT.hpp"

    int main()
    {
       int nobs = 1000;

       // generating random non-stationary data
       urt::Vector<double> data = arma::cumsum(arma::randn<urt::Vector<double>>(nobs));

       // initializing DFGLS test with lag length optimization using BIC and constant term
       urt::DFGLS<double> test(data, "AIC");

       // outputting test results
       test.show();
    
       // switching to test for 3 lags and no constant
       test.lags = 3;
       test.method = "nc";
    
       // outputting test results
       test.show();
    
       // switching to test with optimal lag length selection and maximum lags of 10;
       test.method = "T-STAT";
       test.max_lags = 10;
    
       // outputting test results
       test.show();
    }
    ```

## C++ template class PP
Declared in ./URT/include/PP.hpp, defined in ./URT/src/PP.cpp.

Derived template class from UnitRoot, this class has 2 constructors:

- ### Constructor for computing Phillips-Perron test for a given number of lags

    ```c++
    PP(const Vector<T>& data, int lags, const std::string& trend = "c", const std::string& test_type = "tau", bool regression = false)
    ```

- ### Constructor for computing Phillips-Perron test with a default number of lags (long or short)

    ```c++
    PP(const Vector<T>& data, const std::string& lags_type, const std::string& trend = "c", const std::string& test_type = "tau", bool regression = false)
    ``` 
- ### Code example using Blaze:

    ```c++
    #include "./URT/include/URT.hpp"

    int main()
    {
       int nobs = 1000;

       // generating random stationary data
       urt::Vector<float> data = blaze::Rand<urt::Vector<float>>(nobs);
       
       // initializing Phillips-Perron normalized test with lags of type long and constant term
       urt::PP<float> test(data, "long", "c", "rho");

       // outputting test results
       test.show();
    
       // switching to t-statistic test with 
       test.test_type = "tau";
    
       // outputting test results
       test.show();  
    }
    ```
    
## C++ template class KPSS
Declared in ./URT/include/KPSS.hpp, defined in ./URT/src/KPSS.cpp.

Derived template class from UnitRoot, this class has 2 constructors:

- ### Constructor for computing KPSS test for a given number of lags

    ```c++
    KPSS(const Vector<T>& data, int lags, const std::string& trend = "c")
    ```

- ### Constructor for computing KPSS test with a default number of lags (long or short)

    ```c++
    KPSS(const Vector<T>& data, const std::string lags_type, const std::string& trend = "c")
    ```
    
- ### Code example using Eigen:

    ```c++
    #include "./URT/include/URT.hpp"

    int main()
    {
       int nobs = 1000;

       // generating random stationary data
       urt::Vector<float> data = urt::Vector<float>::Random(nobs);

       // initializing KPSS test with lags of type short and constant trend
       urt::KPSS<double> test(data, "short", "ct");

       // outputting test results
       test.show();
    
       // switching to test with 5 lags and constant term
       test.lags = 5;
       test.trend = "c";
    
       // outputting test results
       test.show();  
    }
    ```  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
