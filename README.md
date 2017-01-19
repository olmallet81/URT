# URT
Fast Unit Root Tests and OLS regression in C++ with wrappers for R and Python.

# Description
URT is a library designed to procure speed while keeping a high level of flexibility for the user when testing for a unit root in a time serie.

URT core code is in C++ and based on three of the most widely used C++ linear algebra libraries: [Armadillo](http://arma.sourceforge.net/), [Blaze](https://bitbucket.org/blaze-lib/blaze) and [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). The user can switch from one library to another and compare performances. While some are faster than other depending on array dimensions all of them have been given a chance as they are under active development and future updates might improve their respective performances. They can all be compiled by linking to external libraries for high-speed [BLAS](http://www.netlib.org/blas/)/[LAPACK](http://www.netlib.org/lapack/) replacements for better performance such as [Intel MKL](https://software.intel.com/en-us/intel-mkl) and [OpenBLAS](http://www.openblas.net/) or by using their own BLAS/LAPACK routines. 

URT can also be used under R and Python. The wrapper for R called RcppURT is currenty using Armadillo and developped under [Rcpp](http://www.rcpp.org/)  and [R6](https://cran.r-project.org/web/packages/R6/vignettes/Introduction.html) using the R package [RcppArmadillo](http://dirk.eddelbuettel.com/code/rcpp.armadillo.html). The wrapper for Python called CyURT is currently using Blaze and developped under [Cython](http://cython.readthedocs.io/en/latest/src/userguide/wrapping_CPlusPlus.html) for C++.

URT contains an [Ordinary Least Squares regression](#c-template-class-ols) (OLS) and four of the most famous unit root tests: the [Augmented Dickey-Fuller test](#c-template-class-adf) (ADF), the [Dickey-Fuller Generalized Least Squares test](#c-template-class-dfgls) (DF-GLS), the [Phillips-Perron test](#c-template-class-pp) and the [Kwiatkowski–Phillips–Schmidt–Shin test](#c-template-class-kpss) (KPSS). ADF and DF-GLS allow for lag length optimization through different methods such as information criterion minimization and t-statistic. Test p-values can be computed via an extension of the method proposed by Cheung and Lai back in 1995 or by bootstrap. 

URT is single-threaded for most of unit root tests but goes parallel for the most time consuming ones by using [OpenMP](http://www.openmp.org/). The unit root tests concerned are ADF and DF-GLS with lag length optimization by information criterion minimization.

URT has been written under Linux but can be easily adapted to run Windows or OSX as all the libraries used in this project exist on these platforms.

# Why such a project and how can you contribute ?
I have been developing algorithmic trading tools for a while and it is no secret that unit root tests are widely used in this domain to decide whether a time serie is (weakly) stationary or not and construct on this idea a profitable mean-reversion strategy. Nowadays you often have to look at smaller and smaller time frames as minute data to find such trading opportunities and that means on the back-testing side using more and more historical data to test whether the strategy can be profitable on the long term or not. I found frustrating that the available libraries under R and Python, interpreted languages commonly used in the first steps of building a trading algorithm, were too slow or did not offer enough flexibility. To that extent I wanted to develop a library that could be used under higher level languages to get a first idea on the profitability of a strategy and also when developping a more serious back-tester on a larger amount of historical data under a lower level language such as C++. 

In algorithmic trading we have to find the right sample size to test for stationarity. If we use a too short sample of historical data on a rolling window the back-testing will be faster but the test precision will be smaller and the results will be less reliable, on the contrary if we use a too large sample the back-testing will be slower but the test precision will be greater and the results will be more reliable. Hence, when testing for stationarity we have to always keep this tradeoff in mind. Sample sizes used are usually between 100 to 5000, leading to relatively small size arrays. I have then decided not to use parallelism for matrix and vector operations as it would not bring any speed improvement and on the contrary would slow down the code when applied on such small dimensions. Although Armadillo does not allow for parallelism yet, Blaze and Eigen do, I made sure to turn off this ability. However, parallelism is used to speed up the lag length optimization by information criterion minimization in ADF and DF-GLS tests by enabling OpenMP. All of these libraries are now using vectorization (from SSE to AVX), activating this feature greatly improves the general performance.

During my experimentations I have tried to find the correct set up for each C++ linear algebra library (Armadillo, Blaze and Eigen compiled with either Intel MKL or OpenBLAS) in order to get the fastest results on a standard sample size of 1000. If anyone can find a faster configuration for any of them, or more generally, if anyone has anything to propose that could make the C++ code or the Cython and Rcpp wrappers faster, he is more than welcome to bring his contribution to this project.

# What is inside this repository ?
- [Ordinary Least Squares regression](#c-template-class-ols) 
- [Augmented Dickey-Fuller test](#c-template-class-adf)
- [Dickey-Fuller Generalized Least Squares test](#c-template-class-dfgls)
- [Phillips-Perron test](#c-template-class-pp) 
- [Kwiatkowski–Phillips–Schmidt–Shin test](#c-template-class-kpss)
- [Lag dependent unit root tests critical values and p-values](#innovation)
- [CyURT: wrapper to run URT under Python](#cyurt-urt-for-python) 
- [RcppURT: wrapper to run URT under R](#rcppurt-urt-for-r)
- [Benchmarks](#benchmarks)

# Innovation
Unit root tests use lags in order to reduce auto-correlation as much as possible in the time serie being tested. The test p-value is lag dependent as the critical values will be different depending on the number of lags, several studies have shown this dependency and it can easily been proved by Monte-Carlo simulations. However, very few unit root tests librairies take this phenomenom into account and return wrong p-values for a large number of lags. The method used in this project is the one explained by Cheung and Lai in "Lag Order and Critical Values of the Augmented Dickey-Fuller Test" (1995). This method has been pushed further and adapted to other unit root tests. 

The method is simple, starting from a chosen set of sample sizes and a chosen set of number of lags, it consists in 3 steps:
- step 1: generate a non-stationary random sample (Wiener process) of a given size for ADF, DF-GLS and Phillips-Perron tests and a stationary random sample (Gaussian noise) of a given size for the KPSS test
- step 2: compute the corresponding test statistic for a given number of lags
- repeat step 1 and 2 many times to get the test statistics for a given pair sample size and number of lags
- step 3: sort the statistics obtained to get their distribution and record the critical value for all required significance levels 
- repeat step 1 to 3 for all possible pairs of sample size and number of lags and fit by OLS regression these critical values for all required significance levels to the equation proposed by Cheung and Lai:

    ![screenshot from 2016-12-16 17-10-54](https://cloud.githubusercontent.com/assets/20603093/21269345/b6abd474-c3b2-11e6-8247-d43163a11b39.png)
    
    where CR(N,k) is the critical value estimate for a sample of size N and a number of lags k (and for a given significance level), T = N - k being the effective number of observations and Epsilon(N,k) the model residuals

In order to increase the precision of the method some terms have been added going further than degree 2 for the first sum and/or the second sum, while trying to get significant heteroskedasticity consistent t-statistics for the regression coefficients obtained. Both sample sizes and number of lags sets proposed by Cheung and Kai have been expanded. For the most important critical values that is the ones at the significance levels 1%, 5% and 10% for ADF, DF-GLS and Phillips-Perron tests and 99%, 95% and 90% for the KPSS test, Monte-Carlo critical values have been computed using a high number of simulations and for reduced sets of sizes and lags to compare and improve the estimated critical values precision by modifying the initial set of sizes and lags and by adding or removing some terms to the original equation proposed by Cheung and Lai.

The coefficients obtained by OLS regression for each unit root test and each significance level are reported in the header files in ./URT/include:
- Coeff_adf.hpp for the ADF test
- Coeff_dfgls.hpp for the DF-GLS test
- Coeff_pp.hpp for the Phillips-Perron tests (t-statistic and normalized statistic)
- Coeff_kpss.hpp for the KPSS test

NB: arrays indexed by 0 contain the asymptotic estimate of the critical value for the corresponding significance level *Tau(0)* and the coefficients of the first term of the equation *Tau(i)*, arrays indexed by 1 contain the coefficients of the second term of the equation *Phi(j)*.

# Requirements
To use the C++ version of URT you will need to install at least one of these three free C++ linear algebra libraries:
- Armadillo version >= 7.600.1
- Blaze version >= 3.0
- Eigen version >= 3.3.1

You will also need to have the C++ libraries Boost (version >= 1.61.0) already installed.

Blaze library requires at least LAPACK to run whereas Armadillo and Eigen have their own internal BLAS/LAPACK routines, however for better performance I recommend installing one of the following high-speed BLAS/LAPACK replacement libraries:
- Intel MKL version >= 2017.0.098
- OpenBLAS version >= 0.2.19

These libraries will need to be on your C++ compiler path. If you decide to link to Intel MKL or OpenBLAS, please use their sequential and not their multi-threaded version. When installing Intel MKL you will get the two versions, however OpenBLAS needs to be built from source as sequential using USE_THREAD=0.

To use the Python wrapper CyURT you will need to install the C++ linear algebra library Blaze and:
- Python version >= 2.7
- Numpy version >= 1.11.3 (preferably built from source with Intel MKL and OpenMP enabled)
- Pandas version >= 0.19.2
- Cython version >= 0.24.1

NB: Pandas is not essential, it will be used in the Python example only.

To use the R wrapper RcppURT you will need to install the C++ linear algebra library Armadillo and:
- R version >= 3.3.2 (preferably built from source with Intel MKL and OpenMP enabled)
- Rcpp version >= 0.12.8
- RcppArmadillo version >= 0.7.600.1.0
- R6 version >= 2.2.0

NB: Some of these tools might work with older versions, I only reported the versions I used to build this project. 

# Compilation
URT is not header only to provide a direct way to be exported as a shared library to Rcpp to be exposed to R and to Cython to be exposed to Python. Build the shared library using the provided makefile located in ./URT/build. The makefile has been written for Linux and GNU/gcc, it can be easily modified to run under Windows or OSX and with another compiler, you are free to adapt this makefile to your own requirements. 

The three C++ linear algebra librairies Armadillo, Blaze and Eigen use vectorization, it will be enabled for all of them when compiling with the flag -march=native. Parallelism is used to speed up ADF and DF-GLS tests only when searching for an optimal lag length by information criterion minimization. Although these linear algebra libraries can use parallelism for matrix and vector operations by calling the multi-threaded version of Intel MKL or OpenBLAS for example, I decided to use their single threaded versions as the usual sample sizes used in unit root tests are relatively small, they are indeed rarely larger than a few thousand. It would not bring any speed improvement, on the contrary it could slow down the code. Moreover, it would induce nested parallism when running ADF test or DF-GLS test with lag length optimization.

To build the shared library, the user can set the following variables:
- USE_OPENMP = 1 to activate parallelism in URT
- USE_BLAZE = 1 to use Blaze library
- USE_EIGEN = 1 to use Eigen library
- USE_MKL = 1 to use Intel MKL library
- USE_BLAS = 1 to use OpenBLAS library

The default configuration when running *make* is OpenMP disabled and Armadillo using its internal BLAS/LAPACK routines. 

Example: *make USE_OPENMP=1 USE_BLAZE=1 USE_BLAS=1* => URT shared library will be built with OpenMP enabled and with the C++ linear algebra library Blaze using OpenBLAS for BLAS/LAPACK routines.

NB: Armadillo does not need any external library for BLAS/LAPACK routines, however it needs to be linked to its shared library. Blaze can run with internal BLAS routines but needs to be linked to an external LAPACK library. Eigen can run without calling any external library.

## Example

- step 1: build URT shared library libURT.so using the makefile under ./URT/build with:
```
$ cd URT/build
$ make USE_OPENMP=1 USE_ARMA=1 USE_MKL=1
$ cd ../../
```

- step 2: export the shared library location with:
```
$ export LD_LIBRARY_PATH=/path/to/URT/lib:$LD_LIBRARY_PATH
```

- step 3: compile example1.cpp in ./URT/examples with:
```
$ g++ -std=c++14 -O3 -march=native -DUSE_ARMA -o run -L./URT/lib ./URT/examples/example1.cpp -lURT
```

- step 4: run executable with: 
```
$ ./run
``` 

You can repeat step 3 and 4 to compile other examples in ./URT/examples.

Instead of exporting URT shared library path, under Linux you can rather add the following lin to th .bashrc file:
```
export LD_LIBRARY_PATH=/path/to/URT/lib:$LD_LIBRARY_PATH 
```

# Design

## Introduction
All URT classes and functions are within the namespace *urt*. As URT allows the use of three different linear algebra libraries, convienent typedefs Vector and Matrix have been defined for manipulating arrays:

- ### with Armadillo
    ```C++
    namespace urt {
       template <typename T>
       using Matrix = arma::Mat<T>;
       template <typename T>
       using Vector = arma::Col<T>;
    }
    ```

- ### with Blaze
    ```C++
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
    ```C++
    namespace urt {
       template <typename T>
       using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
       template <typename T>
       using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    }
    ```
    
NB: It is important to mention the differences of behaviour between these libraries when assigning a matrix to a vector: Armadillo will convert this vector into a matrix, Blaze will return a compilation error and Eigen will assign to this vector the first column of the matrix. 

## C++ template class *OLS*
Declared in ./URT/include/OLS.hpp, defined in ./URT/src/OLS.cpp.

To get fast unit root tests, we need a fast and flexible OLS regression allowing to get the parameters (regressor coefficients) solution of the multiple linear equation y = X.b, as well as their variances to compute the t-statistics. These statistics will be used by the unit-root tests to decide whether the serie is (weakly) stationary or not.

- ### Constructor
    The OLS regression is run by instantiating an *OLS* object using the following constructor:
    ```C++
    OLS<T>::OLS(const Vector<T>& y, const Matrix<T>& x, bool stats = false)
    ```
    with:
    - *y* = vector of the dependent variable
    - *x* = matrix of the independent variables (it can include intercept, constant trend, etc...)
    - *stats* if turned to *true*, additional statistics will be computed such as R-squared, adjusted R-squared, F test statistic and Durbin-Watson test statistic
     
- ### Member variables (public)
    - *param* = regressors coefficients
    - *t_stat* = regressors t-statistics
    - *resid* = regression residuals
    - *var* = regressors variances
    - *MSE* = mean of squares for error
    - *SSR* = sum of squares residuals
    - *R2* = R-squared
    - *adj_R2* = adjusted R-squared
    - *F_stat* = F test statistic
    - *DW_stat* = Durbin-Watson test statistic
    - *IC* = information criterion
    - *nobs* = number of observations
    - *nreg* = number of regressors
    - *ndef* = number of degrees of freedom
    - *lags* = number of lags
    
    NB: *IC* and *lags* are for the case when *OLS* is used by a unit root test class for lag length optimization by information criterion minimization.
    
- ### Member functions (public)
    - *get_stats()* takes same *x* and *y* as arguments as the constructor, detects the presence of an intercept term and computes the additional regression statistics
    - *show()* outputs the results
    
- ### Additionnal tools (not members of *OLS*)
URT provides three functions located in ./URT/include/Tools.hpp, to add quickly constant terms to a *Matrix* object:
    - *add_intercept()* inserts a column of ones at the end of a *Matrix* object as shown in the example below
    - *add_constant_trend()* inserts a constant trend at the end of a *Matrix* object (1,2,3,...)
    - *add_quadratic_trend()* inserts a quadratic trend at the end of a *Matrix* object (1,4,9,...)
    
    In the same header are defined some functions generating random data for testing URT classes:
    - *gaussian_noise()* will generate a *Vector* or *Matrix* object of normally distributed random data (stationary process)
    - *wiener_process()* will generate a a *Vector* or *Matrix* object of Wiener processes (non-stationary process)
    
    The header ./URT/include/CsvManager.hpp contains functions for CSV files management:
    - *WriteToCSV()* will write a *Matrix* or *Vector* object to a CSV file
    - *ReadFromCSV()* will read a CSV file and store the data into a *Matrix* or *Vector* object

    Data can then be exported from C++ to R or Python to validate the results returned by URT wrappers.
    
- ### Exceptions
*OLS* class constructor will throw an exception if *x* and *y* do not have the same number of rows or if at least one of them is empty.
    
- ### Code example:

    ```C++
    // ./URT/examples/example1.cpp
    #include "../include/URT.hpp"

    int main()
    {
       int nrows = 1000;
       int ncols = 10;

       // generating random arrays
       urt::Vector<double> y = urt::wiener_process<double>(nrows);
       urt::Matrix<double> x = urt::wiener_process<double>(nrows, ncols);

       // adding intercept to matrix of independent variables
       urt::add_intercept(x);
       
       // writting data to CSV files
       urt::WriteToCSV("./URT/data/y.csv", y);
       urt::WriteToCSV("./URT/data/x.csv", x);

       // running OLS regression
       urt::OLS<double> fit(y, x, true);

       // outputting regression results
       fit.show();
    }   
    ```
    
- ### Ouput:
    ```
    OLS Regression Results
    ======================

    Coefficients
    --------------------------------------------------
                Estimate  Std Error   t value  P(>|t|)
    Intercept    4.51898     0.4968     9.097    0.000
           x1    0.10240     0.0304     3.366    0.001
           x2   -0.14263     0.0205    -6.959    0.000
           x3   -0.04654     0.0232    -2.003    0.046
           x4    0.00560     0.0367     0.153    0.879
           x5   -0.27260     0.0306    -8.895    0.000
           x6    0.41597     0.0270     15.42    0.000
           x7   -0.02540     0.0247    -1.027    0.305
           x8    0.18289     0.0217     8.433    0.000
           x9    0.07821     0.0303     2.581    0.010
          x10    0.15704     0.0243     6.466    0.000

    Residuals
    ----------------------------------------
         Min      1Q  Median      3Q     Max
      -12.04   -2.19   -0.06    2.40    9.45

    Dimensions
    ----------------------------------------
    Number of observations              1000
    Number of degrees of freedom         989
    Number of variables                   10

    Analysis
    ----------------------------------------
    Residual mean                    1.3e-14
    Residual standard error            3.594
    Multiple R-squared               0.79603
    Adjusted R-squared               0.79397
    F-statistic (p-value)             385.98 (0.00)
    Durbin-Watson statistic            0.106
    ```
    
NB: the choice has been made not to copy Vector and Matrix, arguments of *OLS* class constructor for performance reasons. Indeed, when the Matrix becomes large it can quickly lead to a significative difference in term of performance. Also, if *stats* has not been set to *true* the function *get_stats()* will not be called and the intercept will not be detected in the output.

## C++ template class *UnitRoot*
Declared in ./URT/include/UnitRoot.hpp, defined in ./URT/src/UnitRoot.cpp.

Abstract base class from which all unit root tests will inherit, it contains all the variables and functions the derived classes [*ADF*](#c-template-class-adf), [*DFGLS*](#c-template-class-dfgls), [*PP*](#c-template-class-pp) and [*KPSS*](#c-template-class-kpss) will need.

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
    - *level* = statistic (absolute) threshold for optimal lag length selection, for *T-STAT* method only, set to *1.64* by default
    - *lags* = model number of lags
    - *max_lags* = maximum number of lags, for models with lag length optimization, initialized to Schwert l12-rule by default
    - *niter* = number of iterations when computing test p-value by bootstrap, initialized to 1000 by default
    - *bootstrap* = if set to *true*, test p-value will be computed by bootstrap, set to *false* by default
    - *regression* = if set to *true*, OLS regression results will be outputted when running *show()* method, set to *false* by default
    
- ### Member functions (public)
    - *get_stat()* returns the test statistic
    - *get_pval()* returns the test p-value
    - *get_ols()* returns the OLS regression results
    - *get_trends()* returns the possible trends for the current test
    
    This class has also three pure virtual functions:
    - *statistic()* computes the test statistic
    - *pvalue()* calls *statistic()* and computes the p-value
    - *show()* calls *pvalue()* and outputs the test results
    
    NB: *UnitRoot* template class is designed in a way that test statistic and test p-value will not be re-calculated if all parameters remain identical. For example, if user calls *pvalue()* method and after *show()*, *pvalue()* will not be run again unless the user has modified at least one parameter of the current test.
 
- ### Exceptions
For ADF and DF-GLS tests an exception will be thrown when running the member function *statistic()* or *pvalue()* or *show()* if the serie being tested does not have enough elements for the required number of lags. The number of additional elements to be added will be outputted.

## C++ template class *ADF*
Declared in ./URT/include/ADF.hpp, defined in ./URT/src/ADF.cpp.

Derived template class from [*UnitRoot*](#c-template-class-unitroot), this class has two constructors:

- ### Constructor for computing ADF test for a given number of lags

    ```C++
    ADF(const Vector<T>& data, int lags, const std::string& trend = "c", bool regression = false)
    ```
    
- ### Constructor for computing ADF test with lag length optimization

    ```C++
    ADF(const Vector<T>& data, const std::string& method, const std::string& trend = "c", bool regression = false)
    ```
    
- ### Code example:

    ```C++
    // ./URT/examples/example2.cpp
    #include "../include/URT.hpp"

    int main()
    {
       int nobs = 1000;

       // generating non-stationary random data
       urt::Vector<double> data = urt::wiener_process<double>(nobs);

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
    
- ### Ouput:

    - First ADF test:
    ```
    Augmented Dickey-Fuller Test Results
    ====================================
    Statistic                     -2.949
    P-value                        0.152
    Lags                              10
    Trend                 constant trend
    ------------------------------------

    Test Hypothesis
    ------------------------------------
    H0: The process contains a unit root
    H1: The process is weakly stationary

    Critical Values
    ---------------
     1%      -3.956
     5%      -3.405
    10%      -3.121

    Test Conclusion
    ---------------
    We cannot reject H0
    ```
    - Second ADF test:
    ```
    Augmented Dickey-Fuller Test Results
    ====================================
    Statistic                     -2.761
    P-value                        0.215 (*)
    Optimal Lags                       0
    Criterion                        AIC
    Trend                 constant trend
    ------------------------------------

    Test Hypothesis
    ------------------------------------
    H0: The process contains a unit root
    H1: The process is weakly stationary

    Critical Values (*)
    ---------------
     1%      -4.013
     5%      -3.402
    10%      -3.116

    (*) computed by bootstrap

    Test Conclusion
    ---------------
    We cannot reject H0
    ```
    
## C++ template class *DFGLS*
Declared in ./URT/include/DFGLS.hpp, defined in ./URT/src/DFGLS.cpp.

Derived template class from [*UnitRoot*](#c-template-class-unitroot), this class has two constructors:

- ### Constructor for computing Dickey-Fuller GLS test for a given number of lags
    ```C++
    DFGLS(const Vector<T>& data, int lags, const std::string& trend = "c", bool regression = false)
    ```
    
- ### Constructor for computing Dickey-Fuller GLS test with lag length optimization
    ```C++
    DFGLS(const Vector<T>& data, const std::string& method, const std::string& trend = "c", bool regression = false)
    ```   
    
- ### Code example:
    ```C++
    // ./URT/examples/example3.cpp
    #include "../include/URT.hpp"

    int main()
    {
       int nobs = 1000;

       // generating non-stationary random data
       urt::Vector<double> data = urt::wiener_process<double>(nobs);

       // initializing DFGLS test with lag length optimization using BIC and constant term
       urt::DFGLS<double> test(data, "BIC");

       // outputting test results
       test.show();
    
       // switching to test with 10 lags and constant trend
       test.trend = "ct";
       test.lags = 10;
    
       // outputting test results
       test.show();
    }
    ``` 
    
- ### Ouput:

    - First Dickey-Fuller GLS test:
    ```
        Dickey-Fuller GLS Test Results
    ====================================
    Statistic                     -1.655
    P-value                        0.098
    Optimal Lags                       0
    Criterion                        BIC
    Trend                       constant
    ------------------------------------

    Test Hypothesis
    ------------------------------------
    H0: The process contains a unit root
    H1: The process is weakly stationary

    Critical Values
    ---------------
     1%      -2.586
     5%      -1.963
    10%      -1.642

    Test Conclusion
    ---------------
    We can reject H0 at the 10% significance level
    ```
    
    - Second Dickey-Fuller GLS test:
    ``` 
    Dickey-Fuller GLS Test Results
    ====================================
    Statistic                     -1.914
    P-value                        0.360
    Lags                              10
    Trend                 constant trend
    ------------------------------------

    Test Hypothesis
    ------------------------------------
    H0: The process contains a unit root
    H1: The process is weakly stationary

    Critical Values
    ---------------
     1%      -3.409
     5%      -2.851
    10%      -2.565

    Test Conclusion
    ---------------
    We cannot reject H0
    ```

## C++ template class *PP*
Declared in ./URT/include/PP.hpp, defined in ./URT/src/PP.cpp.

Derived template class from [*UnitRoot*](#c-template-class-unitroot), this class has two constructors:

- ### Constructor for computing Phillips-Perron test for a given number of lags
    ```C++
    PP(const Vector<T>& data, int lags, const std::string& trend = "c", const std::string& test_type = "tau", bool regression = false)
    ```
    
- ### Constructor for computing Phillips-Perron test with a default number of lags (long or short)
    ```C++
    PP(const Vector<T>& data, const std::string& lags_type, const std::string& trend = "c", const std::string& test_type = "tau", bool regression = false)
    ``` 
    
- ### Code example:
    ```C++
    // ./URT/examples/example4.cpp
    #include "./URT/include/URT.hpp"

    int main()
    {
       int nobs = 1000;

       // generating non-stationary random data
       urt::Vector<float> data = urt::wiener_process<float>(nobs);
       
       // initializing Phillips-Perron normalized test with lags of type long and constant term
       urt::PP<float> test(data, "long", "c", "rho");

       // outputting test results
       test.show();
    
       // switching to t-statistic test 
       test.test_type = "tau";
    
       // outputting test results
       test.show();  
    }
    ```
    
- ### Ouput:

    - First Phillips-Perron test:
    ```
    Phillips-Perron Test Results (Z-rho)
    ====================================
    Statistic                     -6.077
    P-value                        0.371
    Lags                              21
    Trend                       constant
    ------------------------------------

    Test Hypothesis
    ------------------------------------
    H0: The process contains a unit root
    H1: The process is weakly stationary

    Critical Values
    ---------------
     1%     -20.548
     5%     -14.058
    10%     -11.225

    Test Conclusion
    ---------------
    We cannot reject H0
    ```
    
    - Second Phillips-Perron test:
    ```
    Phillips-Perron Test Results (Z-tau)
    ====================================
    Statistic                     -1.573
    P-value                        0.497
    Lags                              21
    Trend                       constant
    ------------------------------------

    Test Hypothesis
    ------------------------------------
    H0: The process contains a unit root
    H1: The process is weakly stationary

    Critical Values
    ---------------
     1%      -3.440
     5%      -2.867
    10%      -2.570

    Test Conclusion
    ---------------
    We cannot reject H0
    ```
   
## C++ template class *KPSS*
Declared in ./URT/include/KPSS.hpp, defined in ./URT/src/KPSS.cpp.

Derived template class from [*UnitRoot*](#c-template-class-unitroot), this class has two constructors:

- ### Constructor for computing KPSS test for a given number of lags
    ```C++
    KPSS(const Vector<T>& data, int lags, const std::string& trend = "c")
    ```
   
- ### Constructor for computing KPSS test with a default number of lags (long or short)
    ```C++
    KPSS(const Vector<T>& data, const std::string lags_type, const std::string& trend = "c")
    ```
    
- ### Code example:
    ```C++
    // ./URT/examples/example5.cpp
    #include "../include/URT.hpp"

    int main()
    {
       int nobs = 1000;

       // generating stationary random data
       urt::Vector<float> data = urt::gaussian_noise<float>(nobs);

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
    
- ### Ouput:

    - First KPSS test:
    ```
        KPSS Test Results
    ====================================
    Statistic                      0.029
    P-value                        0.900
    Lags                               7
    Trend                 constant trend
    ------------------------------------

    Test Hypothesis
    ------------------------------------
    H0: The process is weakly stationary
    H1: The process contains a unit root

    Critical Values
    ---------------
     1%       0.213
     5%       0.147
    10%       0.119

    Test Conclusion
    ---------------
    We cannot reject H0
    ```
    
    - Second KPSS test:
    ```
    KPSS Test Results
    ====================================
    Statistic                      0.092
    P-value                        0.648
    Lags                               5
    Trend                       constant
    ------------------------------------

    Test Hypothesis
    ------------------------------------
    H0: The process is weakly stationary
    H1: The process contains a unit root

    Critical Values
    ---------------
     1%       0.732
     5%       0.459
    10%       0.347

    Test Conclusion
    ---------------
    We cannot reject H0
    ```

## CyURT: URT for Python 
URT can be called from Python. The Cython wrapper has been written with the C++ linear algebra library Blaze. This wrapper is located in ./URT/Python.

Before testing CyURT under Python make sure to have built the URT shared library under ./URT/build with Blaze using for example the command:
```
$ make USE_BLAZE=1
```
If you want URT to be multi-threaded and linked to Intel MKL for better performance, use instead the command:
```
$ make USE_OPENMP=1 USE_BLAZE=1 USE_MKL=1
```
To compile the Cython code and build the shared library CyURT.so that will be imported from Python, run setup.py under ./URT/python with the command: 
```
$ python setup.py build_ext --inplace
```
Before running the Python code export the URT C++ shared library path (if you did not add this path already to the .bashrc file) with the command:
```
$ export LD_LIBRARY_PATH=/path/to/URT/lib:$LD_LIBRARY_PATH
```
Before running the example below make sure to have run ./URT/examples/example1.cpp (to write random data to CSV files, the same that were used for the C++ examples).

You are now ready to run ./URT/python/example.py:

```Python
import numpy as np
import pandas as pd

import CyURT as urt
 
if __name__ == "__main__":

    y = pd.read_csv('../data/y.csv', sep=',', header=None)
    x = pd.read_csv('../data/x.csv', sep=',', header=None)

    yd = np.asarray(y).reshape(y.size)
    yf = yd.astype(np.float32)
    xd = np.asarray(x, order='F')
    xf = xd.astype(np.float32)

    # running OLS regression as in ./examples/example1.cpp using double precision type
    fit = urt.OLS_d(yd, xd, True)
    fit.show()
    
    # running OLS regression as in ./examples/example1.cpp using single precision type
    fit = urt.OLS_f(yf, xf, True)
    fit.show()

    # running first ADF test as in ./examples/example2.cpp using double precision type
    test = urt.ADF_d(yd, lags=10, trend='ct')
    test.show()

    # running second ADF test as in ./examples/example2.cpp using double precision type
    test.method = 'AIC'
    test.bootstrap = True
    test.niter = 10000
    test.show()
```

The Cython wrapper behaves the same way than under C++, the only difference being when the user wants single precision instead of double precision, he will have to convert Python data, double precision by default to single precision as shown in the example above with *yf* and *xf* and the URT class name followed by *_f* (*OLS_f* instead of *OLS_d*).

The following changes have been made from the C++ original code to the Cython wrapper:
- *get_stat()* method in C++ has become the class property *stat*
- *get_pval()* method in C++ has become the class property *pval*
- *get_ols()* method in C++ has become the class property *ols*
- *get_trends()* method in C++ has become the class property *trends*

Important: all URT classes accept Numpy arrays only as arguments, *OLS_d* and *OLS_f* classes need a 1-dimension array for the dependent variable vector and a 2-dimension array for the matrix of independent variables. All other classes (unit root tests) need a 1-dimension array. Blaze matrices have been set to be column-major so Numpy arrays need to be Fortran style.

NB: When passing a Numpy array from Python to the Cython wrapper, no copy will be made, the same memory will be re-used for better performance. The same will happen when returning an array from C++ to Python, Numpy under Cython will wrap this array without making any copy.
    
## RcppURT: URT for R  
URT can be called from R. The Rcpp wrapper has been written with the C++ linear algebra library Armadillo and the R package RcppArmadillo. This package is located in ./URT/R.

RcppURT contains two different wrappers for URT C++ classes:
    
- The first wrapper as shown in ./URT/R/example1.R has been written using R6 classes, and behaves the same way than under C++. As for the Python wrapper, the URT C++ class name followed by *_d* is for double precision type and followed by *_f* for single precision type (example: *OLS_d* and *OLS_f* for C++ class *OLS*). 

- The second wrapper as shown in ./URT/R/example2.R has been written as Rcpp functions of URT C++ classes to avoid adding interpreted code as for the first wrapper and to improve the performance in the case the user does not need classes flexibility. These functions names are:

    - *OLSreg_d()* and *OLSreg_f()* for OLS regression
    - *ADFtest_d()* and *ADFtest_f()* for the Augmented Dickey-Fuller test
    - *DFGLStest_d()* and *DFGLStest_f()* for the Dickey-Fuller test with GLS detrending
    - *PPtest_d()* and *PPtest_f()* for the Phillips-Perron test
    - *KPSStest_d()* and *KPSStest_f()* for the Kwiatkowski–Phillips–Schmidt–Shin test

To get URT working under R, I recommend building an R package from URT C++ source files. The R package called RcppURT is already prepared under ./URT/R. All URT headers have been placed into the include directory and all source files into the src directory. Adjust the Makevars file in the src directory whether you want to compile Armadillo with external link to Intel MKL or to OpenBLAS (or any other BLAS/LAPACK library of your choice as long as Armadillo can accept it). 

To build the RcppURT package run under ./URT/R the following command: 
```
$ R CMD build RcppURT
```

Once the package is built, install it with root rights with the following command: 
```
$ R CMD INSTALL --no-lock RcppURT_1.0.tar.gz
```
    
Before running the examples below make sure to have run ./URT/examples/example1.cpp (to write random data to CSV files, the same that were used for the C++ examples).

You are now ready to run ./URT/R/example1.R:

```R
suppressMessages(library(RcppURT))

run <- function()
{
  x = as.matrix(read.table("../data/x.csv", sep=","))
  y = as.matrix(read.table("../data/y.csv", sep=","))

  # running OLS regression as in ./examples/example1.cpp using double precision type
  fit = OLS_d$new(y, x, stats=TRUE)
  fit$show()

  # running OLS regression as in ./examples/example1.cpp using single precision type
  fit = OLS_f$new(y, x, stats=TRUE)
  fit$show()

  # running first ADF test as in ./examples/example2.cpp using double precision type
  test = ADF_d$new(y, lags=10, trend='ct')
  test$show()

  # running second ADF test as in ./examples/example2.cpp using double precision type
  test$method = 'AIC'
  test$bootstrap = TRUE
  test$niter = 10000
  test$show()
}
```
and ./URT/R/example2.R:
    
```R
suppressMessages(library(RcppURT))

run <- function()
{
  x = as.matrix(read.table("../data/x.csv", sep=","))
  y = as.matrix(read.table("../data/y.csv", sep=","))

  # running OLS regression as in ./examples/example1.cpp using double precision type
  fit = OLSreg_d(y, x, stats=TRUE, output=TRUE)

  # running OLS regression as in ./examples/example1.cpp using single precision type
  fit = OLSreg_f(y, x, stats=TRUE, output=TRUE)

  # running first ADF test as in ./examples/example2.cpp using double precision type
  test = ADFtest_d(y, lags=10, trend='ct', output=TRUE)

  # running second ADF test as in ./examples/example2.cpp using double precision type
  test = ADFtest_d(y, method='AIC', trend='ct', output=TRUE, bootstrap=TRUE, niter=10000)
}
```

The choice has been made not to use Rcpp modules to wrap URT C++ classes as the performance was very poor. For unit root test classes we could also have created a base R6 class wrapping C++ *UnitRoot* base class from which all unit root test R6 classes would have inherited but the performance would have been worse than directly including all C++ *UnitRoot* base class variables and methods required into the R6 class wrappers.

The following changes have been made from the C++ original code to the R6 wrapper:
- *get_stat()* method in C++ has become the class property *stat* 
- *get_pval()* method in C++ has become the class property *pval*
- *get_ols()* method in C++ has become the class property *ols*
- *get_trends()* method in C++ has become the class property *trends*

NB: When passing an array from R to the Rcpp wrapper, no copy will be made for double precision type, the same memory will be re-used for better performance. However, as R accepts double precision type only, a copy will be made when passing an array from R to the Rcpp wrapper for single precision type. When returning an array from C++ to R, this array will be first wrapped under Rcpp and then copied in R.
    
# Benchmarks

## C++

In this section we are going to compare URT performance using alternatively Armadillo, Blaze and Eigen, each one of these linear algebra libraries using alternatively Intel MKL, OpenBLAS and their internal BLAS/LAPACK routines (at the exception of Blaze which must be linked at least to an external LAPACK library).
Once the nine URT shared libraries have been built we are ready to proceed by running for each one of them ./URT/benchmark/benchmark.cpp:

```C++
#include "../include/URT.hpp"

#ifndef USE_ARMA
  #include <armadillo>
#endif

// define USE_FLOAT when compiling to switch to single precision
#ifdef USE_FLOAT
  using T = float;
#else
  using T = double;
#endif

int main()
{
   int niter = 0;

   arma::wall_clock timer;

   std::vector<int> sizes = {100,150,200,250,300,350,400,450,500,1000,1500,2000,2500,3000,3500,4000,4500,5000};

   std::cout << std::fixed << std::setprecision(1);

   for (int i = 0; i < sizes.size(); ++i) {
   
      urt::Vector<T> data = urt::wiener_process<T>(sizes[i]);

      (sizes[i] < 1000) ? niter = 10000 : niter = 1000;

      timer.tic();
      for (int k = 0; k < niter; ++k) {
         urt::ADF<T> test(data, "AIC");
         test.statistic();
      }

      auto duration = timer.toc();

      std::cout << std::setw(8) << sizes[i];
      std::cout << std::setw(8) << duration << "\n";
   }
}
```
The benchmark is run on small sample sizes from 100 to 500 (10000 iterations) and large sample sizes from 1000 to 5000 (1000 iterations) for double and single precision types. We choose to compare the linear algebra libraries performances on a time consuming unit root test such as ADF with constant term and lag length optimization by AIC criterion minimization. 

NB: The URT libraries have been built single-threaded for the purpose of this benchmark meaning that the research for the optimal number of lags in the ADF test has been done using one thread only. The machine on which this benchmark has been run is equipped with a processor Intel Core i5-3210M @ 2.50GHz and 6GB of RAM. The Intel Turbo Boost has been turned off to keep the processor running at a constant frequency through the different simulations.

The following graphs show the results obtained.

- ### Results on small sample sizes
![graphs1](https://cloud.githubusercontent.com/assets/20603093/21899692/10af315a-d8e9-11e6-85aa-c9df2df4c2a7.png)

- ### Results on large sample sizes
![graphs2](https://cloud.githubusercontent.com/assets/20603093/21899694/15542800-d8e9-11e6-801b-009ef5279135.png)

- ### Observations
Armadillo without the use of OpenBLAS or Intel MKL wrappers obtains poor performance. However, when the wrappers for OpenBLAS or Intel MKL are enabled within Armadillo, the performance is virtually the same as for the other libraries. For small sample sizes Blaze appears to be the fastest for both double and single precision types and more precisely when enabling Intel MKL. For large sample sizes the performances of the three libraries are quite similar excepted for Armadillo alone that tends to be much slower. We can note the good performance of Eigen in that case with or without external BLAS/LAPACK wrappers.

## Cython wrapper

In this section we are going to compare the performance of CyURT with the original URT in C++ code using the linear algebra library Blaze by running ./URT/Python/benchmark.py:

```Python
import numpy as np
import CyURT as urt
from timeit import default_timer as timer 

if __name__ == "__main__":

    sizes = [100,150,200,250,300,350,400,450,500,1000,1500,2000,2500,3000,3500,4000,4500,5000]

    for i in range(len(sizes)):

        # generating Wiener process
        data = np.cumsum(np.random.normal(size=sizes[i]))
        # uncomment this line and comment the one above to switch to single precision
        #data = np.cumsum(np.random.normal(size=sizes[i])).astype(np.float32)

        if sizes[i] < 1000: niter = 10000
        else: niter = 1000
        
        start = timer()
        for k in range(niter):
            test = urt.ADF_d(data, method='AIC')
            # uncomment this line and comment the one above to switch to single precision
            #test = urt.ADF_f(data, method='AIC')
            test.statistic()
        end = timer()

        print '{:8d}'.format(sizes[i]), '{:8.1f}'.format(end - start)
```
- ### Results on small and large sample sizes
![graphs3](https://cloud.githubusercontent.com/assets/20603093/21984017/d9dc79fc-dbeb-11e6-8632-d4a97ac83a3f.png)

- ### Observations
Although a bit slower than the C++ version of URT for small sample sizes, the Python wrapper performance is almost the same for large sample size and even slightly faster as the sample size increase.

- ### Comparing CyURT to ARCH
The Python package [ARCH](https://pypi.python.org/pypi/arch/3.0) (version 4.0) contains some unit root tests, the same benchmark than above has been run with ARCH package using the same ADF test with constant term and lag length optimization by AIC minimization. CyURT has been built with URT using Blaze library linked to Intel MKL and OpenMP enabled. For a fair comparison I made sure that Numpy was also calling Intel MKL libraries with OpenMP enabled. The table below presents the results obtained (in seconds), the ratio column corresponding to ARCH performance by CyURT performance.  

   ![tab1](https://cloud.githubusercontent.com/assets/20603093/21994520/06dabd8e-dc18-11e6-956a-ec74b10e3b6f.png)
   
- ### Comparing CyURT to STATSMODELS
The same comparison than above has been done with the Python package [STATSMODELS](http://statsmodels.sourceforge.net/devel/generated/statsmodels.tsa.stattools.adfuller.html) (version 0.6.1) that contains some unit root tests. The ratio column corresponding this time to STATSMODELS performance (SM column in the table) by CyURT performance.
ARCH and STATSMODELS performances being very similar (probably due to the fact that they both use the same OLS regression function from STATSMODELS package), the performance factor is almost the same as above.

   ![tab2](https://cloud.githubusercontent.com/assets/20603093/21994747/3bf57fa8-dc19-11e6-823d-192352f90ff1.png)

## Rcpp wrapper

In this section we are going to compare the performance of RccpURT, R6 classes and Rcpp functions with the original URT in C++ code using the linear algebra library Armadillo by running ./URT/R/benchmark.R:

```R
suppressMessages(library(RcppURT))

run <- function()
{
  sizes = c(100,150,200,250,300,350,400,450,500,1000,1500,2000,2500,3000,3500,4000,4500,5000)

  for (i in 1:length(sizes)) {

    # generating Wiener process
    data = cumsum(rnorm(n=sizes[i]))

    if (sizes[i] < 1000) niter = 10000
    else niter = 1000
        
    # with R6 classes
    start1 = Sys.time()
    for (k in 1:niter) {
        test = ADF_d$new(data, method='AIC')
        # uncomment this line and comment the one above to switch to single precision
        #test = ADF_f$new(data, method='AIC')
        test$statistic()
    }
    end1 = Sys.time()

    # with Rcpp functions
    start2 = Sys.time()
    for (k in 1:niter) {
        test = ADFtest_d(data, method='AIC')
        # uncomment this line and comment the one above to switch to single precision
        #test = ADFtest_f(data, method='AIC')
    }
    end2 = Sys.time()

    cat(sprintf("%8d", sizes[i]))
    cat(sprintf("%8.1f", end1 - start1))
    cat(sprintf("%8.1f\n", end2 - start2))
  }
}
```

- ### Results on small and large sample sizes
![graphs4](https://cloud.githubusercontent.com/assets/20603093/21982608/142dcdce-dbe5-11e6-9c2a-26fb6c4e8a75.png)

- ### Observations
We can see that for small sample sizes R6 classes wrapper performance is pretty poor due to the extensive use of interpreted code, however Rcpp functions wrapper performance is close to the performance of the original C++ code. For large sample size this difference tends to disappear. We can also notice that for double precision the performance of the Rcpp functions wrapper is closer to the performance of the C++ version of URT and (even beating it for large sample size) than it is for single precision. This is due to the fact that for double precision no copy is made when passing an array from R to C++, whereas for single precision a copy is made to convert the R array in double precision by default to a single precision C++ array. We can see this phenomenom very clearly for small sample sizes: for double precision Rcpp functions wrapper and C++ code tend to converge as the size increase whereas they tend to diverge for single precision as the size of the array to be copied increases. 

- ### Comparing RcppURT to URCA
The R package [URCA](https://cran.r-project.org/web/packages/urca/index.html) (version 1.3-0) contains some unit root tests, the same benchmark than above has been run with URCA package using the same ADF test with constant term and lag length optimization by AIC minimization. RcppURT has been built with URT using Armadillo library linked to Intel MKL and OpenMP enabled. For a fair comparison I made sure that R and URCA package were also built with Intel MKL libraries with OpenMP enabled. The Rcpp functions have been used for the benchmark and not the R6 classes as URCA is made of functions too. The table below presents the results obtained (in seconds), the ratio column corresponding to URCA performance by RcppURT performance.  

   ![tab3](https://cloud.githubusercontent.com/assets/20603093/21994748/3c0bb7b4-dc19-11e6-9e93-f9c53dd0f80a.png)

# References

This project has been realised using the following articles:

- Denis Kwiatkowski, Peter C.B. Phillips, Peter Schmidt and Yongcheol Shin, "Testing the null hypothesis of stationarity against the alternative of a unit root" (October 1991)
- Yin-Wong Cheun and Kon S. Lai and Tuan Tran, "Finite-Sample Critical Values of the KPSS test" (December 1994)
- Yin-Wong Cheun and Kon S. Lai, "Lag Order and Critical Values of the Augmented Dickey-Fuller Test" (Journal of Business & Economic Statistics, July 1995, Vol. 13, No. 3)
- Graham Elliott, Thomas J. Rothenberg and James H. Stock, "Efficient Tests for an Autoregressive Unit Root" (Econometrica, Vol. 64, No. 4 (July 1996), 813-836)
- Serena Ng and Pierre Perron, "Useful Modifications to some Unit Root Tests with Dependent Errors and their Local Asymptotic Properties" (Review of Economic Studies (1996) 63, 435-463)
- Steven Cook, "Finite-sample critical values of the Augmented Dickey-Fuller statistic: a note on lag order" (Economic Issues, Vol.6, Part 2. September 2001)
- Serena Ng and Pierre Perron, "Lag Length Selection and the Construction of Unit Root Tests with Good Size and Power" (Econometrica, Vol. 69, No. 6 (Nov., 2001), pp. 1519-1554)
- Joon Y. Park, "Bootstrap Unit Root Tests" (October 2002)
- Andrew F. Hayes and Li Cai, "Using heteroskedasticity-consistent stanard error estimators in OLS regression: An introduction and software implementation" (Behavior Research Methods 2007, 39 (4), 709-722)
- Franz C. Palm, Stephan Smeekes, Jean-Pierre Urbain, "Bootstrap Unit Root Tests: Comparison and Extensions" (Journal of Time Series Analysis 29 (2), 371–401 (2008))
- Stephan Smeekes, "Detrending Bootstrap Unit Root Tests" (Econometric Reviews 32 (8), 869-891 (2013))
