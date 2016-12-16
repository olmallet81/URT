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

## C++ template class OLS: 
To get fast unit-root tests, we need a fast and flexible OLS regression allowing to get the parameters (regressor coefficients) solution of the multiple linear equation y = X.b, as well as their variances to compute the t-statistics. These statistics will be used by the unit-root tests to decide whether the serie is (weakly) stationary or not.

- ### Constructor
     The OLS regression is run by declaring an OLS object using the following constructor:
     ```
     OLS<T>::OLS(const Vector<T>& y, const Matrix<T>& x, bool stats = false)
     ```
     with:
     - *y* the vector of the dependent variable
     - *x* the matrix of the independent variables (it can include intercept, constant trend, etc...)
     - *stats* if turned to *true*, additional statistics will be computed as R-squared, adjusted R-squared, F statistic and Durbin-Watson statistic
     
- ### Member variables
    - *param* to get the regressors coefficients
    - *t_stat* to get the theregressors t-statistics
    - *resid* to get the regression residuals
    - *var* to get the regressors variances
    - *R2* to get the R-squared
    - *adj_R2* to get the adjusted R-squared
    - *F_stat* to get the F-statistic
    - *DW_stat* to get the Durbin-Watson statistic
    
- ### Member functions
    - *get_stats()* taking the same Vector *x* and Matrix *y* as arguments, it computes the additional OLS regression statistics
    - *show()* outputs the results
    
- ### Additionnal tools
URT provides 3 functions allowing to add quickly constant terms to a Matrix:
    - *add_intercept()* to insert a column of ones into a Matrix as shown in the example above
    - *add_constant_trend()* to insert a column (1,2,3,...)
    - *add_quadratic_trend()* to insert a column (1,4,9,...)
    
Code example using Armadillo:

```
#include "./URT/include/URT.hpp"

int main()
{
    int nrows = 1000;
    int ncols = 10;

    urt::Vector<double> y = arma::randn<urt::Vector<double>>(nrows);
    urt::Matrix<double> x = arma::randn<urt::Matrix<double>>(nrows, ncols);

    urt::add_intercept(x);

    urt::OLS<double> fit(y, x, true);

    fit.show();
}
    
```
Vector and Matrix are convenient typedefs in the namespace urt, they are alias for vector and matrix representation of the linear algebra in use.

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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
