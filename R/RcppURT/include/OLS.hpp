//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#ifndef OLS_HPP
#define OLS_HPP

namespace urt {

//=================================================================================================

// Multi-Linear Regression by Ordinary Least Squares
// NB: we chose never to make hard copies of x and y for performance reason, indeed in ADF test with lag length optimization we can quickly get a very large x matrix
template <class T>
class OLS
{
 public:
    Vector<T> param;   // regressors parameters
    Vector<T> t_stat;  // regressors statistics
    Vector<T> resid;   // residuals
    Vector<T> var;     // regressors variances

    T MSE;                // mean of squares for error
    T SSR;                // sum of squares residuals
    T R2 = 0;             // R-squared
    T adj_R2 = 0;         // adjusted R-squared
    T F_stat = 0;         // F statistic
    T DW_stat = 0;        // Durbin-Watson statistic
    T IC = 0;             // model information criterion (for Unit-Root tests)
 
    int nobs;             // number of observations
    int nreg;             // number of regressors
    int ndf;              // number of degrees of freedom
    int lags = 0;         // number of lags (for Unit-Root tests)

    // default (nullary) constructor
    OLS() {}
    // parameter constructor
    OLS(const Vector<T>& y, const Matrix<T>& x, bool stats = false);
    // compute OLS statistics (R2, adj_R2, F and DW stats)
    void get_stats(const Vector<T>& y, const Matrix<T>& x);
    // print results
    void show();

 private:
    int j0;                  // intercept index in matrix of independent variables
    int nvar = 0;            // number of independent variables (regressors excluding intercept)
    bool intercept = false;  // control if OLS is with intercept or without
    bool stats = false;      // control if OLS statistics need to be computed
};

//=================================================================================================

}

#endif


