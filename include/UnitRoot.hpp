//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#ifndef UNITROOT_HPP
#define UNITROOT_HPP

namespace urt {

//=================================================================================================

// Abstract base class from which all tests will derive
template <typename T>
class UnitRoot
{
 public:
   std::string lags_type;                 // default lags value long or short
   std::string method;                    // method for lag length optimization  
   std::string test_type;                 // t-statistic or normalized (tau or rho)
   std::string trend;                     // regression trend

   float level = 1.64;                    // statistic (absolute) threshold for optimal lags selection

   int lags = 0;                          // number of lags in model
   int max_lags = 0;                      // maximum number of lags for lag length optimization
   int niter = 1000;                      // number of iterations, required when bootstrap is set to true

   bool bootstrap = false;                // compute critical values and p-value by bootstrap    
   bool regression = false;               // output OLS regression results if set to true

   // get test statistic
   const T& get_stat() const;
   // get test pvalue
   const T& get_pval() const;
   // get test OLS regression results
   const OLS<T>& get_ols() const;
   // get test valid trends
   const std::vector<std::string>& get_trends() const;

   // compute test statistic
   virtual const T& statistic() = 0;
   // compute test statistic p-value
   virtual const T& pvalue() = 0;
   // output test results (can be overriden by derived classes)
   virtual void show() = 0;

   // virtual destructor
   virtual ~UnitRoot() {}

 protected:
   std::shared_ptr<OLS<T>> result;         // test OLS regression result
   std::vector<std::string> valid_trends;  // vector of test valid regression trends
   Vector<T>* ptr = nullptr;               // pointer to data serie

   // pointer to critical values coefficients
   const std::map<float,std::map<int,std::vector<float>>>* coeff_ptr = nullptr; 

   std::string test_name;                 // test name

   T stat = 0.0;                 // test statistic
   T pval = 0.0;                 // test p-value

   int nobs;                     // number of observations
   int npar;                     // number of parameters excluding lag difference terms

   bool new_data = false;        // control if bootstrap has modified original data
   bool new_lags = false;        // control if a new number of lags has been chosen
   bool new_level = false;       // control if new stat level has been chosen
   bool new_method = false;      // control if new method has been chosen
   bool new_test_type = false;   // control if new test type has been chosen
   bool new_trend = false;       // control if a new trend has been chosen

   // constructor for running test for a given number of lags
   UnitRoot(const Vector<T>& data, int lags, const std::string& trend, bool regression = false);
   // constructor for running test with lag length optmization
   UnitRoot(const Vector<T>& data, const std::string& method, const std::string& trend, bool regression = false);

   // set pointer to the original data
   void set_data();
   // set number of lags
   void set_lags();
   // set lags type
   void set_lags_type();
   // set statistic level for optimal lags selection (method T-STAT)
   void set_level();
   // set method for lag length optimization
   void set_method();
   // set number of iterations for bootstrap
   void set_niter();
   // set test type (tau or rho)
   void set_test_type();
   // set regression trend
   void set_trend();
   // OLS demeaning or detrending
   void ols_detrend();
   // GLS demeaning or detrending
   void gls_detrend();
   // run ADF test for a given lag
   void adf_regression();
   // compute ADF test
   void compute_adf();

 private:
   Vector<T> data;               // original data used for the test
   Matrix<T> x;                  // matrix of independent variables for multi-linear regression
   Vector<T> y;                  // vector of dependent variable for multi-linear regression 
   Vector<T> z;                  // vector of detrended data (by OLS or GLS)

   std::string prev_lags_type;   // previous type of lags
   std::string prev_method;      // previous method for lag length optimization
   std::string prev_test_type;   // previous test type (tau or rho)
   std::string prev_trend;       // previous regression trend 
   std::string trend_type;       // regression trend for outputting test results

   T ICcc;                       // information criterion correction coefficient
   float prev_level;             // previous statistic (absolute) value for lag selection

   int nrows;                    // number of rows in matrix x and elements in vector y
   int prev_lags = 0;            // previous number of lags
   int prev_max_lags = 0;        // previous maximum number of lags
   int prev_niter = 1000;        // previous number of iteratons     
           
   bool new_niter = false;       // control if new number of iterations for bootstrap has been chosen
   bool optim = false;           // control if lag length is optimized 
   bool new_test = false;        // control if a new test is run (true) or all parameters remain the same (false)
   bool prev_bootstrap = false;  // previous bootstrap value
   bool prev_regression = false; // previous regression value

   // vector of test results for lag length optimization
   std::vector<std::shared_ptr<OLS<T>>> results; 
   // vector of valid methods for lag length optimization
   const std::vector<std::string> valid_methods{"AIC","BIC","HQC","MAIC","MBIC","MHQC","T-STAT"};
   // array of probabilities for p-value computation
   const float probas[15] = {0.001,0.005,0.01,0.025,0.05,0.10,0.20,0.50,0.80,0.90,0.95,0.975,0.99,0.995,0.999};
   // array containing tne test critical values
   float critical_values[15]; 
   // Information Criterion functor 
   T (UnitRoot<T>::*ICfunc)(const std::shared_ptr<OLS<T>>& res);
   // set Information Criterion function and correction coefficient 
   void set_IC();
   // initialize dependent and independent variables in ADF test regression
   void initialize_adf();
   // optimize ADF test lag length using Information Criterion
   void optimize_lag();
   // select optimal ADF test lag using lags difference terms t-statistics
   void select_lag();
   // compute Information Criterion
   T IC(const std::shared_ptr<OLS<T>>& res);
   // compute Modified Information Criterion
   T MIC(const std::shared_ptr<OLS<T>>& res);
   // compute critical value from probabilities
   void compute_cv();
   // compute p-value by linear interpolation from critical values 
   void compute_pval();
   // bootstrap test residuals for critical values and p-value computation
   void run_bootstrap();
};

//=================================================================================================

}

#endif





