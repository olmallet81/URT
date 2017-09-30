//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#include "../include/URT.hpp"

namespace urt {

//=================================================================================================

// parameter constructor
template <class T>
OLS<T>::OLS(const Vector<T>& y, const Matrix<T>& x, bool stats)
{ 
   // controlling x and y dimensions
   #ifdef USE_ARMA
   if (y.n_elem != x.n_rows) {
   #elif defined(USE_BLAZE) || defined(USE_EIGEN)
   if (y.size() != x.rows()) {
   #endif
      throw std::invalid_argument("\n  Error: in OLS::OLS(), dependent and independent variables in OLS regression must have same number of rows.\n\n");
   }
   #ifdef USE_ARMA
   else if (!y.n_elem || !x.n_rows) {
   #elif defined(USE_BLAZE) || defined(USE_EIGEN)
   else if (!y.size() || !x.rows()) {
   #endif
      throw std::invalid_argument("\n  Error: in OLS::OLS(), y vector and x matrix must not be empty.\n\n");
   }
 
   #ifdef USE_ARMA
   // number of observations
   this->nobs = x.n_rows;
   // number of regressors
   this->nreg = x.n_cols;
   #elif defined(USE_BLAZE) || defined(USE_EIGEN)
   // number of observations
   this->nobs = x.rows();
   // number of regressors
   #ifdef USE_BLAZE
   this->nreg = x.columns(); 
   #elif USE_EIGEN
   this->nreg = x.cols();
   #endif    
   #endif

   // error number of degrees of freedom 
   this->ndf = nobs - nreg;

   try
   {
      #ifdef USE_ARMA
      // inverting matrix x square
      Matrix<T> ix2 = (x.t() * x).i();
      // computing vector b solution of y = x * b
      this->param = ix2 * (x.t() * y);
      #elif USE_BLAZE
      // inverting matrix x square
      Matrix<T> ix2 = blaze::inv(blaze::trans(x) * x);
      // computing vector b solution of y = x * b
      this->param = ix2 * (blaze::trans(x) * y);
      #elif USE_EIGEN
      // inverting matrix x square
      Matrix<T> tmp = Matrix<T>::Zero(nreg, nreg);
      tmp.template selfadjointView<Eigen::Lower>().rankUpdate(x.transpose());
      tmp.template triangularView<Eigen::Upper>() = tmp.transpose().template triangularView<Eigen::Lower>();
      Matrix<T> ix2 = tmp.inverse();
      // computing vector b solution of y = x * b
      this->param = ix2 * (x.transpose() * y);
      #endif

      // computing residuals 
      this->resid = y - x * param;

      // computing sum of squared residuals
      #ifdef USE_ARMA
      this->SSR = arma::as_scalar(resid.t() * resid);
      #elif USE_BLAZE
      this->SSR = blaze::trans(resid) * resid;
      #elif USE_EIGEN
      this->SSR = resid.dot(resid);
      #endif

      // computing error variance (mean squares error or mse)
      this->MSE = SSR / ndf;

      #ifdef USE_ARMA
      // computing regressors variances
      this->var = diagvec(ix2) * MSE;
      // computing regressors t-statistics
      this->t_stat = param / sqrt(var);
      #elif USE_BLAZE
      // computing regressors variances
      this->var.resize(nreg);
      for (int i = 0; i < nreg; ++i) {
         this->var[i] = ix2(i, i) * MSE;
      }
      // computing regressors t-statistics
      this->t_stat = param * invsqrt(var);
      #elif USE_EIGEN
      // computing regressors variances
      this->var = ix2.diagonal() * MSE;
      // computing regressors t-statistics
      this->t_stat = param.cwiseProduct(var.cwiseSqrt().cwiseInverse());
      #endif
   }
   // catching any exception caused by arma::inv()
   catch (std::exception& e) {
      std::cout << e.what() << "\n";  
      throw std::invalid_argument("\n  Error: in OLS::OLS(), matrix inversion failed in OLS regression, please check x matrix.\n\n");
   }

   this->stats = stats;  
   // computing regression statistics
   if (stats) {
      this->get_stats(y, x);
   }
}

//*************************************************************************************************

// compute OLS statistics (R2, adj_R2, F and DW stats)
template <typename T>
void OLS<T>::get_stats(const Vector<T>& y, const Matrix<T>& x)
{
   // number of variables
   this->nvar = nreg;
        
   // testing if x contains an intercept (some stats are computed differently wether the regression contains an intercept or not)
   for (int i = 0; i < nreg; ++i) {
      // testing if the sum of elements if the ith column is equal to the number of rows
      #ifdef USE_ARMA
      if (arma::sum(x.col(i)) == nobs) {
      #elif USE_BLAZE
      if (x(0, i) == 1 && blaze::isUniform(column(x, i))) {
      #elif USE_EIGEN
      if ((x.col(i)).sum() == nobs) {
      #endif
         // setting OLS intercept to true
         this->intercept = true;
         // recording intercept column index
         this->j0 = i;
         // adjusting number of variables
         this->nvar -= 1;

         break;
      }
   }

   T SST;
   // computing x * b
   Vector<T> z = y - resid;

   if (intercept) {
      #ifdef USE_ARMA
      T my = mean(y);
      Vector<T> diff = y - my;
      SST = as_scalar(diff.t() * diff);
      z -= my;
      #elif USE_BLAZE
      T my = std::accumulate(&y[0], &y[y.size()], 0.0) / y.size();
      Vector<T> diff = forEach(y, [&my](const T& val){ return val - my; });
      SST = blaze::trans(diff) * diff;
      z = forEach(z, [&my](const T& val){ return val - my; });
      #elif USE_EIGEN
      T my = y.mean();
      Vector<T> diff = y - Vector<T>::Constant(y.size(), 1, my);;
      SST = diff.dot(diff);
      z -= Vector<T>::Constant(z.size(), 1, my);
      #endif
   } else {
      #ifdef USE_ARMA
      SST = arma::as_scalar(y.t() * y);
      #elif USE_BLAZE
      SST = blaze::trans(y) * y;
      #elif USE_EIGEN
      SST = y.dot(y);
      #endif
   }

   // computing R-squared
   this->R2 = 1 - SSR / SST;
   // computing adjusted R-squared
   this->adj_R2 = R2 - (1.0 - R2) * nvar / ndf;

   // computing mean of squares for model
   #ifdef USE_ARMA
   T MSM = arma::as_scalar(z.t() * z) / nvar;
   #elif USE_BLAZE
   T MSM = (blaze::trans(z) * z) / nvar;
   #elif USE_EIGEN
   T MSM = z.dot(z) / nvar;
   #endif

   // computing F-statistic
   this->F_stat = MSM / MSE;

   // computing Durbin-Watson statistic
   #ifdef USE_ARMA
   Vector<T> temp = Vector<T>(resid.memptr() + 1, resid.n_elem - 1, false) - Vector<T>(resid.memptr(), resid.n_elem - 1, false);
   this->DW_stat = arma::as_scalar(temp.t() * temp) / SSR;
   #elif USE_BLAZE
   Vector<T> temp = subvector(resid, 1, resid.size() - 1) - subvector(resid, 0, resid.size() - 1);
   this->DW_stat = (blaze::trans(temp) * temp) / SSR; 
   #elif USE_EIGEN
   Vector<T> temp = resid.segment(1, resid.size() - 1) - resid.segment(0, resid.size() - 1);
   this->DW_stat = temp.dot(temp) / SSR;
   #endif

   this->stats = true;
}

//*************************************************************************************************

// print results
template <class T>
void OLS<T>::show()
{
   Vector<T> pval(nreg);

   // declaring Student's distribution
   boost::math::students_t tdistrib(ndf);

   // declaring function computing t-statistic p-value
   std::function<T(T)> P = [&tdistrib](const T& z){ return 2 * (1 - boost::math::cdf(tdistrib, fabs(z))); };

   for (int j = 0; j < nreg; ++j) {
      if (!std::isnan(t_stat[j])) {
         pval[j] = P(t_stat[j]);
      } 
      else {
         pval[j] = -std::nan("");
      }
   }
    
   int gap = 2;

   (intercept) ? gap += 9 : gap += 1 + std::to_string(nreg).length();

   std::string line1("  ");
   for (int i = 0; i < gap + 39; ++i) {
      line1 += "-";
   }

   std::string line2("  ");
   for (int i = 0; i < 40; ++i) {
      line2 += "-";
   }

   // outputting OLS results
   std::cout << "\n";
   std::cout << "  OLS Regression Results\n";
   std::cout << "  ======================\n";

   std::cout << "\n";
   std::cout << "  Coefficients\n";
   std::cout << line1 << "\n";

   // outputting table column names
   std::cout << std::setw(gap + 42) << "Estimate  Std Error   t value  P(>|t|)\n";

   // constructing vector of column indexes
   std::vector<int> idx(nreg);

   // initializing idx with intercept column index
   int i, n, k = 0;

   (intercept) ? idx[k++] = j0, i = 0 : i = 1;
    
   for (int j = 0; j < nreg; ++j) {
      if (!intercept || (intercept && j != j0)) {
         idx[k++] = j;
      }
   }

   for (const auto& j : idx)
   { 
      if (intercept && !i) {
         std::cout << std::setw(gap) << "Intercept";
         ++i;
      } else {
         std::cout << std::setw(gap) << "x" + std::to_string(i++);
      }
      // outputting parameter estimates
      n = std::to_string(int(fabs(param[j]))).length();

      if (fabs(param[j]) >= 1e6 || fabs(param[j]) < 1e-5) {
         std::cout << std::scientific; 
         std::cout << std::setprecision(1);
      } else {
         std::cout << std::fixed; 
         std::cout << std::setprecision(6 - n);
      }
      std::cout << std::setw(11) << param[j];
      std::cout.unsetf(std::ios_base::scientific);
      std::cout.unsetf(std::ios_base::fixed);

      // outputting parameter standard errors
      n = std::to_string(int(sqrt(var[j]))).length();

      if (sqrt(var[j]) >= 1e6 || sqrt(var[j]) < 1e-4) {
         std::cout << std::scientific;
         std::cout << std::setprecision(1); 
      } else {
         std::cout << std::fixed; 
         std::cout << std::setprecision(std::max(5 - n, 0));
      }
      std::cout << std::setw(11) << sqrt(var[j]);
      std::cout.unsetf(std::ios_base::scientific);
      std::cout.unsetf(std::ios_base::fixed);

      // outputting parameter t-statistics
      n = std::to_string(int(fabs(t_stat[j]))).length();

      if (fabs(t_stat[j]) >= 1e6 || fabs(t_stat[j]) < 1e-3) {
         std::cout << std::scientific;
         std::cout << std::setprecision(1); 
      } else {
         std::cout << std::fixed; 
         std::cout << std::setprecision(std::max(4 - n, 0));
      }
      std::cout << std::setw(10) << t_stat[j];
      std::cout.unsetf(std::ios_base::scientific);
      std::cout.unsetf(std::ios_base::fixed);

      // outputting parameter p-values
      std::cout << std::fixed << std::setprecision(3) << std::setw(9) << pval[j];
      std::cout << "\n";
   }

   std::cout << "\n";
   std::cout << "  Residuals\n";
   std::cout << line2 << "\n";

   // outputting residuals results
   #ifdef USE_ARMA
   int nres = resid.n_elem;
   #elif defined(USE_BLAZE) || defined(USE_EIGEN)
   int nres = resid.size();
   #endif

   std::cout << std::setw(10) << "Min" 
             << std::setw(8) << "1Q" 
             << std::setw(8) << "Median" 
             << std::setw(8) << "3Q" 
             << std::setw(8) << "Max";
   std::cout << "\n";

   #if defined(USE_ARMA) || defined(USE_BLAZE)
   n = std::to_string(int(max(abs(resid)))).length();
   if (max(abs(resid)) >= 1e3 || max(abs(resid)) < 1e-3) {
   #elif USE_EIGEN
   n = std::to_string(int(resid.cwiseAbs().maxCoeff())).length();
   if (resid.cwiseAbs().maxCoeff() >= 1e3 || resid.cwiseAbs().maxCoeff() < 1e-3) {
   #endif
      std::cout << std::scientific;
      std::cout << std::setprecision(1); 
   } else {
      std::cout << std::fixed; 
      std::cout << std::setprecision(std::max(4 - n, 0));
   }

   // sorting residuals
   Vector<T> sorted_resid(resid);
   std::sort(&sorted_resid[0], &sorted_resid[nres]);

   std::cout << std::setw(10) << sorted_resid[0];
   std::cout << std::setw(8)  << sorted_resid[round(0.25 * nres) - 1];
   std::cout << std::setw(8)  << sorted_resid[round(0.50 * nres) - 1];
   std::cout << std::setw(8)  << sorted_resid[round(0.75 * nres) - 1];
   std::cout << std::setw(8)  << sorted_resid[nres - 1];

   std::cout << "\n\n";
   std::cout << "  Dimensions\n";
   std::cout << line2 << "\n";
   std::cout << "  Number of observations        " << std::setw(10) << nobs << "\n";
   std::cout << "  Number of degrees of freedom  " << std::setw(10) << ndf << "\n";

   if (stats) {
      std::cout << "  Number of variables           " << std::setw(10) << nvar << "\n";
   }

   std::cout << "\n";
   std::cout << "  Analysis\n";
   std::cout << line2 << "\n";

   // computing residual mean
   #ifdef USE_ARMA
   T mu = mean(resid);
   #elif USE_BLAZE
   T mu = std::accumulate(&resid[0], &resid[resid.size()], 0.0) / resid.size();
   #elif USE_EIGEN
   T mu = resid.mean();
   #endif

   // outputting residual mean
   n = std::to_string(int(fabs(mu))).length();

   if (fabs(mu) >= 1e3 || fabs(mu) < 1e-3) {
      std::cout << std::scientific;
      std::cout << std::setprecision(1); 
   } else {
      std::cout << std::fixed; 
      std::cout << std::setprecision(std::max(4 - n, 0));
   }

   std::cout << "  Residual mean                 " << std::setw(10) << mu << "\n";

   // outputting residual standard error
   n = std::to_string(int(sqrt(MSE))).length();

   if (sqrt(MSE) >= 1e3 || sqrt(MSE) < 1e-3) {
      std::cout << std::scientific;
      std::cout << std::setprecision(1); 
   } else {
      std::cout << std::fixed; 
      std::cout << std::setprecision(std::max(4 - n, 0));
   }

   std::cout << "  Residual standard error       " << std::setw(10) << sqrt(MSE) << "\n";

   // if OLS regression statistics have been computed previously
   if (stats) {
      std::cout << std::fixed;
      // outputting multiple R-squared and adjusted R-squared
      std::cout << "  Multiple R-squared            " << std::setw(10) << std::setprecision(5) << R2 << "\n";
      std::cout << "  Adjusted R-squared            " << std::setw(10) << std::setprecision(5) << adj_R2 << "\n";
      
      // outputting F-statistic
      n = std::to_string(int(F_stat)).length();

      if (sqrt(MSE) >= 1e4 || sqrt(MSE) < 1e-3) {
         std::cout << std::scientific;
         std::cout << std::setprecision(1); 
      } else {
         std::cout << std::fixed; 
         std::cout << std::setprecision(std::max(4 - n, 0));
      }

      std::cout << "  F-statistic (p-value)         " << std::setw(10) << std::setprecision(2) << F_stat << " ";
      std::cout << "(" << 1 - boost::math::ibeta(nvar/2, ndf/2, F_stat*nvar/(F_stat*nvar + ndf)) << ")\n";
      
      // outputting Durbin-Watson statistic
      std::cout << std::fixed;
      std::cout << "  Durbin-Watson statistic       " << std::setw(10) << std::setprecision(3) << DW_stat << "\n";  
   }
   std::cout << "\n";
}

//*************************************************************************************************

// cout overload for printing results from a multi-linear regression
template <class T>
std::ostream& operator<<(std::ostream& out, urt::OLS<T>& result)
{
    result.show();

    return out;
}

//=================================================================================================

}


