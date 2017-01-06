//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#ifndef ADF_HPP
#define ADF_HPP

#include "Coeff_adf.hpp"

namespace urt {

//=================================================================================================

// Augmented Dickey-Fuller test 
template <typename T>
class ADF : public UnitRoot<T> 
{
 public:
   // parameter constructor for computing ADF test for a given number of lags
   ADF(const Vector<T>& data, int lags, const std::string& trend = "c", bool regression = false);
   // parameter constructor for computing ADF test with lag length optimization
   ADF(const Vector<T>& data, const std::string& method, const std::string& trend = "c", bool regression = false);
   // compute test statistic
   const T& statistic() override;
   // compute test p-value
   const T& pvalue() override;
   // output test results
   void show() override;

 private:
   using ur = UnitRoot<T>;
   const char* _test_name = "Augmented Dickey-Fuller";
   const std::vector<std::string> _valid_trends{"nc","c","ct","ctt"};
}; 

//=================================================================================================

}

#endif
