//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#ifndef PP_HPP
#define PP_HPP

#include "Coeff_pp.hpp"

namespace urt { 

//=================================================================================================

// Phillips-Perron test
template <class T>
class PP : public UnitRoot<T>
{
 public:
   // PP has its own lags variable
   int lags = 0;    
   // parameter constructor for computing PP test for a given number of lags
   PP(const Vector<T>& data, int lags, const std::string& trend = "c", const std::string& test_type = "tau", bool regression = false);
   // parameter constructor for computing PP test for a default number of lags (long or short) 
   PP(const Vector<T>& data, const std::string& lags_type, const std::string& trend = "c", const std::string& test_type = "tau", bool regression = false);
   // compute PP test statistic
   const T& statistic() override;
   // compute PP test p-value
   const T& pvalue() override;
   // output PP test results
   void show() override;

 private:
   using ur = UnitRoot<T>;
   const char* _test_name = "Phillips-Perron";
   const std::vector<std::string> _valid_trends{"nc","c","ct"};       
   // compute test statistic
   void compute_stat();
}; 

//=================================================================================================

}

#endif


