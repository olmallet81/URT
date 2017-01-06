//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

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
