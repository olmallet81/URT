//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#include "../include/URT.hpp"

int main()
{
   int nobs = 1000;

   // generating stationary random data
   urt::Vector<double> data = urt::gaussian_noise<double>(nobs);

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
