//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

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
