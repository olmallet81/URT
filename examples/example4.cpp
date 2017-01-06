//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#include "../include/URT.hpp"

int main()
{
   // generating non-stationary random data
   urt::Vector<double> data = urt::wiener_process<double>(1000);

   // initializing Phillips-Perron normalized test with lags of type long and constant term
   urt::PP<double> test(data, "long", "c", "rho");

   // outputting test results
   test.show();
    
   // switching to t-statistic test 
   test.test_type = "tau";
    
   // outputting test results
   test.show();  
}
