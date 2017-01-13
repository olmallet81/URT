//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#include "../include/URT.hpp"

#ifndef USE_ARMA
  #include <armadillo>
#endif

// define USE_FLOAT when compiling to switch to single precision
#ifdef USE_FLOAT
  using T = float;
#else
  using T = double;
#endif

int main()
{
   int niter = 0;

   arma::wall_clock timer;

   std::vector<int> sizes = {100,150,200,250,300,350,400,450,500,1000,1500,2000,2500,3000,3500,4000,4500,5000};

   std::cout << std::fixed << std::setprecision(1);

   for (int i = 0; i < sizes.size(); ++i) {

      urt::Vector<T> data = urt::wiener_process<T>(sizes[i]);

      (sizes[i] < 1000) ? niter = 10000 : niter = 1000;

      timer.tic();
      for (int k = 0; k < niter; ++k) {
         urt::ADF<T> test(data, "AIC");
         test.statistic();
      }

      auto duration = timer.toc();

      std::cout << std::setw(8) << sizes[i];
      std::cout << std::setw(8) << duration << "\n";
   }
}
