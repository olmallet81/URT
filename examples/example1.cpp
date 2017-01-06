//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#include "../include/URT.hpp"

int main()
{
   int nrows = 1000;
   int ncols = 10;

   // generating random arrays
   urt::Vector<double> y = urt::wiener_process<double>(nrows);
   urt::Matrix<double> x = urt::wiener_process<double>(nrows, ncols);

   // adding intercept to matrix of independent variables
   urt::add_intercept(x);

   // writting data to CSV files
   urt::WriteToCSV("./URT/data/y.csv", y);
   urt::WriteToCSV("./URT/data/x.csv", x);

   // running OLS regression
   urt::OLS<double> fit(y, x, true);

   // outputting regression results
   fit.show();
}
