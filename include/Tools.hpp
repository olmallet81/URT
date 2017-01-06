//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#ifndef TOOLS_HPP
#define TOOLS_HPP

// unseeded random number generator 
static boost::mt19937 rng;

namespace urt {

//=================================================================================================

// generate Vector of normally distributed random data 
template <typename T>
Vector<T> gaussian_noise(int n, T mu = 0.0, T sigma = 1.0)
{
   boost::normal_distribution<T> ndistrib(mu, sigma);
   boost::variate_generator<boost::mt19937&, boost::normal_distribution<T>> rnorm(rng, ndistrib);

   Vector<T> x(n);

   for (int i = 0; i < n; ++i)
      x[i] = rnorm();

   return x;
}

//*************************************************************************************************

// generate Matrix of normally distributed random data
template <typename T>
Matrix<T> gaussian_noise(int nr, int nc, T mu = 0.0, T sigma = 1.0)
{
   boost::normal_distribution<T> ndistrib(mu, sigma);
   boost::variate_generator<boost::mt19937&, boost::normal_distribution<T>> rnorm(rng, ndistrib);

   Matrix<T> x(nr, nc);

   for (int i = 0; i < nr; ++i)
      for (int j = 0; j < nc; ++j)
         x(i, j) = rnorm();

   return x;
}

//*************************************************************************************************

// generate Wiener process
template <typename T>
Vector<T> wiener_process(int n, T mu = 0.0, T sigma = 1.0)
{
   Vector<T> x(n);

   Vector<T> z = gaussian_noise<T>(n, mu, sigma);

   std::partial_sum(&z[0], &z[n], &x[0]);

   return x;
}

//*************************************************************************************************

// generate Matrix of Wiener processes
template <typename T>
Matrix<T> wiener_process(int nr, int nc, T mu = 0.0, T sigma = 1.0)
{
   Matrix<T> x(nr, nc);

   Matrix<T> z = gaussian_noise<T>(nr, nc, mu, sigma);

   for (int j = 0; j < nc; ++j)
      std::partial_sum(&z(0, j), &z(nr, j), &x(0, j));

   return x;
}

//*************************************************************************************************

// << operator overload for outputting std::vector
template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& x)
{
   for (auto it = x.cbegin(), end = x.cend(); it != end; ++it)
      std::cout << *it << " ";

   std::cout << "\n";

   return out;
}

//*************************************************************************************************

// modify matrix x by adding an intercept
template <typename T>
void add_intercept(Matrix<T>& x)
{
   #ifdef USE_ARMA
   x.insert_cols(x.n_cols, arma::ones<arma::Col<T>>(x.n_rows));
   #elif USE_BLAZE
   x.resize(x.rows(), x.columns() + 1);
   column(x, x.columns() - 1) = forEach(column(x, x.columns() - 1), [](T val){ return 1; });
   #elif USE_EIGEN
   x.conservativeResize(Eigen::NoChange, x.cols() + 1);
   x.col(x.cols() - 1) = Vector<T>::Ones(x.rows());
   #endif
}

//*************************************************************************************************

// modify matrix x by adding a constant trend
template <typename T>
void add_constant_trend(Matrix<T>& x)
{
   add_intercept(x);

   #ifdef USE_ARMA
   x.insert_cols(x.n_cols, arma::cumsum(x.col(x.n_cols - 1)));
   #elif USE_BLAZE
   x.resize(x.rows(), x.columns() + 1);
   std::partial_sum(&x(0, x.columns() - 2), &x(x.rows() - 1, x.columns() - 2), &x(0, x.columns() - 1));
   #elif USE_EIGEN
   x.conservativeResize(Eigen::NoChange, x.cols() + 1);
   x.col(x.cols() - 1) = Vector<T>::LinSpaced(x.rows(), 0, x.rows() - 1);
   #endif
}

//*************************************************************************************************

// modify matrix x by adding a quadratic trend
template <typename T>
void add_quadratic_trend(Matrix<T>& x)
{
   add_constant_trend(x);

   #ifdef USE_ARMA
   Vector<T> z(x.col(x.n_cols - 1));
   x.insert_cols(x.n_cols, z % z);
   #elif USE_BLAZE
   column(x, x.columns() - 1) = forEach(column(x, x.columns() - 2), [](T val){ return val * val; });
   #elif USE_EIGEN
   Vector<T> z(x.col(x.cols() - 1));
   x.conservativeResize(Eigen::NoChange, x.cols() + 1);
   x.col(x.cols() - 1) = z.cwiseProduct(z);
   #endif
}

//=================================================================================================

}

#endif
