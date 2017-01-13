//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#ifndef CSVMANAGER_HPP
#define CSVMANAGER_HPP

namespace urt {

//=================================================================================================

// write Vector to csv file
template <typename T>
void WriteToCSV(const std::string& filename, const Vector<T>& x)
{
   std::ofstream myfile(filename);

   #ifdef USE_ARMA
   int nrows = x.n_rows;
   #elif defined(USE_BLAZE) || defined(USE_EIGEN)
   int nrows = x.size();
   #endif

   for (int i = 0; i < nrows; ++i)
      myfile << x[i] << "\n";

   myfile.close();
}

//*************************************************************************************************

// write Matrix to csv file
template <typename T>
void WriteToCSV(const std::string& filename, const Matrix<T>& x)
{
   std::ofstream myfile(filename);

   #ifdef USE_ARMA
   int nrows = x.n_rows;
   int ncols = x.n_cols;
   #elif USE_BLAZE
   int nrows = x.rows();
   int ncols = x.columns(); 
   #elif USE_EIGEN
   int nrows = x.rows();
   int ncols = x.cols(); 
   #endif

   for (int i = 0; i < nrows; ++i) { 
      for (int j = 0; j < ncols; ++j) {

         myfile << x(i, j);

         if (j != ncols - 1) {
            myfile << ",";
         } else {
            myfile << "\n";
         }
      }
   }

   myfile.close();
}

//*************************************************************************************************

// read data from csv file and store them into a Vector 
template <typename T>
void ReadFromCSV(const std::string& filename, Vector<T>& x)
{
   std::ifstream myfile(filename);

   if (!myfile.good())  {
      std::cerr << "\n  Error: " << filename << " cannot be found!\n\n";  
      return;
   }

   std::string line;

   int nrows = 0;

   // reading data
   while (std::getline(myfile, line)) {
      #if defined(USE_ARMA) || defined(USE_BLAZE)
      x.resize(nrows + 1);
      #elif USE_EIGEN
      x.conservativeResize(nrows + 1, 1);
      #endif

      std::stringstream lineStream(line);
      std::string cell;
      std::getline(lineStream, cell, ',');

      x[nrows] = std::stod(cell);
      ++nrows;
   }
   myfile.close();
}

//*************************************************************************************************

// read data from csv file and store them into a Matrix 
template <typename T>
void ReadFromCSV(const std::string& filename, Matrix<T>& x)
{
   std::ifstream myfile(filename);

   if (!myfile.good())  {
      std::cerr << "\n  Error: " << filename << " cannot be found!\n\n";  
      return;
   }

   std::string line;

   int nrows = 0;

   // reading data
   while (std::getline(myfile, line)) {

      std::stringstream lineStream(line);
      std::string cell;

      #ifdef USE_ARMA
      arma::Row<T> z;
      #elif USE_BLAZE
      blaze::DynamicVector<T,blaze::rowVector> z;
      #elif USE_EIGEN
      Vector<T> z;
      #endif

      int ncols = 0;

      while (std::getline(lineStream, cell, ',')) {
         #if defined(USE_ARMA) || defined(USE_BLAZE)
         z.resize(ncols + 1);
         #elif USE_EIGEN
         z.conservativeResize(ncols + 1, Eigen::NoChange);
         #endif
         z[ncols] = std::stod(cell);
         ++ncols;
      }
      #ifdef USE_ARMA
      x.insert_rows(nrows, z);
      #elif USE_BLAZE
      x.resize(nrows + 1, ncols);
      row(x, nrows) = z;
      #elif USE_EIGEN
      x.conservativeResize(nrows + 1, ncols);
      x.row(nrows) = z;
      #endif
      ++nrows;
   }
   myfile.close();
}

//=================================================================================================

}

#endif
