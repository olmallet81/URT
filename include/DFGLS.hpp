//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#ifndef DFGLS_HPP
#define DFGLS_HPP

#include "Coeff_dfgls.hpp"

namespace urt {

//=================================================================================================

// Dickey-Fuller Generalized Least-Squares test
template <typename T>
class DFGLS : public UnitRoot<T>
{
 public:   
    // parameter constructor for computing DF-GLS test for a given number of lags
    DFGLS(const Vector<T>& data, int lags, const std::string& trend = "c", bool regression = false);
    // parameter constructor for computing ADF test with lag length optimization
    DFGLS(const Vector<T>& data, const std::string& method, const std::string& trend = "c", bool regression = false);
    // compute test statistic
    const T& statistic() override;
    // compute test p-value
    const T& pvalue() override;
    // output test results
    void show() override;

 private:
    using ur = UnitRoot<T>;
    const char* _test_name = "Dickey-Fuller GLS";
    const std::vector<std::string> _valid_trends{"c","ct"};
};

//=================================================================================================

}

#endif


