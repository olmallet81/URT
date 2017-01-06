//=================================================================================================
//                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#ifndef KPSS_HPP
#define KPSS_HPP

#include "Coeff_kpss.hpp"

namespace urt {

//=================================================================================================

// Kwiatkowski–Phillips–Schmidt–Shin test
template <class T>
class KPSS : public UnitRoot<T>
{
 public:
    // parameter constructor for computing KPSS test for a given number of lags
    KPSS(const Vector<T>& data, int lags, const std::string& trend = "c");
    // parameter constructor for computing PP test for a default number of lags (long or short)
    KPSS(const Vector<T>& data, const std::string lags_type, const std::string& trend = "c");
    // compute KPSS test statistic
    const T& statistic() override;
    // compute KPSS test p-value
    const T& pvalue() override;
    // output KPSS test results
    void show() override;

 private:
    using ur = UnitRoot<T>;
    const char* _test_name = "KPSS";
    const std::vector<std::string> _valid_trends{"c","ct"};
    // compute test statistic
    void compute_stat();
};

//=================================================================================================

}

#endif


