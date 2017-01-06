cimport cython

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from CyBlaze cimport _unaligned
from CyBlaze cimport _unpadded
from CyBlaze cimport _columnVector
from CyBlaze cimport _columnMajor

from CyBlaze cimport CustomVector
from CyBlaze cimport CustomMatrix
from CyBlaze cimport DynamicVector

cdef extern from "../include/URT.hpp" namespace "urt":

    cdef cppclass OLS[T]:
        # members
        DynamicVector[T] param
        DynamicVector[T] t_stat
        DynamicVector[T] resid
        DynamicVector[T] var
        T MSE
        T R2
        T adj_R2
        T F_stat
        T DW_stat
        T IC
        # constructors
        OLS()
        OLS(const CustomVector[T,_unaligned,_unpadded,_columnVector]&, const CustomMatrix[T,_unaligned,_unpadded,_columnMajor]&, bool) except +
        # method
        void show() except +

    cdef cppclass UnitRoot[T]:
        # members
        string lags_type
        string method
        string test_type
        string trend
        T level
        int lags
        int max_lags
        int niter
        bool bootstrap
        bool regression
        # method
        const T& get_stat() const
        const T& get_pval() const
        const OLS[T]& get_ols() const
        const vector[string]& get_trends() const
        const T& statistic() except +
        const T& pvalue() except +
        void show() except +

    cdef cppclass ADF[T](UnitRoot[T]):
        # constructors
        ADF(const CustomVector[T,_unaligned,_unpadded,_columnVector]&, int, const char*, bool) except +    
        ADF(const CustomVector[T,_unaligned,_unpadded,_columnVector]&, const char*, const char*, bool) except +  

    cdef cppclass DFGLS[T](UnitRoot[T]):
        # constructors
        DFGLS(const CustomVector[T,_unaligned,_unpadded,_columnVector]&, int, const char*, bool) except +   
        DFGLS(const CustomVector[T,_unaligned,_unpadded,_columnVector]&, const char*, const char*, bool) except +  

    cdef cppclass PP[T](UnitRoot[T]):
        # constructorss
        PP(const CustomVector[T,_unaligned,_unpadded,_columnVector]&, int, const char*, const char*, bool) except +      
        PP(const CustomVector[T,_unaligned,_unpadded,_columnVector]&, const char*, const char*, const char*, bool) except +

    cdef cppclass KPSS[T](UnitRoot[T]):
        # constructors
        KPSS(const CustomVector[T,_unaligned,_unpadded,_columnVector]&, int, const char*) except +     
        KPSS(const CustomVector[T,_unaligned,_unpadded,_columnVector]&, const char*, const char*) except +
