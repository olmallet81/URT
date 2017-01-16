import numpy as np
cimport numpy as np 
cimport cython

from CyBlaze cimport VectorCd   
from CyBlaze cimport MatrixCd
from CyBlaze cimport VectorCf
from CyBlaze cimport MatrixCf    

from CyBlaze cimport blaze_vec_to_numpy
from CyBlaze cimport blaze_fvec_to_numpy

from libcpp.memory cimport unique_ptr
from cython.operator cimport dereference as deref

# essential to initialize numpy array, not doing so will cause seg fault when trying to wrap Blaze arrays
np.import_array()

#-------------------------------------------------------------------------------------------------------#

class pyOLS:
    def __init__(self, param, t_stat, resid, var, MSE, R2, adj_R2, F_stat, DW_stat, IC):
        self.param = np.copy(param)
        self.t_stat = np.copy(t_stat)
        self.resid = np.copy(resid)
        self.var = np.copy(var)
        self.MSE = MSE
        self.R2 = R2
        self.adj_R2 = adj_R2
        self.F_stat = F_stat
        self.DW_stat = DW_stat
        self.IC = IC

#-------------------------------------------------------------------------------------------------------#

# C++ class OLS<double> exposed to Python
cdef class OLS_d(object):

    # member
    cdef unique_ptr[OLS[double]] ols

    def __init__(self, np.ndarray[np.double_t, ndim = 1] y not None, np.ndarray[np.double_t, ndim = 2] x not None, bool stats = False):

        # converting numpy arrays to Blaze CustomVector and CustomMatrix (no copy)
        cdef VectorCd yd = VectorCd(<double*> y.data, y.shape[0])
        cdef MatrixCd xd = MatrixCd(<double*> x.data, x.shape[0], x.shape[1])

        # running OLS regression
        self.ols = unique_ptr[OLS[double]](new OLS[double](yd, xd, stats))

    # methods
    @property
    def param(self):
        return blaze_vec_to_numpy(deref(self.ols).param)

    @property
    def t_stat(self):
        return blaze_vec_to_numpy(deref(self.ols).t_stat)

    @property
    def resid(self):
        return blaze_vec_to_numpy(deref(self.ols).resid)

    @property
    def var(self):
        return blaze_vec_to_numpy(deref(self.ols).var)

    @property
    def MSE(self):
        return deref(self.ols).MSE

    @property
    def R2(self):
        return deref(self.ols).R2

    @property
    def adj_R2(self):
        return deref(self.ols).adj_R2

    @property
    def F_stat(self):
        return deref(self.ols).F_stat

    @property
    def DW_stat(self):
        return deref(self.ols).DW_stat

    def show(self):
        deref(self.ols).show() 


#-------------------------------------------------------------------------------------------------------#

# C++ class UnitRoot<double> exposed to Python
cdef class UnitRoot_d:

    # members
    cdef UnitRoot[double]* ur

    # contructor
    def __init__(self):
        pass

    def statistic(self):
        return self.ur.statistic()

    def pvalue(self):
        return self.ur.pvalue()

    def show(self):
        self.ur.show()

    @property 
    def stat(self):
        return self.ur.get_stat();

    @property 
    def pval(self):
        return self.ur.get_pval();

    @property
    def ols(self):
        cdef OLS[double] result = self.ur.get_ols()
        return pyOLS(blaze_vec_to_numpy(result.param),
                     blaze_vec_to_numpy(result.t_stat),
                     blaze_vec_to_numpy(result.resid),
                     blaze_vec_to_numpy(result.var),
                     result.MSE,
                     result.R2,
                     result.adj_R2,
                     result.F_stat,
                     result.DW_stat,
                     result.IC)

    @property
    def trends(self):
        return np.array(self.ur.get_trends())

    # methods
    @property
    def bootstrap(self):
        return self.ur.bootstrap

    @bootstrap.setter
    def bootstrap(self, bool value):
        if self.ur.bootstrap != value:
            self.ur.bootstrap = value

    @property
    def lags(self):
        return self.ur.lags

    @lags.setter
    def lags(self, int value):
        if self.ur.lags != value:
            self.ur.lags = value

    @property
    def lags_type(self):
        return self.ur.lags_type

    @lags_type.setter
    def lags_type(self, bytes value):
        if self.ur.lags_type != value:
            self.ur.lags_type = value

    @property
    def level(self):
        return self.ur.level

    @level.setter
    def level(self, float value):
        if self.ur.level != value:
            self.ur.level = value

    @property
    def max_lags(self):
        return self.ur.max_lags

    @max_lags.setter
    def max_lags(self, int value):
        if self.ur.max_lags != value:
            self.ur.max_lags = value

    @property
    def method(self):
        return self.ur.method

    @method.setter
    def method(self, bytes value):
        if self.ur.method != value:
            self.ur.method = value

    @property
    def niter(self):
        return self.ur.niter

    @niter.setter
    def niter(self, int value):
        if self.ur.niter != value:
            self.ur.niter = value

    @property
    def regression(self):
        return self.ur.regression

    @regression.setter
    def regression(self, bool value):
        if self.ur.regression != value:
            self.ur.regression = value

    @property
    def test_type(self):
        return self.ur.test_type

    @test_type.setter
    def test_type(self, bytes value):
        if self.ur.test_type != value:
            self.ur.test_type = value

    @property
    def trend(self):
        return self.ur.trend

    @trend.setter
    def trend(self, bytes value):
        if self.ur.trend != value:
            self.ur.trend = value

#-------------------------------------------------------------------------------------------------------#

# C++ class ADF<double> exposed to Python
cdef class ADF_d(UnitRoot_d):

    cdef unique_ptr[ADF[double]] adf

    # constructor
    def __init__(self, np.ndarray[np.double_t, ndim = 1] x not None, lags = None, str method = 'AIC', str trend = 'c', bool regression = False):
        # converting numpy array to Blaze CustomVector (no copy)
        cdef VectorCd xd = VectorCd(<double*> x.data, x.shape[0])
 
        if lags is None:
            self.adf = unique_ptr[ADF[double]](new ADF[double](xd, method, trend, regression))
            self.ur = self.adf.get()
        else:
            self.adf = unique_ptr[ADF[double]](new ADF[double](xd, <int>lags, trend, regression))
            self.ur = self.adf.get()

#-------------------------------------------------------------------------------------------------------#

# C++ class DFGLS<double> exposed to Python
cdef class DFGLS_d(UnitRoot_d):

    cdef unique_ptr[DFGLS[double]] dfgls

    # constructor
    def __init__(self, np.ndarray[np.double_t, ndim = 1] x not None, lags = None, str method = 'AIC', str trend = 'c', bool regression = False):

        # converting numpy array to Blaze CustomVector (no copy)
        cdef VectorCd xd = VectorCd(<double*> x.data, x.shape[0])

        if lags is None:
            self.dfgls = unique_ptr[DFGLS[double]](new DFGLS[double](xd, method, trend, regression))
            self.ur = self.dfgls.get()
        else:
            self.dfgls = unique_ptr[DFGLS[double]](new DFGLS[double](xd, <int>lags, trend, regression))
            self.ur = self.dfgls.get()

#-------------------------------------------------------------------------------------------------------#

# C++ class PP<double> exposed to Python
cdef class PP_d(UnitRoot_d):

    cdef unique_ptr[PP[double]] pp

    # constructor
    def __init__(self, np.ndarray[np.double_t, ndim = 1] x not None, lags = None, str lags_type = 'long', str trend = 'c', str test_type = 'tau', bool regression = False):

        # converting numpy array to Blaze CustomVector (no copy)
        cdef VectorCd xd = VectorCd(<double*> x.data, x.shape[0])

        if lags is None:
            self.pp = unique_ptr[PP[double]](new PP[double](xd, lags_type, trend, test_type, regression))
            self.ur = self.pp.get()
        else:
            self.pp = unique_ptr[PP[double]](new PP[double](xd, <int>lags, trend, test_type, regression))
            self.ur = self.pp.get()

#-------------------------------------------------------------------------------------------------------#

# C++ class KPSS<double> exposed to Python
cdef class KPSS_d(UnitRoot_d):

    cdef unique_ptr[KPSS[double]] kpss

    # constructor
    def __init__(self, np.ndarray[np.double_t, ndim = 1] x not None, lags = None, str lags_type = 'long', str trend = 'c'):

        # converting numpy array to Blaze CustomVector (no copy)
        cdef VectorCd xd = VectorCd(<double*> x.data, x.shape[0])

        if lags is None:
            self.kpss = unique_ptr[KPSS[double]](new KPSS[double](xd, lags_type, trend))
            self.ur = self.kpss.get()
        else:
            self.kpss = unique_ptr[KPSS[double]](new KPSS[double](xd, <int>lags, trend))
            self.ur = self.kpss.get()

#-------------------------------------------------------------------------------------------------------#
#########################################################################################################
#-------------------------------------------------------------------------------------------------------#

# C++ class OLS<float> exposed to Python
cdef class OLS_f(object):

    # member
    cdef unique_ptr[OLS[float]] ols

    def __init__(self, np.ndarray[np.float32_t, ndim = 1] y not None, np.ndarray[np.float32_t, ndim = 2] x not None, bool stats = False):

        # converting numpy arrays to Blaze CustomVector and CustomMatrix (no copy)
        cdef VectorCf yf = VectorCf(<float*> y.data, y.shape[0])
        cdef MatrixCf xf = MatrixCf(<float*> x.data, x.shape[0], x.shape[1])

        # running OLS regression
        self.ols = unique_ptr[OLS[float]](new OLS[float](yf, xf, stats))

    # methods
    @property
    def param(self):
        return blaze_fvec_to_numpy(deref(self.ols).param)

    @property
    def t_stat(self):
        return blaze_fvec_to_numpy(deref(self.ols).t_stat)

    @property
    def resid(self):
        return blaze_fvec_to_numpy(deref(self.ols).resid)

    @property
    def var(self):
        return blaze_fvec_to_numpy(deref(self.ols).var)

    @property
    def MSE(self):
        return deref(self.ols).MSE

    @property
    def R2(self):
        return deref(self.ols).R2

    @property
    def adj_R2(self):
        return deref(self.ols).adj_R2

    @property
    def F_stat(self):
        return deref(self.ols).F_stat

    @property
    def DW_stat(self):
        return deref(self.ols).DW_stat

    def show(self):
        deref(self.ols).show()

#-------------------------------------------------------------------------------------------------------#

# C++ class UnitRoot<float> exposed to Python
cdef class UnitRoot_f:

    # members
    cdef UnitRoot[float]* ur

    # contructor
    def __init__(self):
        pass

    def statistic(self):
        return self.ur.statistic()

    def pvalue(self):
        return self.ur.pvalue()

    def show(self):
        self.ur.show()

    @property 
    def stat(self):
        return self.ur.get_stat();

    @property 
    def pval(self):
        return self.ur.get_pval();

    @property
    def ols(self):
        cdef OLS[float] result = self.ur.get_ols()
        return pyOLS(blaze_fvec_to_numpy(result.param),
                     blaze_fvec_to_numpy(result.t_stat),
                     blaze_fvec_to_numpy(result.resid),
                     blaze_fvec_to_numpy(result.var),
                     result.MSE,
                     result.R2,
                     result.adj_R2,
                     result.F_stat,
                     result.DW_stat,
                     result.IC)

    @property
    def trends(self):
        return np.array(self.ur.get_trends())

    # methods
    @property
    def bootstrap(self):
        return self.ur.bootstrap

    @bootstrap.setter
    def bootstrap(self, bool value):
        if self.ur.bootstrap != value:
            self.ur.bootstrap = value

    @property
    def lags(self):
        return self.ur.lags

    @lags.setter
    def lags(self, int value):
        if self.ur.lags != value:
            self.ur.lags = value

    @property
    def lags_type(self):
        return self.ur.lags_type

    @lags_type.setter
    def lags_type(self, bytes value):
        if self.ur.lags_type != value:
            self.ur.lags_type = value

    @property
    def level(self):
        return self.ur.level

    @level.setter
    def level(self, float value):
        if self.ur.level != value:
            self.ur.level = value

    @property
    def max_lags(self):
        return self.ur.max_lags

    @max_lags.setter
    def max_lags(self, int value):
        if self.ur.max_lags != value:
            self.ur.max_lags = value

    @property
    def method(self):
        return self.ur.method

    @method.setter
    def method(self, bytes value):
        if self.ur.method != value:
            self.ur.method = value

    @property
    def niter(self):
        return self.ur.niter

    @niter.setter
    def niter(self, int value):
        if self.ur.niter != value:
            self.ur.niter = value

    @property
    def regression(self):
        return self.ur.regression

    @regression.setter
    def regression(self, bool value):
        if self.ur.regression != value:
            self.ur.regression = value

    @property
    def test_type(self):
        return self.ur.test_type

    @test_type.setter
    def test_type(self, bytes value):
        if self.ur.test_type != value:
            self.ur.test_type = value

    @property
    def trend(self):
        return self.ur.trend

    @trend.setter
    def trend(self, bytes value):
        if self.ur.trend != value:
            self.ur.trend = value

#-------------------------------------------------------------------------------------------------------#

# C++ class ADF<float> exposed to Python
cdef class ADF_f(UnitRoot_f):

    cdef unique_ptr[ADF[float]] adf

    # constructor
    def __init__(self, np.ndarray[np.float32_t, ndim = 1] x not None, lags = None, str method = 'AIC', str trend = 'c', bool regression = False):
        # converting numpy array to Blaze CustomVector (no copy)
        cdef VectorCf xf = VectorCf(<float*> x.data, x.shape[0])

        if lags is None:
            self.adf = unique_ptr[ADF[float]](new ADF[float](xf, method, trend, regression))
            self.ur = self.adf.get()
        else:
            self.adf = unique_ptr[ADF[float]](new ADF[float](xf, <int> lags, trend, regression))
            self.ur = self.adf.get()

#-------------------------------------------------------------------------------------------------------#

# C++ class DFGLS<float> exposed to Python
cdef class DFGLS_f(UnitRoot_f):

    cdef unique_ptr[DFGLS[float]] dfgls

    # constructor
    def __init__(self, np.ndarray[np.float32_t, ndim = 1] x not None, lags = None, str method = 'AIC', str trend = 'c', bool regression = False):

        # converting numpy array to Blaze CustomVector (no copy)
        cdef VectorCf xf = VectorCf(<float*> x.data, x.shape[0])

        if lags is None:
            self.dfgls = unique_ptr[DFGLS[float]](new DFGLS[float](xf, method, trend, regression))
            self.ur = self.dfgls.get()
        else:
            self.dfgls = unique_ptr[DFGLS[float]](new DFGLS[float](xf, <int>lags, trend, regression))
            self.ur = self.dfgls.get()

#-------------------------------------------------------------------------------------------------------#

# C++ class PP<float> exposed to Python
cdef class PP_f(UnitRoot_f):

    cdef unique_ptr[PP[float]] pp

    # constructor
    def __init__(self, np.ndarray[np.float32_t, ndim = 1] x not None, lags = None, str lags_type = 'long', str trend = 'c', str test_type = 'tau', bool regression = False):

        # converting numpy array to Blaze CustomVector (no copy)
        cdef VectorCf xf = VectorCf(<float*> x.data, x.shape[0])

        if lags is None:
            self.pp = unique_ptr[PP[float]](new PP[float](xf, lags_type, trend, test_type, regression))
            self.ur = self.pp.get()
        else:
            self.pp = unique_ptr[PP[float]](new PP[float](xf, <int>lags, trend, test_type, regression))
            self.ur = self.pp.get()

#-------------------------------------------------------------------------------------------------------#

# C++ class KPSS<float> exposed to Python
cdef class KPSS_f(UnitRoot_f):

    cdef unique_ptr[KPSS[float]] kpss

    # constructor
    def __init__(self, np.ndarray[np.float32_t, ndim = 1] x not None, lags = None, str lags_type = 'long', str trend = 'c'):

        # converting numpy array to Blaze CustomVector (no copy)
        cdef VectorCf xf = VectorCf(<float*> x.data, x.shape[0])

        if lags is None:
            self.kpss = unique_ptr[KPSS[float]](new KPSS[float](xf, lags_type, trend))
            self.ur = self.kpss.get()
        else:
            self.kpss = unique_ptr[KPSS[float]](new KPSS[float](xf, <int>lags, trend))
            self.ur = self.kpss.get()

#-------------------------------------------------------------------------------------------------------#


