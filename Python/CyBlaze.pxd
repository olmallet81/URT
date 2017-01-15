import numpy as np
cimport numpy as np
cimport cython

from libc.string cimport memcpy
from libcpp cimport bool
from libcpp.memory cimport unique_ptr
from cython.operator cimport dereference as deref

#--------------------------------------------------------------------------------------------------

# Blaze C++ constant unaligned exposed to Cython
cdef extern from "blaze/math/AlignmentFlag.h" namespace "blaze":
    ctypedef bool _unaligned "blaze::unaligned"

# Blaze C++ constant unpadded exposed to Cython
cdef extern from "blaze/math/PaddingFlag.h" namespace "blaze":
    ctypedef bool _unpadded "blaze::unpadded"

# Blaze C++ constant columnVector exposed to Cython
cdef extern from "blaze/math/TransposeFlag.h" namespace "blaze":
    ctypedef bool _columnVector "blaze::columnVector"

# Blaze C++ constant columnMajor exposed to Cython
cdef extern from "blaze/math/StorageOrder.h" namespace "blaze":
    ctypedef bool _columnMajor "blaze::columnMajor"

#--------------------------------------------------------------------------------------------------

# Blaze C++ class CustomVector exposed to Cython
cdef extern from "blaze/math/CustomVector.h" namespace "blaze" nogil:
    cdef cppclass CustomVector[T,_unaligned,_unpadded,_columnVector]:
        CustomVector() nogil
        CustomVector(T*, int ) nogil

# convenient typedef
ctypedef CustomVector[double,_unaligned,_unpadded,_columnVector] VectorCd
ctypedef CustomVector[float,_unaligned,_unpadded,_columnVector] VectorCf

#--------------------------------------------------------------------------------------------------

# Blaze C++ class CustomMatrix exposed to Cython
cdef extern from "blaze/math/CustomMatrix.h" namespace "blaze" nogil:
    cdef cppclass CustomMatrix[T,_unaligned,_unpadded,_columnMajor]:
        CustomMatrix() nogil
        CustomMatrix(T*, int, int ) nogil

# convenient typedef
ctypedef CustomMatrix[double,_unaligned,_unpadded,_columnMajor] MatrixCd
ctypedef CustomMatrix[float,_unaligned,_unpadded,_columnMajor] MatrixCf

#--------------------------------------------------------------------------------------------------

# Blaze C++ class DynamicVector exposed to Cython
cdef extern from "blaze/math/DynamicVector.h" namespace "blaze" nogil:
    cdef cppclass DynamicVector[T]:
        DynamicVector() nogil
        int size() nogil
        const T& operator[](int ) nogil

# convenient typedef
ctypedef DynamicVector[double] VectorDd
ctypedef DynamicVector[float] VectorDf

#--------------------------------------------------------------------------------------------------

# Blaze DynamicVector[double] wrapped into a Numpy 1-dim array without copying the data
cdef inline np.ndarray[np.double_t,ndim=1] blaze_vec_to_numpy(VectorDd& xd):
    cdef np.npy_intp shape[1]
    shape[0] = <np.npy_intp> xd.size()
    cdef np.ndarray[np.double_t,ndim=1] x = np.PyArray_SimpleNewFromData(1,shape,np.NPY_DOUBLE,&xd[0])
    return x

#--------------------------------------------------------------------------------------------------

# Blaze DynamicVector[float] wrapped into a Numpy 1-dim array without copying the data
cdef inline np.ndarray[np.float32_t,ndim=1] blaze_fvec_to_numpy(VectorDf& xf):
    cdef np.npy_intp shape[1]
    shape[0] = <np.npy_intp> xf.size()
    cdef np.ndarray[np.float32_t,ndim=1] x = np.PyArray_SimpleNewFromData(1,shape,np.NPY_FLOAT32,&xf[0])
    return x

#--------------------------------------------------------------------------------------------------

