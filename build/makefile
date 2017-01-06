# makefile for building URT as a shared library
# example:
# $ make USE_OPENMP=1 USE_BLAZE=1 USE_BLAS=1
# will build URT with parallelism active and with Blaze C++ linear algebra library using OpenBLAS for BLAS/LAPACK routines

# NB: when using Intel MKL, MKL environment for bash has to be set up to allow the compiler to find the path to MKL include and lib directories by adding: . /opt/intel/mkl/bin/mklvars.sh intel64 mod lp64

# compiler
CXX = g++
# source directory
SRC_DIR = ../src
# objects directory
OBJ_DIR = ./obj
# shared library to be created
TARGET = ../lib/libURT.so
# remove files
RM = rm -rf core

# add OpenMP if required
ifeq ($(USE_OPENMP), 1)
  CXX += -fopenmp
endif

# compiler flags
CXXFLAGS = -std=c++14 -Wall -fpic -O3 -march=native

ifeq ($(USE_EIGEN), 1)
  CXXFLAGS += -DUSE_EIGEN
else ifeq ($(USE_BLAZE), 1)
  CXXFLAGS += -DUSE_BLAZE 
else
  USE_ARMA = 1
  CXXFLAGS += -DUSE_ARMA
endif

# setting up MKL libs path
MKLPATH = /opt/intel/mkl/lib/intel64/
# setting up OpenBLAS lib path
BLASPATH = /usr/local/lib/

# linker flags
ifeq ($(USE_MKL), 1)
  CXXFLAGS += -DUSE_MKL
  LIBFLAGS = -L${MKLPATH} -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -ldl
  #LIBFLAGS = -Wl,--start-group ${MKLPATH}libmkl_gf_lp64.a ${MKLPATH}libmkl_core.a ${MKLPATH}libmkl_sequential.a -Wl,--end-group -ldl
else ifeq ($(USE_BLAS), 1) 
  CXXFLAGS += -DUSE_BLAS
  LIBFLAGS = -L$(BLASPATH) -lopenblas
  #LIBFLAGS = $(BLASPATH)libopenblas.a
else
  ifeq ($(USE_ARMA), 1)
    LIBFLAGS = -larmadillo
  else ifeq ($(USE_BLAZE), 1)
    LIBFLAGS = -llapack
  endif
endif

# NB: if neither Intel MKL nor OpenBLAS are called, Armadillo needs to link to its shared library to get its internal BLAS/LAPACK wrappers and Blaze needs to link to LAPACK library, only Eigen does not need to link to any library to work

# get .cpp source files and name .o object files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))

# build shared library
$(TARGET): $(OBJS)
	@$(CXX) -shared -o $(TARGET) $^ $(LIBFLAGS)
	@echo "----- built shared library "$(notdir $(TARGET))" successfully!"

# compile source files into object files
$(OBJS): $(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	@$(CXX) $(CXXFLAGS) -c $^ -o $@
	@echo "----- compiled "$(notdir $<)" into "$(notdir $@)" successfully!"

.PHONY: clean deep_clean

# remove obj directory and object files
clean:
	@$(RM) $(OBJ_DIR)

# remove obj directory, its content and shared library
deep_clean:
	@$(RM) $(OBJ_DIR) $(TARGET)
