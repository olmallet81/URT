# makefile for building URT as a shared library
# example:
# $ make USE_OPENMP=1 USE_BLAZE=1 USE_BLAS=1
# This command will build URT with OpenMP enabled and Blaze C++ linear algebra library using OpenBLAS for BLAS/LAPACK routines
# Important: when using Intel MKL, MKL environment for bash under Linux has to be set up to allow the compiler to find the path to MKL include and lib directories by adding in the .bashrc file: . /path/to/intel/mkl/bin/mklvars.sh intel64 mod lp64 (to be adapted to your own path depending on where your Intel MKL installation is located)

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
CXXFLAGS = -std=c++11 -Wall -fpic -O3 -march=native

ifeq ($(USE_EIGEN), 1)
  CXXFLAGS += -DUSE_EIGEN
else ifeq ($(USE_BLAZE), 1)
  CXXFLAGS += -DUSE_BLAZE 
else
  USE_ARMA = 1
  CXXFLAGS += -DUSE_ARMA
endif

# adjust the following paths to your own paths, only necessary when linking dynamically and libraries are not on the compiler path, or when linking statically

# setting up Intel MKL libraries path
MKLPATH = /opt/intel/mkl/lib/intel64/
# setting up OpenBLAS library path
BLASPATH = /usr/local/lib/

# linker flags
ifeq ($(USE_MKL), 1)
  CXXFLAGS += -DUSE_MKL
  # with dynamic link to Intel MKL
  LIBFLAGS = -L${MKLPATH} -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -ldl
  # or with static link to Intel MKL
  #LIBFLAGS = -Wl,--start-group ${MKLPATH}libmkl_gf_lp64.a ${MKLPATH}libmkl_core.a ${MKLPATH}libmkl_sequential.a -Wl,--end-group -lpthread -ldl
else ifeq ($(USE_BLAS), 1) 
  CXXFLAGS += -DUSE_BLAS
  # with dynamic link to OpenBLAS
  LIBFLAGS = -L$(BLASPATH) -lopenblas
  # or with static link to OpenBLAS
  #LIBFLAGS = $(BLASPATH)libopenblas.a
else
  ifeq ($(USE_ARMA), 1)
    LIBFLAGS = -larmadillo
  else ifeq ($(USE_BLAZE), 1)
    LIBFLAGS = -llapack
  endif
endif

# NB: if neither Intel MKL nor OpenBLAS are enabled, Armadillo needs to be linked to its shared library to get its own BLAS/LAPACK routines and Blaze needs to be linked to LAPACK library, only Eigen does not need to be linked to any library to work

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
	@$(CXX) $(OPT) $(CXXFLAGS) -c $^ -o $@
	@echo "----- compiled "$(notdir $<)" into "$(notdir $@)" successfully!"

.PHONY: clean deep_clean

# remove obj directory and object files
clean:
	@$(RM) $(OBJ_DIR)

# remove obj directory, its content and shared library
deep_clean:
	@$(RM) $(OBJ_DIR) $(TARGET)
