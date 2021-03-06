#
# The following environment variables should be predefined:
#
# OP2_COMPILER (gnu,intel,etc)
#
# ... and optionally the following. If not defined, then 
# they must be locatable via PATH variables:
#
# OP2_INSTALL_PATH
# CUDA_INSTALL_PATH
# PARMETIS_INSTALL_PATH
# PTSCOTCH_INSTALL_PATH
# HDF5_INSTALL_PATH
#

include ../../../../make-common.inc

ifeq ($(OP2_COMPILER),gnu)
  CPP		= g++
  CPPFLAGS	= -g -fPIC -DUNIX -Wall
#  CPPFLAGS	= -O3 -fPIC -DUNIX -Wall -Wextra
  OMPFLAGS	= -fopenmp
  MPICPP	= $(MPICPP_PATH)
  MPIFLAGS	= $(CCFLAGS)
else
ifeq ($(OP2_COMPILER),intel)
  CPP		= icpc
  CPPFLAGS	= -O3 -xHost #-DMPICH_IGNORE_CXX_SEEK#-DCOMM_PERF #-DDEBUG #-vec-report
  OMPFLAGS	= -qopenmp #-openmp-report2
  MPICPP	= $(MPICXX_PATH)
  MPIFLAGS	= $(CPPFLAGS)
else
ifeq ($(OP2_COMPILER),clang)
  CPP       	= clang++
  CCFLAGS  	= -O3 #-ffast-math
  CPPFLAGS 	= $(CCFLAGS)
  OMPFLAGS 	= -I$(OMPTARGET_LIBS)/../tools/openmp/runtime/src/ -fopenmp=libomp -Rpass-analysis
  OMPOFFLOAD 	= -fopenmp-targets=nvptx64-nvidia-cuda -fopenmp-nonaliased-maps -ffp-contract=fast -Xcuda-ptxas -v #-Xclang -target-feature -Xclang +ptx35
  MPICC   	= $(MPICC_PATH)
  MPICPP	= $(MPICXX_PATH)
  MPIFLAGS 	= $(CPPFLAGS)
else
print:
	@echo "unrecognised value for OP2_COMPILER"
endif
endif
endif

#
# set flags for NVCC compilation and linking
#

ifndef NV_ARCH
  MESSAGE=select an NVIDA device to compile in CUDA, e.g. make NV_ARCH=KEPLER
  NV_ARCH=Kepler
endif
ifeq ($(NV_ARCH),Fermi)
  CODE_GEN_CUDA=-gencode arch=compute_20,code=sm_21
else
ifeq ($(NV_ARCH),Kepler)
  CODE_GEN_CUDA=-gencode arch=compute_35,code=sm_35
endif
endif

NVCCFLAGS += -O3 $(CODE_GEN_CUDA) -m64 -Xptxas -dlcm=ca -Xptxas=-v -use_fast_math #-g -G -O0

#
# master to make all versions
#

all: clean aero_seq aero_openmp aero_cuda aero_mpi aero_mpi_openmp aero_mpi_cuda

#
# simple sequential version
#

aero_seq: aero.cpp Makefile
	     $(MPICPP) $(CPPFLAGS) aero.cpp $(OP2_INC) $(HDF5_INC) $(OP2_LIB) -lop2_seq -lop2_hdf5 $(HDF5_LIB) -o aero_seq


#
# x86 version using kernel files generated by op2.m
#

aero_openmp:	aero_op.cpp openmp/aero_kernels.cpp \
                Makefile
		$(MPICPP) $(VAR) $(CPPFLAGS) $(OMPFLAGS) $(OP2_INC) $(HDF5_INC) -Iopenmp -I. \
		$(OP2_LIB) aero_op.cpp openmp/aero_kernels.cpp -lm -lop2_openmp -lop2_hdf5 $(HDF5_LIB) -o aero_openmp

#
# CUDA version using kernel files generated by op2.m
#

aero_cuda:	aero_op.cpp cuda/aero_kernels_cu.o Makefile
		$(MPICPP) $(VAR) $(CPPFLAGS) aero_op.cpp cuda/aero_kernels_cu.o \
		$(CUDA_INC) $(OP2_INC) $(HDF5_INC) \
		$(OP2_LIB) $(CUDA_LIB) -lcudart -lop2_cuda -lop2_hdf5 $(HDF5_LIB) -o aero_cuda -DCUDA

cuda/aero_kernels_cu.o:	cuda/aero_kernels.cu      \
                Makefile
		nvcc  $(INC) $(NVCCFLAGS) $(OP2_INC) $(HDF5_INC) \
		-Icuda -I.  -c -o cuda/aero_kernels_cu.o \
		cuda/aero_kernels.cu -DCUDA

#
# mpi with sequential-nodes version
#

aero_mpi: aero.cpp Makefile
	$(MPICPP) $(MPIFLAGS) aero.cpp $(OP2_INC) $(PARMETIS_INC) $(PTSCOTCH_INC) $(HDF5_INC) \
	$(OP2_LIB) -lop2_mpi $(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) -o aero_mpi


#
# mpi with openmp-nodes version
#
aero_mpi_openmp:	aero_op.cpp openmp/aero_kernels.cpp \
                Makefile
		$(MPICPP) $(VAR) $(MPIFLAGS) $(OMPFLAGS) $(OP2_INC) $(PARMETIS_INC) $(PTSCOTCH_INC) $(HDF5_INC) -Iopenmp -I. \
		$(OP2_LIB) aero_op.cpp openmp/aero_kernels.cpp -lm -lop2_mpi $(HDF5_LIB) $(PARMETIS_LIB) $(PTSCOTCH_LIB) \
		-o aero_mpi_openmp

#
# mpi with CUDA version
#
aero_mpi_cuda:	aero_op.cpp cuda/aero_kernels_mpi_cu.o Makefile
		$(MPICPP) $(VAR) $(MPIFLAGS) aero_op.cpp cuda/aero_kernels_mpi_cu.o \
		$(OP2_INC) $(PARMETIS_INC) $(PTSCOTCH_INC) $(HDF5_INC) \
		$(OP2_LIB) -lop2_mpi_cuda $(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) \
		$(CUDA_LIB) -lcudart -o aero_mpi_cuda

cuda/aero_kernels_mpi_cu.o:	cuda/aero_kernels.cu \
                Makefile
		nvcc  $(INC) $(NVCCFLAGS) $(OP2_INC) -I $(MPI_INSTALL_PATH)/include \
		-Icuda -I. -c -o cuda/aero_kernels_mpi_cu.o \
		cuda/aero_kernels.cu


#
# OpenMP4 version
#
aero_openmp4: aero_op.cpp openmp4/aero_omp4kernels.cpp openmp4/aero_omp4kernel_funcs.cpp Makefile
	$(CPP) $(VAR) $(CPPFLAGS) $(OMPFLAGS) $(OMPOFFLOAD) $(OP2_INC) $(HDF5_INC) $(OMP4_REG_COUNT)\
		-Iopenmp4 -I. -c openmp4/aero_omp4kernel_funcs.cpp  -o openmp4/aero_omp4kernel_funcs.o
	$(CPP) $(VAR) $(CPPFLAGS) $(OMPFLAGS) $(OP2_INC) $(HDF5_INC) $(OMP4_REG_COUNT)\
		-Iopenmp4 -I. -c openmp4/aero_omp4kernels.cpp -o openmp4/aero_omp4kernels.o
	$(MPICPP) $(VAR) $(CPPFLAGS) $(OMPFLAGS) $(OMPOFFLOAD) $(OP2_INC) $(OP2_LIB) $(HDF5_INC)  \
		aero_op.cpp openmp4/aero_omp4kernels.o openmp4/aero_omp4kernel_funcs.o -o aero_openmp4 \
		-lop2_openmp4 -lop2_hdf5 $(HDF5_LIB) -L$(CUDA_INSTALL_PATH)/lib64 -lcudart


write_hdf5: write_hdf5.cpp
	$(MPICPP) $(MPIFLAGS) write_hdf5.cpp $(OP2_INC) $(PARMETIS_INC) $(PTSCOTCH_INC) $(HDF5_INC) \
	$(OP2_LIB) -lop2_mpi $(PARMETIS_LIB) $(PTSCOTCH_LIB) $(HDF5_LIB) -o write_hdf5


#
# cleanup
#

clean:
		rm -f aero_seq aero_openmp aero_cuda aero_mpi write_hdf5 aero_mpi_openmp aero_mpi_cuda aero_openmp4 *.o
