include ../../../common/make.config
targets=$(BINDIR)/mummergpu
all: mummergpu
.SUFFIXES : .cu .cu_dbg_o .c_dbg_o .cpp_dbg_o .cu_rel_o .c_rel_o .cpp_rel_o .cubin

include ../experiments/test_rule.mk

# Warning flags
CXXWARN_FLAGS := \
	-W -Wall \
	-Wimplicit \
	-Wswitch \
	-Wformat \
	-Wchar-subscripts \
	-Wparentheses \
	-Wmultichar \
	-Wtrigraphs \
	-Wpointer-arith \
	-Wcast-align \
	-Wreturn-type \
	-Wno-unused-function

CWARN_FLAGS := $(CXXWARN_FLAGS) \
	-Wstrict-prototypes \
	-Wmissing-prototypes \
	-Wmissing-declarations \
	-Wnested-externs \
	-Wmain 

#CXX_FLAGS  += -m64 $(CXXWARN_FLAGS)
#CC_FLAGS   += -m64 $(CWARN_FLAGS)
OPENACC_FLAGS += -I$(CUDA_DIR)/include

# Compilers
# Add source files here
STATIC_LIB := libmummergpu.a
# Cuda source files (compiled with cudacc)
CUFILES := mummergpu.cu
# C/C++ source files (compiled with gcc / c++)
CXXFILES := mummergpu_gold.cpp suffix-tree.cpp PoolMalloc.cpp

################################################################################
# Rules and targets

# Basic directory setup for SDK
# (override directories only if they are not already defined)
ROOTDIR    ?= ..
ROOTBINDIR ?= $(ROOTDIR)/bin
BINDIR     ?= $(ROOTBINDIR)
ROOTOBJDIR ?= obj
LIBDIR     := $(ROOTDIR)/lib
BINNAME    ?= mummergpu

# Libs
LIB       := -L$(LIBDIR) -lGL -lGLU -lcuda -L$(CUDA_DIR)/lib64 -lcudart

OBJDIR := $(ROOTOBJDIR)

$(OBJDIR)/%.o : %.c
	$(OPENACC_CXX) -o $@ -c $<

$(OBJDIR)/%.o : %.cpp
	$(OPENACC_CXX) -o $@ -c $<

$(OBJDIR)/%.o : %.cu
	$(OPENACC_CXX) -o $@ -c $<

%.cubin : mummergpu.cu mummergpu_kernel.cu
	$(OPENACC_CXX) -o $(BINNAME).cubin -cubin $<

%.ptx : %.cu mummergpu_kernel.cu
	$(OPENACC_CXX) -o $@ -ptx -Xopencc -LIST:source=on $<

makedirectories:
	@mkdir -p $(LIBDIR)
	@mkdir -p $(OBJDIR)
	@mkdir -p $(BINDIR)

LIBOBJS= \
    $(OBJDIR)/mummergpu_gold.o \
    $(OBJDIR)/suffix-tree.o \
    $(OBJDIR)/PoolMalloc.o \
    $(OBJDIR)/mummergpu.o \

$(LIBDIR)/libmummergpu$(LIBSUFFIX).a: $(LIBOBJS)
	ar qv $@ $(LIBOBJS)

BINOBJS= \
    $(OBJDIR)/mummergpu_main.o

$(BINDIR)/mummergpu: makedirectories $(LIBDIR)/libmummergpu$(LIBSUFFIX).a $(BINOBJS)
	$(OPENACC_LINK) -o $(BINDIR)/$(BINNAME)  $(BINOBJS) -lmummergpu$(LIBSUFFIX) $(LIB) 

clean: 
	rm -f $(LIBOBJS)
	rm -f $(LIBDIR)/libmummergpu$(LIBSUFFIX).a
	rm -f $(BINOBJS)
	rm -f $(BINDIR)/$(BINNAME)
	rm -f *.linkinfo
	rm -rf $(ROOTBINDIR)
	rm -rf $(ROOTOBJDIR)
	rm -rf $(LIBDIR)

mummergpu: $(BINDIR)/mummergpu

cubin: mummergpu.cubin

ptx: mummergpu.ptx

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
