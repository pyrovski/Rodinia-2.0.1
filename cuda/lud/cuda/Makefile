include ../../../common/make.config

NVCC_FLAGS += -use_fast_math -I../common

DEFS += -DGPU_TIMER 

CC_FLAGS += -I../common -Wall $(DEFS) 
# Add source files here
EXECUTABLE  := lud_cuda
# Cuda source files (compiled with cudacc)
CUFILES     := lud_kernel.cu
# C/C++ source files (compiled with gcc / c++)
CCFILES     := lud.c lud_cuda.c ../common/common.c

targets=$(EXECUTABLE)

OBJS = ../common/common.o lud.o lud_kernel.o

.PHONY: all clean 
all : $(EXECUTABLE)

.c.o : 
	$(NVCC_OBJ) -o $@ -c $<

%.o:	%.cu 
	$(NVCC_OBJ) -o $@ -c $<


$(EXECUTABLE) : $(OBJS)
	$(NVCC_LINK) -o $@  $?

clean:
	rm -f $(EXECUTABLE) $(OBJS) *.linkinfo

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
