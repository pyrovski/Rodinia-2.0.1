include common/make.config

all: CUDA OMP OPENCL

CUDA: 
	$(MAKE) -C cuda

OMP:
	$(MAKE) -C openmp

OPENCL:
	$(MAKE) -C opencl

.PHONY: CUDA_clean OMP_clean OCL_clean

clean: CUDA_clean OMP_clean OCL_clean


CUDA_clean:
	rm -f $(CUDA_BIN_DIR)/*
	$(MAKE) -C cuda clean

OMP_clean:
	rm -f $(OMP_BIN_DIR)/*
	$(MAKE) -C openmp clean

OCL_clean:
	rm -f $(OCL_BIN_DIR)/*
	$(MAKE) -C opencl clean
