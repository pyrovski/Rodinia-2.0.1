include ../../common/make.config
targets=euler3d  euler3d_double  pre_euler3d  pre_euler3d_double 

NVCC_LINK += $(CUTIL_LIB)

all: $(targets)

euler3d: euler3d.cu
	$(NVCC_LINK) $^ -o $@

euler3d_double: euler3d_double.cu
	$(NVCC_LINK) $^ -o $@

pre_euler3d: pre_euler3d.cu
	$(NVCC_LINK) $^ -o $@

pre_euler3d_double: pre_euler3d_double.cu
	$(NVCC_LINK) $^ -o $@

clean:
	rm -f $(targets) *.linkinfo

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
