include ../../common/make.config

TARGET_G = sc_gpu
targets=$(TARGET_G)

$(TARGET_G): streamcluster_cuda_cpu.cpp streamcluster_cuda.cu streamcluster_header.cu
	$(NVCC) $(NVCC_FLAGS) streamcluster_cuda_cpu.cpp streamcluster_cuda.cu streamcluster_header.cu -o $(TARGET_G) -lcuda


clean:
	rm -f *.o *~ *.txt $(TARGET_G) *.linkinfo

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
