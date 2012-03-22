include ../../common/make.config

SRC = gaussian.cu
EXE = gaussian
targets=$(EXE)

release: $(SRC)
	$(NVCC) $(SRC) -o $(EXE) -L$(CUDA_LIB_DIR) 

clean:
	rm -f gaussian

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
