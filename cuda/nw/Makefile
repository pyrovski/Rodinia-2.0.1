include ../../common/make.config

SRC = needle.cu

EXE = needle
targets=$(EXE)

$(EXE): $(SRC)
	$(NVCC) $(SRC) -o $(EXE) -L$(CUDA_LIB_DIR) 

clean: $(SRC)
	rm -f $(EXE) $(EXE).linkinfo result.txt

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
