include ../../../common/make.config

SRC = srad.cu

EXE = srad
targets=$(EXE)

$(EXE): $(SRC)
	$(NVCC_LINK) $^ -o $@

clean: $(SRC)
	rm -f $(EXE) $(EXE).linkinfo result.txt

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
