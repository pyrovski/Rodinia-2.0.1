include ../../common/make.config

SRC = gaussian.cu
EXE = gaussian
targets=$(EXE)

$(EXE): $(SRC)
	$(NVCC_LINK) $^ -o $@

clean:
	rm -f gaussian

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
