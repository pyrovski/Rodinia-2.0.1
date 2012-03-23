include ../../common/make.config

SRC = pathfinder.cu

EXE = pathfinder

targets=$(EXE)

$(EXE): $(SRC)
	$(NVCC_LINK) $^ -o $@

clean:
	rm -f $(targets)

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
