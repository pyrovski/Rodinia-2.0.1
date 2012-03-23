include ../../common/make.config

SRC = clutils.cpp gaussianElim.cpp utils.cpp

EXE = gaussian
targets=$(EXE)

$(EXE): $(SRC)
	$(OCL_LINK) -o $@ $^

clean: 
	rm -f $(EXE) 

install: $(targets)
	install -m $(install_perms) -t $(OCL_BIN_DIR) $(targets)
