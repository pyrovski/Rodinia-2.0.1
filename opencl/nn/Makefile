include ../../common/make.config

SRC = nearestNeighbor.cpp clutils.cpp utils.cpp

EXE = nn
targets=nn

$(EXE): $(SRC)
	$(OCL_LINK) $^ -o $@

clean: 
	rm -f $(EXE) 

install: $(targets)
	install -m $(install_perms) -t $(OCL_BIN_DIR) $(targets)
