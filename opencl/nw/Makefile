include ../../common/make.config

EXE = nw
targets=$(EXE)

$(EXE): nw.cpp
	$(OCL_LINK) $^ -o $@

clean:
	rm -f $(EXE)

install: $(targets)
	install -m $(install_perms) -t $(OCL_BIN_DIR) $(targets)
