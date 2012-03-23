include ../../common/make.config
targets=needle
CC=$(CXX)

needle: needle.cpp
	$(CC_OMP_LINK) -o $@ $<

clean:
	rm -f needle

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)

