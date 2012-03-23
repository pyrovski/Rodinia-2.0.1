include ../../../common/make.config
targets=srad_v2
CC=$(CXX)

srad_v2: 
	$(CC_OMP_LINK) srad.cpp -o srad_v2

clean:
	rm -f $(targets)

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
