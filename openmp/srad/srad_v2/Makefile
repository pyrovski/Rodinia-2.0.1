include ../../../common/make.config
targets=srad_v2

srad_v2: 
	$(CC) $(CC_FLAGS) srad.cpp -o srad_v2

clean:
	rm -f $(targets)

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
