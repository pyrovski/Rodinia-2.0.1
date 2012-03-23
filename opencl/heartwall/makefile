include ../../common/make.config
targets=heartwall

# link objects(binaries) together
heartwall: main.cpp OpenCL.cpp avilib.c avimod.c
	$(OCL_LINK) $^ -o $@

clean:
	@rm -f *.o *.out *~ *\#

install: $(targets)
	install -m $(install_perms) -t $(OCL_BIN_DIR) $(targets)
