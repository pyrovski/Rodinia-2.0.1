include ../../common/make.config
EXE = hotspot
targets=$(EXE)

$(EXE): hotspot.c hotspot.h OpenCL_helper_library.c OpenCL_helper_library.h
	$(OCL_LINK) $(^:%.h=) -o $@

clean:
	rm -f $(EXE)

install: $(targets)
	install -m $(install_perms) -t $(OCL_BIN_DIR) $(targets)
