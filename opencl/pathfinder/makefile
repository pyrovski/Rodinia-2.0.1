include ../../common/make.config
targets=pathfinder

# link objects(binaries) together
pathfinder: main.cpp OpenCL.cpp
	$(OCL_LINK) $^ -o $@


clean:
	@rm -f pathfinder *.o *.out *~ *\#

install: $(targets)
	install -m $(install_perms) -t $(OCL_BIN_DIR) $(targets)
