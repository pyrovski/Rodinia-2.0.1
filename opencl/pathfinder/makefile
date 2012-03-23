include ../../common/make.config
targets=pathfinder

# link objects(binaries) together
pathfinder: main.o OpenCL.o
	$(OCL_LINK) $^ -o $@

# compile main function file into object (binary)
main.o: main.cpp OpenCL.cpp kernels.cl
	$(OCL_LINK) main.cpp OpenCL.cpp -o $@

clean:
	@rm -f pathfinder *.o *.out *~ *\#

install: $(targets)
	install -m $(install_perms) -t $(OCL_BIN_DIR) $(targets)
