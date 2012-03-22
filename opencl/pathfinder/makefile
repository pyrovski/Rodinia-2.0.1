include ../../common/make.config
targets=pathfinder

# link objects(binaries) together
pathfinder: main.o OpenCL.o
	$(OCL_LINK) $^ -o $@

# compile main function file into object (binary)
main.o: main.cpp OpenCL.cpp kernels.cl
	$(CXX) $(CXXFLAGS) main.cpp OpenCL.cpp

clean:
	@rm -f pathfinder *.o *.out *~ *\#

