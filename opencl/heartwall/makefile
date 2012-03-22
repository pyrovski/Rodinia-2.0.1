include ../../common/make.config

CC=g++

# link objects(binaries) together
a.out: main.o OpenCL.o avilib.o avimod.o
	$(OCL_LINK) $^ -o $@

# compile main function file into object (binary)
main.o: main.cpp OpenCL.cpp kernels.cl avilib.c avimod.c
	$(OCL_OBJ) -o $@ $^

clean:
	@rm -f *.o *.out *~ *\#
