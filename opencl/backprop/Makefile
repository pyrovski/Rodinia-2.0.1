include ../../common/make.config

# C compiler
CC = $(CXX)
targets=backprop

backprop: backprop.o imagenet.o facetrain.o backprop_ocl.cpp
	$(OCL_LINK) $^ -o $@

%.o: %.[ch]
	$(CC) $(CC_FLAGS) $< -c

backprop.o: backprop.c 
	$(CC) $(CC_FLAGS) backprop.c -c

facetrain.o: facetrain.c
	$(CC) $(CC_FLAGS) facetrain.c -c

imagenet.o: imagenet.c
	$(CC) $(CC_FLAGS) imagenet.c -c

clean:
	rm -f *.o *~ backprop *.linkinfo

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
