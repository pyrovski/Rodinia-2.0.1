include ../../common/make.config
targets=backprop

CC_LINK_FLAGS += -lm

backprop: backprop.o facetrain.o imagenet.o backprop_kernel.o 
	$(CC_OMP_LINK) $^ -o $@

%.o: %.[ch]
	$(CC_OMP) $< -c

facetrain.o: facetrain.c backprop.h
	$(CC_OMP) facetrain.c -c

backprop.o: backprop.c backprop.h
	$(CC_OMP) backprop.c -c

backprop_kernel.o: backprop_kernel.c backprop.h
	$(CC_OMP) backprop_kernel.c -c

imagenet.o: imagenet.c backprop.h
	$(CC_OMP) imagenet.c -c

clean:
	rm -f *.o *~ backprop backprop_cuda.linkinfo

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
