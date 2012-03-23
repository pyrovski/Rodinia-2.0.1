include ../../common/make.config
targets=backprop

CC_FLAGS += -I.
CC_LINK_FLAGS += -lm

backprop: backprop.o facetrain.o imagenet.o backprop_cuda.o 
	$(NVCC_LINK) $^ -o $@

%.o: %.[ch]
	$(NVCC_OBJ) $< -c

facetrain.o: facetrain.c backprop.h
	$(NVCC_OBJ) facetrain.c -c

backprop.o: backprop.c backprop.h
	$(NVCC_OBJ) backprop.c -c

backprop_cuda.o: backprop_cuda.cu backprop.h
	$(NVCC_OBJ) -c backprop_cuda.cu

imagenet.o: imagenet.c backprop.h
	$(NVCC_OBJ) imagenet.c -c

clean:
	rm -f *.o *~ backprop backprop_cuda.linkinfo

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
