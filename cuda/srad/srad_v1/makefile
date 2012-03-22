include ../../../common/make.config

targets = srad

# link objects(binaries) together
srad:	main.o
	$(CC)	main.o \
	-L$(CUDA_LIB_DIR)  \
	-lm -lcuda -lcudart -o srad

# compile main function file into object (binary)
main.o: 	main.cu \
	define.c \
	graphics.c \
	extract_kernel.cu \
	prepare_kernel.cu \
	reduce_kernel.cu \
	srad_kernel.cu \
	srad2_kernel.cu \
	compress_kernel.cu
	$(NVCC)	main.cu \
	-c -O3 -arch sm_13

# delete all object files
clean:
	rm -f *.o srad

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
