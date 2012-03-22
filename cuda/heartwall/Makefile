include ../../common/make.config
targets=heartwall

# link objects(binaries) together
heartwall: main.o ./AVI/avilib.o ./AVI/avimod.o 
	$(NVCC) $(NVCC_FLAGS) main.o ./AVI/avilib.o ./AVI/avimod.o -I/usr/local/cuda/include -lm -o heartwall

# compile main function file into object (binary)
main.o: main.cu kernel.cu define.c
	$(NVCC) $(NVCC_FLAGS) main.cu -I./AVI -c

./AVI/avilib.o ./AVI/avimod.o:
	$(MAKE) -C AVI

# delete all object files
clean:
	rm -f *.o AVI/*.o heartwall *.linkinfo

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
