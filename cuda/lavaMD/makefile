include ../../common/make.config

OMP_LIB = -lgomp
OMP_FLAG = -fopenmp
targets=lavaMD

# link objects (binaries) together
lavaMD:	main.o \
	./kernel/kernel_gpu_cuda_wrapper.o \
	./util/num/num.o \
	./util/timer/timer.o \
	./util/device/device.o
	$(CC)	main.o \
	./kernel/kernel_gpu_cuda_wrapper.o \
	./util/num/num.o \
	./util/timer/timer.o \
	./util/device/device.o \
	-lm \
	-L$(CUDA_LIB_DIR) \
	-lcuda -lcudart \
	$(OMP_LIB) \
	-o lavaMD

# compile function files into objects (binaries)
main.o:	main.h \
	main.c \
	./kernel/kernel_gpu_cuda_wrapper.h \
	./kernel/kernel_gpu_cuda_wrapper.cu \
	./util/num/num.h \
	./util/num/num.c \
	./util/timer/timer.h \
	./util/timer/timer.c \
	./util/device/device.h \
	./util/device/device.cu
	$(CC)	main.c \
	-c \
	-o main.o \
	-O3

./kernel/kernel_gpu_cuda_wrapper.o:	./kernel/kernel_gpu_cuda_wrapper.h \
	./kernel/kernel_gpu_cuda_wrapper.cu
	$(NVCC)	./kernel/kernel_gpu_cuda_wrapper.cu \
	-c \
	-o ./kernel/kernel_gpu_cuda_wrapper.o \

./util/num/num.o:	./util/num/num.h \
	./util/num/num.c
	$(CC)	./util/num/num.c \
	-c \
	-o ./util/num/num.o \
	-O3

./util/timer/timer.o:	./util/timer/timer.h \
	./util/timer/timer.c
	$(CC)	./util/timer/timer.c \
	-c \
	-o ./util/timer/timer.o \
	-O3

./util/device/device.o:	./util/device/device.h \
	./util/device/device.cu
	$(NVCC)	./util/device/device.cu \
	-c \
	-o ./util/device/device.o \
	-O3

# delete all object and executable files
clean:
	rm -f	*.o \
	./kernel/*.o \
	./util/num/*.o \
	./util/timer/*.o \
	./util/device/*.o \
	lavaMD

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
