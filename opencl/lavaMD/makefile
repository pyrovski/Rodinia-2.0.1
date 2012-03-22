include ../../common/make.config
targets=lavaMD

# link objects together
lavaMD:	main.o \
	./kernel/kernel_gpu_opencl_wrapper.o \
	./util/num/num.o \
	./util/timer/timer.o \
	./util/opencl/opencl.o
	$(OCL_LINK) $^\
	-o $@

# compile function files into objects
main.o:	main.h \
	main.c \
	./kernel/kernel_gpu_opencl_wrapper.h \
	./kernel/kernel_gpu_opencl_wrapper.c \
	./util/num/num.h \
	./util/num/num.c \
	./util/timer/timer.h \
	./util/timer/timer.c \
	./util/opencl/opencl.h \
	./util/opencl/opencl.c
	$(OCL_OBJ)	main.c \
	-c \
	-o $@

./kernel/kernel_gpu_opencl_wrapper.o:	./kernel/kernel_gpu_opencl_wrapper.h \
	./kernel/kernel_gpu_opencl_wrapper.c
	$(OCL_OBJ)	./kernel/kernel_gpu_opencl_wrapper.c \
	-c \
	-o $@

./util/num/num.o:	./util/num/num.h \
	./util/num/num.c
	$(OCL_OBJ)	./util/num/num.c \
	-c \
	-o $@

./util/timer/timer.o:	./util/timer/timer.h \
	./util/timer/timer.c
	$(OCL_OBJ)	./util/timer/timer.c \
	-c \
	-o $@

./util/opencl/opencl.o:	./util/opencl/opencl.h \
	./util/opencl/opencl.c
	$(OCL_OBJ) ./util/opencl/opencl.c \
	-c \
	-o $@
# delete all object and executable files
clean:
	rm -f	*.o lavaMD \
	./kernel/*.o \
	./util/num/*.o \
	./util/timer/*.o \
	./util/opencl/*.o

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
