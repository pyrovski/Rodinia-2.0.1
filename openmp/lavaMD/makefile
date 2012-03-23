include ../../common/make.config
CC_LINK_FLAGS += -lm

targets=lavaMD

# link objects (binaries) together
lavaMD:	main.o \
	./kernel/kernel_cpu.o \
	./util/num/num.o \
	./util/timer/timer.o
	$(CC_OMP_LINK)	$^ -o $@

# compile function files into objects (binaries)
main.o:	main.h \
	main.c \
	./kernel/kernel_cpu.h \
	./kernel/kernel_cpu.c \
	./util/num/num.h \
	./util/num/num.c \
	./util/timer/timer.h \
	./util/timer/timer.c
	$(CC_OMP) main.c -c

kernel/kernel_cpu.o: kernel/kernel_cpu.h kernel/kernel_cpu.c
	$(CC_OMP) kernel/kernel_cpu.c -c -o $@

util/num/num.o:	./util/num/num.h \
	./util/num/num.c
	$(CC_OMP)	./util/num/num.c -c -o $@

util/timer/timer.o:	./util/timer/timer.h \
	./util/timer/timer.c
	$(CC_OMP) ./util/timer/timer.c -c -o $@

# delete all object and executable files
clean:
	rm -f	*.o \
	./kernel/*.o \
	./util/num/*.o \
	./util/timer/*.o \
	lavaMD

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
