include ../../common/make.config
CC_LINK_FLAGS += -lm

targets=lavaMD

# link objects (binaries) together
lavaMD:	main.o \
	./kernel/kernel_cpu.o \
	./util/num/num.o \
	./util/timer/timer.o
	$(CC)	main.o \
	./kernel/kernel_cpu.o \
	./util/num/num.o \
	./util/timer/timer.o \
	$(CC_LINK_FLAGS) \
	-o lavaMD
# compile function files into objects (binaries)
main.o:	main.h \
	main.c \
	./kernel/kernel_cpu.h \
	./kernel/kernel_cpu.c \
	./util/num/num.h \
	./util/num/num.c \
	./util/timer/timer.h \
	./util/timer/timer.c
	$(CC)	main.c \
	-c \
	-o main.o \
	$(CC_FLAGS)

./kernel/kernel_cpu.o:	./kernel/kernel_cpu.h \
	./kernel/kernel_cpu.c
	$(CC)	./kernel/kernel_cpu.c \
	-c \
	-o ./kernel/kernel_cpu.o \
	$(CC_FLAGS)\
	$(OMP_FLAG)

./util/num/num.o:	./util/num/num.h \
	./util/num/num.c
	$(CC)	./util/num/num.c \
	-c \
	-o ./util/num/num.o \
	$(CC_FLAGS)

./util/timer/timer.o:	./util/timer/timer.h \
	./util/timer/timer.c
	$(CC)	./util/timer/timer.c \
	-c \
	-o ./util/timer/timer.o \
	$(CC_FLAGS)

# delete all object and executable files
clean:
	rm -f	*.o \
	./kernel/*.o \
	./util/num/*.o \
	./util/timer/*.o \
	lavaMD

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
