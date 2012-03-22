include ../../../common/make.config
CC_LINK_FLAGS += -lm

targets=srad

# link objects(binaries) together
srad:	main.o
	$(CC) $(CC_FLAGS) $(OMP_FLAG) main.o $(CC_LINK_FLAGS) -o $@

# compile main function file into object (binary)
main.o: 	main.c \
	define.c \
	graphics.c
	$(CC) $(CC_FLAGS) $(OMP_FLAG) main.c -c

# delete all object files
clean:
	rm -f *.o srad

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
