include ../../common/make.config
targets=heartwall

CC_LINK_FLAGS += -lm

# link objects(binaries) together
heartwall: main.o ./AVI/avilib.o ./AVI/avimod.o
	$(CC_OMP_LINK) -o $@ $^

# compile main function file into object (binary)
main.o: main.c define.c kernel.c
	$(CC_OMP) main.c -I./AVI -c

./AVI/avilib.o ./AVI/avimod.o:
	$(MAKE) -C AVI

# delete all object files
clean:
	rm -f *.o AVI/*.o heartwall

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
