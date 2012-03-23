include ../../common/make.config
CC_FLAGS += -Wall
CC_LINK_FLAGS += -lm
targets=nn

all : nn

clean :
	rm -rf *.o nn

nn : nn_openmp.c
	$(CC_OMP_LINK) -o $@ $<

hurricane_gen : hurricane_gen.c
	$(CC_OMP_LINK) -o $@ $<

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
