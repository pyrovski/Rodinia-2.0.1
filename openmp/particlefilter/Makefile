include ../../common/make.config
targets=particle_filter

CC_FLAGS += -ffast-math
CC_LINK_FLAGS += -lm

particle_filter: ex_particle_OPENMP_seq.c
	$(CC_OMP_LINK) $^ -o $@

clean:
	rm -f particle_filter

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
