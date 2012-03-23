include ../../common/make.config
targets=OCL_particlefilter_naive OCL_particlefilter_float
ALL: $(targets)

OCL_particlefilter_naive: ex_particle_OCL_naive_seq.o
	$(OCL_LINK) -o $@ $^

OCL_particlefilter_float: ex_particle_OCL_float_seq.o
	$(OCL_LINK) -o $@ $^

%.o: %.[ch]
	$(OCL_OBJ) -o $@ $< -c

ex_particle_OCL_naive_seq.o: ex_particle_OCL_naive_seq.cpp
	$(OCL_OBJ) -o $@ $^ -c

ex_particle_OCL_float_seq.o: ex_particle_OCL_float_seq.cpp
	$(OCL_OBJ) -o $@ $^ -c

clean:
	rm -f *.o *~ OCL_particlefilter_naive OCL_particlefilter_float *.linkinfo

install: $(targets)
	install -m $(install_perms) -t $(OCL_BIN_DIR) $(targets)
