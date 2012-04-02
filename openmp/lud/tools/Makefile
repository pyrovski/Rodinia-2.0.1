include ../../../common/make.config
gen_input: gen_input.c
	$(CC_OMP_LINK) -o $@ -O3 $^

clean: 
	rm -f gen_input *~
