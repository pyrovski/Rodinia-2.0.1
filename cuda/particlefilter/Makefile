include ../../common/make.config

NVCC_FLAGS += -use_fast_math -lm

targets=particlefilter_naive particlefilter_float

all: $(targets)

particlefilter_naive: ex_particle_CUDA_naive_seq.cu
	$(CC) -I$(INCLUDE) -L$(CUDA_LIB_DIR) $^ -o $@

particlefilter_float: ex_particle_CUDA_float_seq.cu
	$(NVCC) -L$(CUDA_LIB_DIR) $^ -o $@

clean:
	rm -f particlefilter_naive particlefilter_float

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
