include ../../common/make.config

NVCC_FLAGS += -use_fast_math
CC_LINK_FLAGS += -lm

targets=particlefilter_naive particlefilter_float

all: $(targets)

particlefilter_naive: ex_particle_CUDA_naive_seq.cu
	$(NVCC_LINK) $^ -o $@

particlefilter_float: ex_particle_CUDA_float_seq.cu
	$(NVCC_LINK) $^ -o $@

clean:
	rm -f $(targets)

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
