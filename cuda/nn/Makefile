include ../../common/make.config

CC_LINK_FLAGS += -lm

targets=nn hurricane_gen

all : $(targets)

clean :
	rm -rf *.o $(targets)

nn_cuda.cu.cpp: nn_cuda.cu
	$(NVCC_OBJ) -cuda $^ -o $@

nn: nn_cuda.cu.cpp
	$(NVCC_LINK) -o $@ $^

hurricane_gen : hurricane_gen.c
	$(NVCC_LINK) -o $@ $<

#data :
#	mkdir data
#	./gen_dataset.sh

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
