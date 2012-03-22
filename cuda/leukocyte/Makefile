include ../../common/make.config
targets=CUDA/leukocyte

all: CUDA/leukocyte 

CUDA/leukocyte:
	make -C CUDA

clean:
	make -C CUDA clean
	rm -f meschach_lib/meschach.a

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
