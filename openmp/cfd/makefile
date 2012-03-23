include ../../common/make.config

CC=$(CXX)

CC_FLAGS += -Dblock_length=$(OMP_NUM_THREADS)

targets=euler3d_cpu euler3d_cpu_double pre_euler3d_cpu pre_euler3d_cpu_double
all: $(targets)

euler3d_cpu: euler3d_cpu.cpp
	$(CC_OMP_LINK) $^ -o $@

euler3d_cpu_double: euler3d_cpu_double.cpp
	$(CC_OMP_LINK) $^ -o $@

pre_euler3d_cpu: pre_euler3d_cpu.cpp
	$(CC_OMP_LINK) $^ -o $@

pre_euler3d_cpu_double: pre_euler3d_cpu_double.cpp
	$(CC_OMP_LINK) $^ -o $@

clean:
	rm -f euler3d euler3d_cpu euler3d_double euler3d_cpu_double pre_euler3d pre_euler3d_cpu pre_euler3d_double pre_euler3d_cpu_double

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
