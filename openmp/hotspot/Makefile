include ../../common/make.config
CC=g++
targets=hotspot

hotspot: hotspot_openmp.cpp
	$(CC_OMP_LINK) -o $@ $^

clean:
	rm -f hotspot

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
