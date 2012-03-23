include ../../common/make.config
targets=bfs

CC=g++

bfs: bfs.cpp
	$(CC_OMP_LINK) $^ -o $@

clean:
	rm -f $(targets) result.txt

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
