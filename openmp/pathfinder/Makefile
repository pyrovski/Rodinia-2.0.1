include ../../common/make.config
CC = $(CXX)
SRC = pathfinder.cpp
EXE = pathfinder
targets=$(EXE)

$(EXE): $(SRC)
	$(CC_OMP_LINK) $^ -o $@

clean:
	rm -f pathfinder

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
