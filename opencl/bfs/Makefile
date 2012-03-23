include ../../common/make.config

#C compiler
CC = $(CXX)

SRC = bfs.cpp timer.cc

EXE = bfs
targets=$(EXE)

$(EXE):$(SRC)
	$(OCL_LINK) $^ -o $@

errmsg:$(SRC)
	$(CC) $(CC_FLAGS) -lOpenCL -fopenmp $(SRC) -o $(EXE) -I$(OPENCL_INC) -L$(OPENCL_LIB) -D ERRMSG PTX_MSG

ptx:$(SRC)
	$(CC) $(CC_FLAGS) -lOpenCL -fopenmp $(SRC) -o $(EXE) -I$(OPENCL_INC) -L$(OPENCL_LIB) -D PTX_MSG

profile:$(SRC)
	$(CC) $(CC_FLAGS) -lOpenCL -fopenmp $(SRC) -o $(EXE) -I$(OPENCL_INC) -L$(OPENCL_LIB) -D PROFILING

res:$(SRC)
	$(CC) $(CC_FLAGS) -lOpenCL -fopenmp $(SRC) -o $(EXE) -I$(OPENCL_INC) -L$(OPENCL_LIB) -D RES_MSG

clean:
	rm -f $(EXE) $(EXE).linkinfo result*

install: $(targets)
	install -m $(install_perms) -t $(OCL_BIN_DIR) $(targets)
