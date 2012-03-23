include ../../common/make.config
CC =$(CXX)

SRC = euler3d.cpp

EXE = euler3d
targets=$(EXE)

$(EXE):$(SRC)
	$(OCL_LINK) $^ -o $@

errmsg:$(SRC)
	$(CC) -lOpenCL $(SRC) -o $(EXE) -I$(OPENCL_INC) -L$(OPENCL_LIB) -D ERRMSG $(FLAGS)

ptx:$(SRC)
	$(CC) -lOpenCL $(SRC) -o $(EXE) -I$(OPENCL_INC) -L$(OPENCL_LIB) -D PTX_MSG $(FLAGS)

res:$(SRC)
	$(CC) -lOpenCL $(SRC) -o $(EXE) -I$(OPENCL_INC) -L$(OPENCL_LIB) -D RES_MSG $(FLAGS)

clean: $(SRC)
	rm -f $(EXE) $(EXE).linkinfo result*

install: $(targets)
	install -m $(install_perms) -t $(OCL_BIN_DIR) $(targets)
