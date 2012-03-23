include ../../../common/make.config

EXE = lud
targets=$(EXE)
CC_FLAGS += -I../common

$(EXE): lud.cpp ../common/common.c
	$(OCL_LINK) -o $@ $^ 

clean:
	rm -f $(EXE)

install: $(targets)
	install -m $(install_perms) -t $(OCL_BIN_DIR) $(targets)
