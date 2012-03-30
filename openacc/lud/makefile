include ../../common/make.config

EXECUTABLE      := lud_openacc
targets=$(EXECUTABLE)

SOURCES         := lud.c lud_openacc.c common/common.c 

OPENACC_FLAGS += -Icommon

BASENAMES       = $(basename $(SOURCES))

# ------------  generate the names of the object files  ------------------------
OBJECTS         = $(addsuffix .o,$(BASENAMES))

# ------------  make the executable (the default goal)  ------------------------
$(EXECUTABLE):	$(OBJECTS)
	$(OPENACC_LINK) $^ -o $@

# ------------  make the objects  ----------------------------------------------
%.o:%.c
	$(OPENACC_CC) -c $< -o $@

clean:
	-rm --force $(EXECUTABLE) $(OBJECTS) *~

.PHONY: clean

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
