include ../../common/make.config

EXECUTABLE      := lud_openacc
targets=$(EXECUTABLE)

SOURCES         := lud.c lud_openacc.c common/common.c lud_base.c lud_lapack.c lud_cula.c

OPENACC_FLAGS += -Icommon -I/usr/local/atlas/include -I/usr/local/cula/include

BASENAMES       = $(basename $(SOURCES))

# ------------  generate the names of the object files  ------------------------
OBJECTS         = $(addsuffix .o,$(BASENAMES))

# ------------  make the executable (the default goal)  ------------------------
$(EXECUTABLE):	$(OBJECTS)
	$(OPENACC_LINK) $^ -o $@ -L/usr/local/cula/lib64 -lcula_core -lcula_lapack -L/usr/local/atlas/lib -llapack -lptf77blas -lptcblas -latlas 

# ------------  make the objects  ----------------------------------------------
%.o:%.c
	$(OPENACC_CC) -c $< -o $@

clean:
	-rm --force $(EXECUTABLE) $(OBJECTS) *~

.PHONY: clean

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
