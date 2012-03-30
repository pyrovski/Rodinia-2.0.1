include ../../../common/make.config

EXECUTABLE      := lud_omp
targets=$(EXECUTABLE)

SOURCES         := lud.c lud_omp.c ../common/common.c 

CC_FLAGS += -I../common

C_SOURCES       = $(filter     %.c, $(SOURCES))
BASENAMES       = $(basename $(SOURCES))

# ------------  generate the names of the object files  ------------------------
OBJECTS         = $(addsuffix .o,$(BASENAMES))

# ------------  make the executable (the default goal)  ------------------------
$(EXECUTABLE):	$(OBJECTS)
	$(CC_OMP_LINK) $^ -o $@

# ------------  make the objects  ----------------------------------------------
%.o:%.c
	$(CC_OMP) $(CC_FLAGS) -c $< -o $@

clean:
	-rm --force $(EXECUTABLE) $(OBJECTS) *~

.PHONY: clean

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
