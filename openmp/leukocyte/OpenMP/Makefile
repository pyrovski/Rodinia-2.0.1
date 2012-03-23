include ../../../common/make.config

CC_FLAGS += -I$(MATRIX_DIR)
CC_LINK_FLAGS += -lm
targets=leukocyte

# Matrix library path
MATRIX_DIR = ../meschach_lib

leukocyte: detect_main.o avilib.o find_ellipse.o track_ellipse.o misc_math.o $(MATRIX_DIR)/meschach.a
	$(CC_OMP_LINK) $^ -o $@

%.o: %.[ch]
	$(CC) $(CC_FLAGS) $< -c

detect_main.o: find_ellipse.h track_ellipse.h avilib.h
find_ellipse.o: avilib.h
track_ellipse.o: find_ellipse.h track_ellipse.h avilib.h

$(MATRIX_DIR)/meschach.a:
	cd $(MATRIX_DIR); ./configure --with-all; make all; make clean

clean:
	rm -f *.o leukocyte

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
