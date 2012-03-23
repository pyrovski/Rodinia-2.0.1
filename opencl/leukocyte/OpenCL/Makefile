include ../../../common/make.config

# Matrix library path
MATRIX_DIR = ../meschach_lib

CC_FLAGS += -I$(MATRIX_DIR)
CC_LINK_FLAGS += -lm
targets=leukocyte

leukocyte: detect_main.o avilib.o find_ellipse.o find_ellipse_opencl.o track_ellipse.o track_ellipse_opencl.o misc_math.o OpenCL_helper_library.o $(MATRIX_DIR)/meschach.a
	$(OCL_LINK) $^ -o $@  

%.o: %.[ch]
	$(OCL_OBJ) $< -c

detect_main.o: detect_main.c find_ellipse.h track_ellipse.h avilib.h

find_ellipse.o: avilib.h

track_ellipse.o: find_ellipse.h avilib.h

$(MATRIX_DIR)/meschach.a:
	cd $(MATRIX_DIR); ./configure --with-all; make all; make clean

clean:
	rm -f *.o *~ leukocyte

install: $(targets)
	install -m $(install_perms) -t $(OCL_BIN_DIR) $(targets)
