include ../../common/make.config

TARGET_C = sc_serial

TARGET_G = sc_openacc

TARGET_O = sc_openmp

targets=$(TARGET_C) $(TARGET_G) $(TARGET_O)

all: $(targets)

$(TARGET_C): streamcluster_original.cpp
	$(PGCXX_LINK) -Mnoopenmp $^ -o $@

$(TARGET_O): streamcluster_original.cpp
	$(PGCXX_LINK) -DOMP -mp -o $@ $^

$(TARGET_G): streamcluster_original.cpp
	$(OPENACC_CXX_LINK) -o $@ $^

clean:
	rm -f *.o *~  $(targets)

install: $(targets)
	install -m $(install_perms) -t $(OAC_BIN_DIR) $(targets)
