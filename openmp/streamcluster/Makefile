include ../../common/make.config

TARGET_C = sc_cpu

TARGET_O = sc_omp

targets=$(TARGET_C) $(TARGET_O)

#CXXFLAGS += -DENABLE_THREADS

all: $(targets)

$(TARGET_C): streamcluster_original.cpp
	$(CXX_LINK) $^ -o $@ -DENABLE_THREADS -lpthread

$(TARGET_O): streamcluster_omp.cpp
	$(CXX_OMP_LINK) -o $@ $^


clean:
	rm -f *.o *~ *.txt $(targets)

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
