include ../../common/make.config
PREFIX=${PARSECDIR}/pkgs/kernels/streamcluster/inst/${PARSECPLAT}

TARGET_C = sc_cpu

TARGET_O = sc_omp

targets=$(TARGET_C) $(TARGET_O)
ifdef version
  ifeq "$(version)" "parallel"
    CXXFLAGS :=	$(CXXFLAGS) -DENABLE_THREADS -pthread
  endif
endif

all: $(targets)

sc_cpu:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) streamcluster_original.cpp -o $(TARGET_C) $(LIBS) -DENABLE_THREADS -pthread

sc_omp:
	g++ -O3 -fopenmp -o $(TARGET_O) streamcluster_omp.cpp


clean:
	rm -f *.o *~ *.txt sc_cpu sc_omp

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
