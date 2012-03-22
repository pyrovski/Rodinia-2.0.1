include ../common/make.config

OCL_DIRS  := backprop bfs cfd gaussian heartwall hotspot kmeans lavaMD leukocyte lud nn	nw srad streamcluster particlefilter pathfinder

all:
	mkdir -p $(OPENCL_BIN_DIR)
	$(MAKE) -C backprop install
	$(MAKE) -C bfs install
	$(MAKE) -C cfd install
	$(MAKE) -C hotspot install
	$(MAKE) -C kmeans install
	$(MAKE) -C lavaMD install
	$(MAKE) -C leukocyte install
	$(MAKE) -C lud/ocl install
	$(MAKE) -C nw install
	$(MAKE) -C srad install
	$(MAKE) -C streamcluster install
	$(MAKE) -C pathfinder install
	$(MAKE) -C particlefilter install
	$(MAKE) -C gaussian install
	$(MAKE) -C nn install
	$(MAKE) -C heartwall install

clean
	$(MAKE) -C backprop clean
	$(MAKE) -C bfs clean
	$(MAKE) -C cfd clean
	$(MAKE) -C hotspot clean
	$(MAKE) -C kmeans clean
	$(MAKE) -C lavaMD clean
	$(MAKE) -C leukocyte clean
	$(MAKE) -C lud/ocl clean
	$(MAKE) -C nw clean
	$(MAKE) -C srad clean
	$(MAKE) -C streamcluster clean
	$(MAKE) -C pathfinder clean
	$(MAKE) -C particlefilter clean
	$(MAKE) -C gaussian clean
	$(MAKE) -C nn clean
	$(MAKE) -C heartwall clean
