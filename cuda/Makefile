include ../common/make.config

CUDA_DIRS := backprop bfs cfd gaussian heartwall hotspot kmeans lavaMD leukocyte lud nn	nw srad streamcluster particlefilter pathfinder mummergpu

all: 
	mkdir -p $(CUDA_BIN_DIR)
	$(MAKE) -C backprop install
	$(MAKE) -C bfs install
	$(MAKE) -C cfd install
	$(MAKE) -C gaussian install
	$(MAKE) -C heartwall install
	$(MAKE) -C hotspot install
	$(MAKE) -C kmeans install
	$(MAKE) -C lavaMD install
	$(MAKE) -C leukocyte install
	$(MAKE) -C lud install
	$(MAKE) -C nn install
	$(MAKE) -C nw install
	$(MAKE) -C srad/srad_v1 install
	$(MAKE) -C srad/srad_v2 install
	$(MAKE) -C streamcluster install
	$(MAKE) -C particlefilter install
	$(MAKE) -C pathfinder install
	$(MAKE) -C mummergpu install


clean:
	$(MAKE) -C backprop clean
	$(MAKE) -C bfs clean
	$(MAKE) -C cfd clean
	$(MAKE) -C gaussian clean
	$(MAKE) -C heartwall clean
	$(MAKE) -C hotspot clean
	$(MAKE) -C kmeans clean
	$(MAKE) -C lavaMD clean
	$(MAKE) -C leukocyte clean
	$(MAKE) -C lud clean
	$(MAKE) -C nn clean
	$(MAKE) -C nw clean
	$(MAKE) -C srad/srad_v1 clean
	$(MAKE) -C srad/srad_v2 clean
	$(MAKE) -C streamcluster clean
	$(MAKE) -C particlefilter clean
	$(MAKE) -C pathfinder clean
	$(MAKE) -C mummergpu clean
