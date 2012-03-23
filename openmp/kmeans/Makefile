all: OPENMP/kmeans SERIAL/kmeans

OPENMP/kmeans:
	$(MAKE) -C kmeans_openmp install

SERIAL/kmeans:
	$(MAKE) -C kmeans_serial install

clean:
	$(MAKE) -C kmeans_serial clean
	$(MAKE) -C kmeans_openmp clean

install: all
