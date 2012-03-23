include ../../../common/make.config
targets=kmeans

kmeans: cluster.o getopt.o kmeans.o kmeans_clustering.o 
	$(CC_OMP_LINK) -o $@ $^

%.o: %.[ch]
	$(CC_OMP) $< -c

kmeans_clustering.o: kmeans_clustering.c kmeans.h
	$(CC_OMP) $< -c kmeans_clustering.c -c

clean:
	rm -f *.o *~ kmeans 

install: $(targets)
	install -m $(install_perms) -t $(OMP_BIN_DIR) $(targets)
