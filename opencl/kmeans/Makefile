include ../../common/make.config
targets=kmeans

kmeans: cluster.o getopt.o read_input.o kmeans_clustering.o rmse.o kmeans.o
	$(OCL_LINK) $^ -o $@

%.o: %.c
	$(OCL_OBJ) $< -c

%.o: %.cpp
	$(OCL_OBJ) $< -c

clean:
	rm -f *.o *~ kmeans *.linkinfo

install: $(targets)
	install -m $(install_perms) -t $(OCL_BIN_DIR) $(targets)
