include ../common/make.config

OAC_DIRS := streamcluster lud mummergpu

all: $(OAC_BIN_DIR) $(OAC_DIRS)

$(OAC_BIN_DIR):
	mkdir -p $(OAC_BIN_DIR)

.PHONY: $(OAC_DIRS)

$(OAC_DIRS): 
	$(MAKE) -C $@ install

clean:
	@$(foreach dir, $(OAC_DIRS), $(MAKE) -C $(dir) clean;)
