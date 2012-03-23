include ../../common/make.config
targets=srad

srad: $(wildcard *.cpp)
	$(OCL_LINK) $^ -o $@

clean:
	@rm -f srad *.o *.out *~ *\#

install: $(targets)
	install -m $(install_perms) -t $(OCL_BIN_DIR) $(targets)
