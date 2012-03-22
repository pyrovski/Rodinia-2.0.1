all: SRAD_V1 SRAD_V2

SRAD_V1:
	$(MAKE) -C srad_v1 install

SRAD_V2:
	$(MAKE) -C srad_v2 install

clean: SRAD_V1_clean SRAD_V2_clean

SRAD_V1_clean:
	$(MAKE) -C srad_v1 clean

SRAD_V2_clean:
	$(MAKE) -C srad_v2 clean

install: all
