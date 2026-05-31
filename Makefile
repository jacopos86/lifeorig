ROOT = $(shell pwd)
VENV = $(ROOT)/env
PYTHON = $(VENV)/bin/python3
PIP = $(VENV)/bin/pip

PHREEQC_DIR = $(ROOT)/external/phreeqc
PHREEQC_DB = $(PHREEQC_DIR)/database/phreeqc.dat

configure : $(ROOT)/requirements.txt
	python3 -m venv $(VENV); \
	. $(VENV)/bin/activate; \
	$(PIP) install -r $(ROOT)/requirements.txt
install :
	. $(VENV)/bin/activate ; \
	$(PIP) install .
install-phreeqc-db :
	mkdir -p $(PHREEQC_DIR)/database ; \
	wget -O $(PHREEQC_DIR)/database/phreeqc.dat https://raw.githubusercontent.com/usgs-coupled/phreeqc/master/database/phreeqc.dat ; \
	wget -O $(PHREEQC_DIR)/database/pitzer.dat https://raw.githubusercontent.com/usgs-coupled/phreeqc/master/database/pitzer.dat
test-ocean :
	PHREEQC_DATABASE=$(PHREEQC_DB) $(PYTHON) ./src/ocean_chem/test_ocean_chem_0.py
.PHONY :
	clean install-phreeqc-db test-ocean
clean :
	rm -rf $(ROOT)/src/*~ ; \
	if [ -d $(ROOT)/src/__pycache__ ] ; \
	then \
		rm -rf $(ROOT)/src/__pycache__ ; \
	fi ; \
	if [ -d $(ROOT)/build ] ; \
	then \
		rm -rf $(ROOT)/build ; \
	fi ; \
	if [ -d $(ROOT)/__pycache__ ] ; \
	then \
		rm -rf $(ROOT)/__pycache__ ; \
	fi ; \
	if [ -d $(ROOT)/lifeorig.egg-info ] ; \
	then \
		rm -rf $(ROOT)/lifeorig.egg-info ; \
	fi ; \
	if [ -d $(VENV) ] ; \
	then \
		rm -rf $(VENV) ; \
	fi ;
