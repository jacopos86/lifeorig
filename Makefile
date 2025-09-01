ROOT = $(shell pwd)
VENV = $(ROOT)/env
PYTHON = $(VENV)/bin/python3
PIP = $(VENV)/bin/pip

configure : $(ROOT)/requirements.txt
	python3 -m venv $(VENV); \
	. $(VENV)/bin/activate; \
	$(PIP) install -r $(ROOT)/requirements.txt
install :
	. $(VENV)/bin/activate ; \
	$(PIP) install .
.PHONY :
	clean
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
