VENV = env
PYTHON = $(VENV)/bin/python3
PIP = $(VENV)/bin/pip

configure : requirements.txt
	python3 -m venv $(VENV); \
	. $(VENV)/bin/activate; \
	$(PIP) install -r requirements.txt
install :
	. $(VENV)/bin/activate ; \
	$(PIP) install .
.PHONY :
	clean
clean :
	rm -rf ./src/*~ ; \
	if [ -d ./src/__pycache__ ] ; \
	then \
		rm -rf ./src/__pycache__ ; \
	fi ; \
	if [ -d ./build ] ; \
	then \
		rm -rf ./build ; \
	fi ; \
	if [ -d ./__pycache__ ] ; \
	then \
		rm -rf ./__pycache__ ; \
	fi ; \
	if [ -d ./lifeorig.egg-info ] ; \
	then \
		rm -rf ./lifeorig.egg-info ; \
	fi ; \
	if [ -d $(VENV) ] ; \
	then \
		rm -rf $(VENV) ; \
	fi ;
