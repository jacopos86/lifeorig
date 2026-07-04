include config.mk

export ROOT
export GPU_ACTIVE
export PHREEQC_DATABASE := $(PHREEQC_DB)
export LOG_LEVEL
export COLORED_LOGGING
export LOGFILE

# ========================
# Base conda dependencies
# ========================

CONDA_BASE_DEPS = \
  $(PYTHON_VERSION) \
  $(NUMPY_VERSION) \
  scipy \
  colorlog \
  matplotlib \
  pytest \
  pint \
  periodictable \
  mpi4py \
  pip

# ========================
# PETSc dependencies
# ========================

ifeq ($(GPU_ACTIVE),1)
  PETSC_DEPS = \
    "petsc=*=*cuda*" \
    "petsc4py"
  GPU_DEPS = \
    cuda-version=$(CUDA_VERSION)
else
  PETSC_DEPS = \
    "petsc" \
    "petsc4py"
  GPU_DEPS =
endif

# ================================
# Generate conda-environment.yml
# ================================

environment:
	@echo "Generating $(CONDA_ENV_FILE)"
	@rm -f $(CONDA_ENV_FILE)
	@echo "name: $(CONDA_ENV_NAME)" >> $(CONDA_ENV_FILE)
	@echo "channels:" >> $(CONDA_ENV_FILE)
	@echo "  - conda-forge" >> $(CONDA_ENV_FILE)
	@echo "dependencies:" >> $(CONDA_ENV_FILE)
	@for pkg in $(CONDA_BASE_DEPS); do \
		echo "  - $$pkg" >> $(CONDA_ENV_FILE); \
	done
	@for pkg in $(PETSC_DEPS); do \
		echo "  - $$pkg" >> $(CONDA_ENV_FILE); \
	done
	@for pkg in $(GPU_DEPS); do \
		echo "  - $$pkg" >> $(CONDA_ENV_FILE); \
	done
	@echo "  - pip:" >> $(CONDA_ENV_FILE)
	@echo "      - -r $(ROOT)/requirements.txt" >> $(CONDA_ENV_FILE)
	@echo "      - -e $(ROOT)" >> $(CONDA_ENV_FILE)

$(CONDA_ENV_FILE): environment

# ===================
#  configuration
# ===================

configure : $(CONDA_ENV_FILE)
	@echo "Checking conda environment $(CONDA_ENV_NAME)..."
	@if ! $(CONDA) env list | grep -qw $(CONDA_ENV_NAME); then \
		echo "Creating conda environment $(CONDA_ENV_NAME) ..."; \
		$(CONDA) env create -f $(CONDA_ENV_FILE); \
	else \
		echo "Updating existing conda environment $(CONDA_ENV_NAME) ..."; \
		$(CONDA) env update -f $(CONDA_ENV_FILE) --prune; \
	fi

install :
	@echo "Installing package into conda environment $(CONDA_ENV_NAME)"
	@$(CONDA) run -n $(CONDA_ENV_NAME) python -m pip install -e .

build :
	@echo "No compiled LIFEORIG extension build is required yet."

test :
	@echo "Running tests in conda environment $(CONDA_ENV_NAME)"
	@$(CONDA) run -n $(CONDA_ENV_NAME) python -m pytest $(UNIT_TEST_DIR)

smoke-earth-volcanic :
	@PYTHON_BIN=python \
	OUTPUT_DIR="$(ROOT)/TEST_OUTPUT/Earth_volcanic_rock" \
	$(CONDA) run -n $(CONDA_ENV_NAME) bash $(ROOT)/TEST_SCRIPTS/run_preset_planet.sh

smoke-earth-hydro :
	@PYTHON_BIN=python \
	OUTPUT_DIR="$(ROOT)/TEST_OUTPUT/Earth_hydro_vent" \
	$(CONDA) run -n $(CONDA_ENV_NAME) bash $(ROOT)/TEST_SCRIPTS/run_hydro_vent_network.sh

install-phreeqc-db :
	mkdir -p $(PHREEQC_DIR)/database ; \
	wget -O $(PHREEQC_DIR)/database/phreeqc.dat https://raw.githubusercontent.com/usgs-coupled/phreeqc/master/database/phreeqc.dat ; \
	wget -O $(PHREEQC_DIR)/database/pitzer.dat https://raw.githubusercontent.com/usgs-coupled/phreeqc/master/database/pitzer.dat
test-ocean :
	PHREEQC_DATABASE=$(PHREEQC_DB) $(CONDA) run -n $(CONDA_ENV_NAME) python ./src/hydro_solver/test_ocean_chem_0.py
.PHONY :
	environment configure build install clean install-phreeqc-db test-ocean test smoke-earth-volcanic smoke-earth-hydro
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
	fi
	@if $(CONDA) env list | grep -qw $(CONDA_ENV_NAME); then \
		echo "Removing conda environment $(CONDA_ENV_NAME) ..."; \
		$(CONDA) env remove -n $(CONDA_ENV_NAME); \
	fi
	@if [ -f "$(CONDA_ENV_FILE)" ]; then \
		rm -f "$(CONDA_ENV_FILE)"; \
	fi
