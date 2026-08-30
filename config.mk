# ===================
#  Project paths
# ===================

ROOT := $(shell pwd)
CONDA_HOME ?= $(HOME)/miniforge3
CONDA ?= $(CONDA_HOME)/bin/conda
CONDA_ENV_NAME ?= lifeorig_env
CONDA_ENV_FILE := $(ROOT)/conda-environment.yml

# ===================
#  Python version
# ===================

PYTHON_VERSION ?= python=3.11
NUMPY_VERSION ?= "numpy>=2,<3"

# ===================
# Build mode
# ===================

BUILD_MODE ?= local

# ===================
# Logging
# ===================

LOG_LEVEL ?= INFO
COLORED_LOGGING ?= 1
LOGFILE ?= lifeorig.log

# ===================
# GPU section
# ===================

GPU_ACTIVE ?= 0
CUDA_VERSION ?= 12.8
ifeq ($(GPU_ACTIVE),1)
	CONDA_ENV_NAME := lifeorig_gpu_env
endif

# ===================
# MPI launcher
# ===================

MPI_LAUNCHER ?= mpirun
ifeq ($(BUILD_MODE),nersc)
	MPI_LAUNCHER := srun
endif

# ===================
#  unit tests
# ===================

NP_MAX ?= 1
UNIT_TEST_DIR := $(ROOT)/tests

# ===================
#  PHREEQC
# ===================

PHREEQC_DIR := $(ROOT)/external/phreeqc
PHREEQC_DB := $(PHREEQC_DIR)/database/phreeqc.dat

# ===================
#  Reaction MySQL DB
# ===================

LIFEORIG_REACTION_DB_HOST ?= localhost
LIFEORIG_REACTION_DB_NAME ?= lifeorig_reactions
LIFEORIG_REACTION_DB_USER ?= lifeorig
LIFEORIG_REACTION_DB_PASSWORD ?=
