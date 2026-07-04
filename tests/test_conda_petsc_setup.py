import importlib

import pytest


def test_petsc_driver_imports_when_petsc4py_is_available():
    pytest.importorskip("petsc4py")

    module = importlib.import_module("src.dynamics.petsc_driver")

    assert hasattr(module, "PETScDriver")


def test_mpi_and_petsc_are_available_in_conda_environment():
    mpi4py = pytest.importorskip("mpi4py")
    petsc4py = pytest.importorskip("petsc4py")

    assert mpi4py is not None
    assert petsc4py is not None
