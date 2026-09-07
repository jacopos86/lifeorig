import numpy as np
from src.parallelization.config import GPU_ACTIVE
from src.parallelization.gpu_runtime import initialize

"""Process-wide NumPy/CuPy array backend."""

runtime = initialize()

if GPU_ACTIVE:
    try:
        import cupy as cp
    except ImportError as exc:
        raise RuntimeError(
            "GPU_ACTIVE=1, but CuPy is not available in the environment"
        ) from exc
    xp = cp
else:
    cp = None
    xp = np

#  convert array

def asarray(array, dtype=None):
    """Convert an array to the active backend."""
    return xp.asarray(array, dtype=dtype)

#  get numpy array representation

def to_host(array):
    """Return a NumPy representation of a backend array."""
    if GPU_ACTIVE:
        return cp.asnumpy(array)
    return np.asarray(array)

#  get scalar object

def scalar(value):
    """Convert a NumPy/CuPy scalar or size-one array to a Python scalar."""
    return value.item()

# is GPU array

def is_gpu_array(array):
    """Return whether an object is a CuPy array."""
    return GPU_ACTIVE and isinstance(array, cp.ndarray)

# synchronize GPU work

def synchronize():
    """Wait for outstanding GPU work; no-op on CPU."""
    if GPU_ACTIVE:
        cp.cuda.get_current_stream().synchronize()