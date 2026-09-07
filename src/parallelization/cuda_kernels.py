from functools import lru_cache
from pathlib import Path
from src.parallelization.config import (
    CUDA_SOURCE_DIR,
    GPU_ACTIVE
)

"""Compile and access custom CUDA kernels through CuPy."""

if GPU_ACTIVE:
    from src.parallelization.array_backend import cp
else:
    cp = None

#  require GPU

def _require_gpu():
    if not GPU_ACTIVE:
        raise RuntimeError(
            "Custom CUDA kernels require GPU_ACTIVE=1"
        )

# source path

def _source_path(source_file):
    path = Path(CUDA_SOURCE_DIR) / source_file
    if not path.is_file():
        raise FileNotFoundError(
            f"CUDA source file not found: {path}"
        )
    return path

#  load module

@lru_cache(maxsize=None)
def load_module(source_file, kernel_names):
    """
    Compile one CUDA source file and cache the resulting CuPy module.

    Parameters
    ----------
    source_file
        Filename relative to CUDA_SOURCE_DIR.
    kernel_names
        Tuple containing the kernel entry-point names.
    """
    _require_gpu()

    source_path = _source_path(source_file)
    names = tuple(kernel_names)

    return cp.RawModule(
        code=source_path.read_text(),
        backend="nvrtc",
        options=(
            "--std=c++11",
            f"-I{CUDA_SOURCE_DIR}",
        ),
        name_expressions=names
    )

#  get kernel

def get_kernel(source_file, kernel_name, kernel_names=None):
    """Return one compiled CUDA kernel."""
    if kernel_names is None:
        kernel_names = (kernel_name,)
    else:
        kernel_names = tuple(kernel_names)
    module = load_module(source_file, kernel_names)
    return module.get_function(kernel_name)

#  launch cupy kernel

def launch_kernel(
    kernel,
    work_size,
    args,
    shared_mem=0,
    stream=None
):
    """Launch one CUDA thread for each work item."""
    _require_gpu()

    work_size = int(work_size)
    if work_size <= 0:
        raise ValueError("work_size must be positive")

    threads_per_block = min(256, work_size)
    number_of_blocks = (
        work_size + threads_per_block - 1
    ) // threads_per_block
    
    grid = (number_of_blocks,)
    block = (threads_per_block,)

    kernel(
        grid,
        block,
        args,
        shared_mem=shared_mem,
        stream=stream
    )