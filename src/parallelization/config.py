import os
from pathlib import Path

GPU_ACTIVE = (
    os.environ.get("GPU_ACTIVE", "0") == "1"
)

MPI_ROOT = 0

CUDA_SOURCE_DIR = (
    Path(__file__).parent / "cuda"
)