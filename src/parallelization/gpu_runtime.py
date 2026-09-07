"""CUDA device selection."""

import os
from dataclasses import dataclass
from src.parallelization.config import GPU_ACTIVE
from src.parallelization.mpi import mpi

@dataclass(frozen=True)
class GPURuntime:
    enabled: bool
    device_id: int | None
    device_count: int

#
#   set up GPU environment
#

def initialize() -> GPURuntime:
    if not GPU_ACTIVE:
        return GPURuntime(
            enabled=False,
            device_id=None,
            device_count=0
        )
    # if GPU active
    try:
        import cupy as cp
    except ImportError as exc:
        raise RuntimeError(
            "GPU_ACTIVE=1, but no cuda devices are visible"
        ) from exc
    # number devices
    try:
        device_count = cp.cuda.runtime.getDeviceCount()
    except cp.cuda.runtime.CUDARuntimeError as exc:
        raise RuntimeError(
            "GPU_ACTIVE=1, but CUDA initialization failed"
        ) from exc
    if device_count == 0:
        raise RuntimeError(
            "GPU_ACTIVE=1, but no CUDA devices are visible"
        )
    explicit_device = os.environ.get("GPU_DEVICE_ID")
    # explicit device
    if explicit_device is not None:
        device_id = int(explicit_device)
    else:
        local_rank = int(
            os.environ.get(
                "OMPI_COMM_WORLD_LOCAL_RANK",
                os.environ.get(
                    "SLURM_LOCALID",
                    os.environ.get(
                        "MPI_LOCALRANKID",
                        mpi.rank,
                    ),
                ),
            )
        )
        device_id = local_rank % device_count
    if not 0 <= device_id < device_count:
        raise RuntimeError(
            f"GPU_DEVICE_ID={device_id} is invalid; "
            f"{device_count} CUDA device(s) are visible"
        )
    cp.cuda.Device(device_id).use()
    return GPURuntime(
        enabled=True,
        device_id=device_id,
        device_count=device_count
    )