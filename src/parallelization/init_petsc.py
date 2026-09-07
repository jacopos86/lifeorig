import os
import subprocess
from mpi4py import MPI
from src.parallelization.config import GPU_ACTIVE
from src.parallelization.mpi import mpi

# Central PETSc initialization point.
# All code should import PETSc from this module, not directly from
# petsc4py.

args = ["lifeorig"]

if GPU_ACTIVE:
    if "-use_gpu_aware_mpi" not in args:
        args += ["-use_gpu_aware_mpi", "0"]
    if "-mat_type" not in args:
        args += ["-mat_type", "aijcusparse"]
    if "-vec_type" not in args:
        args += ["-vec_type", "cuda"]

import petsc4py
petsc4py.init(args)
from petsc4py import PETSc

#
#   print PETSc info
#

def _print_petsc_info():
    opts = PETSc.Options()
    mat_type = opts.getString("mat_type", "default")
    vec_type = opts.getString("vec_type", "default")
    comm = PETSc.COMM_WORLD
    rank = comm.getRank()
    size = comm.getSize()
    if comm.getRank() == 0:
        PETSc.Sys.Print(
            "[LIFEORIG PETSc] "
            f"MPI ranks={size}, "
            f"GPU_ACTIVE={GPU_ACTIVE}, "
            f"mat_type={mat_type}, "
            f"vec_type={vec_type}"
        )
        if GPU_ACTIVE:
            PETSc.Sys.Print(
                "[LIFEORIG PETSc] PETSc GPU mode requested. "
                "Actual GPU use requires CUDA-enabled PETSc and "
                "Mat/Vec objects calling setFromOptions()."
            )
    if GPU_ACTIVE:
        visible = os.environ.get("CUDA_VISIBLE_DEVICES", "all")
        local_rank = os.environ.get(
            "OMPI_COMM_WORLD_LOCAL_RANK",
            os.environ.get("MPI_LOCALRANKID", "unknown")
        )
        device_report = comm.tompi4py().gather(
            (rank, local_rank, visible),
            root=0,
        )
        if rank == 0:
            PETSc.Sys.Print("[LIFEORIG PETSc] CUDA visibility by MPI rank:")
            for mpi_rank, mpi_local_rank, cuda_visible in device_report:
                PETSc.Sys.Print(
                    "[LIFEORIG PETSc] "
                    f"rank={mpi_rank}, "
                    f"local_rank={mpi_local_rank}, "
                    f"CUDA_VISIBLE_DEVICES={cuda_visible}"
                )
            try:
                nvidia_smi = subprocess.run(
                    ["nvidia-smi", "-L"],
                    check=False,
                    capture_output=True,
                    text=True,
                )
                if nvidia_smi.returncode == 0:
                    ngpu = sum(
                        line.startswith("GPU")
                        for line in nvidia_smi.stdout.splitlines()
                    )
                    PETSc.Sys.Print(
                        f"[LIFEORIG PETSc] nvidia-smi reports {ngpu} GPU(s)"
                    )
                    if visible == "all" and size > ngpu:
                        PETSc.Sys.Print(
                            "[LIFEORIG PETSc] WARNING: more MPI ranks than visible "
                            "GPUs. Multiple ranks may share a GPU."
                        )
            except FileNotFoundError:
                PETSc.Sys.Print(
                    "[LIFEORIG PETSc] nvidia-smi not found; cannot count GPUs."
                )

#
#   compare with mpi communicator
#

petsc_comm = PETSc.COMM_WORLD.tompi4py()
comparison = MPI.Comm.Compare(mpi.comm, petsc_comm)

if comparison not in (MPI.IDENT, MPI.CONGRUENT):
    raise RuntimeError(
        "mpi4py and PETSc communicators are not aligned"
    )

assert mpi.rank == PETSc.COMM_WORLD.getRank()
assert mpi.size == PETSc.COMM_WORLD.getSize()

#
#    print PETSc info
#

_print_petsc_info()
