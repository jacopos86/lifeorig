

# ============================================================
#  PETSc driver
# ============================================================

class PETScDriver:
    """
    PETSc-facing wrapper.

    For a pure particle code, PETSc is mainly coordinating:
    - MPI communicator
    - options database
    - future global diagnostics / load balance
    - optional TS-like stepping wrapper

    Later you can:
    - add DM for protocell ownership
    - add Vec for summary observables
    - couple to TS if hybrid field + particle model is used
    """
    def __init__(self, solver: DynamicsSolver, comm=PETSc.COMM_WORLD):
        self.solver = solver
        self.comm = comm
        self.rank = comm.getRank()
        self.size = comm.getSize()
    def solve(self, t0, tfinal, dt, drift_callback=None, n_species=None, log_every=10):
        nsteps = int(np.ceil((tfinal - t0) / dt))
        if self.rank == 0:
            PETSc.Sys.Print(
                f"Starting particle dynamics: t0={t0}, tfinal={tfinal}, dt={dt}, nsteps={nsteps}"
            )
        self.solver.setup()
        self.solver.run(
            nsteps=nsteps,
            dt=dt,
            drift_callback=drift_callback,
            n_species=n_species,
            log_every=log_every,
        )