from contextlib import contextmanager
import numpy as np
import pytest
from src.common.units import Q_
from src.parallelization.array_backend import (
    GPU_ACTIVE,
    cp,
    to_host,
    xp,
)
from src.parallelization.init_petsc import PETSc
from src.reactions.reaction_class import (
    ChemicalReaction,
    ColliderSet,
    ReactionRateParameters,
)
from src.reactions.reaction_rate_models import (
    DensityPowerLawRate,
    ModifiedArrheniusRate,
    TwoBodyRate,
)

pytest.importorskip("petsc4py")

@contextmanager
def vector_array(vector, shape, mode):
    if not GPU_ACTIVE:
        yield vector.getArray(readonly=mode == "r").reshape(shape)
        return

    handle = vector.getCUDAHandle(mode)
    try:
        size = int(np.prod(shape))
        memory = cp.cuda.UnownedMemory(
            handle,
            size * np.dtype(PETSc.ScalarType).itemsize,
            vector,
        )
        pointer = cp.cuda.MemoryPointer(memory, 0)
        yield cp.ndarray(
            shape,
            dtype=PETSc.ScalarType,
            memptr=pointer,
        )
    finally:
        vector.restoreCUDAHandle(handle, mode)


def test_petsc_evolves_implemented_rates_on_small_vertical_grid():
    """Evolve three independent gas reactions in three atmospheric layers."""
    species = ("H2O", "O1D", "OH", "H2CO", "H", "H2", "HCO")
    species_to_index = {name: index for index, name in enumerate(species)}
    nz = 3

    reactions = (
        ChemicalReaction(
            reaction_id="VPL0002",
            react_indices=xp.asarray([0, 1]),
            react_coeff=xp.asarray([1, 1]),
            prod_indices=xp.asarray([2]),
            prod_coeff=xp.asarray([2]),
            rate_model=TwoBodyRate(A=2.2e-10, tfac=0.0),
        ),
        ChemicalReaction(
            reaction_id="VPL0043",
            react_indices=xp.asarray([3, 4]),
            react_coeff=xp.asarray([1, 1]),
            prod_indices=xp.asarray([5, 6]),
            prod_coeff=xp.asarray([1, 1]),
            rate_model=ModifiedArrheniusRate(
                alpha=2.14e-12,
                beta=1.62,
                gamma=1090.0,
                T_ref=298.0,
            ),
        ),
        ChemicalReaction(
            reaction_id="VPL0044",
            react_indices=xp.asarray([4]),
            react_coeff=xp.asarray([2]),
            prod_indices=xp.asarray([5]),
            prod_coeff=xp.asarray([1]),
            rate_model=DensityPowerLawRate(
                A=8.85e-33,
                T_ref=287.0,
                temperature_power=-0.6,
            ),
            collider_set=ColliderSet(
                species_indices=xp.arange(len(species)),
                efficiencies=xp.ones(len(species)),
            ),
        ),
    )

    stoichiometry = xp.zeros((len(species), len(reactions)))
    for column, reaction in zip(stoichiometry.T, reactions):
        reaction.add_stoichiometry(column)

    initial = xp.zeros((nz, len(species)))
    initial[:, species_to_index["H2O"]] = xp.asarray(
        [2.0e8, 1.5e8, 1.0e8]
    )
    initial[:, species_to_index["O1D"]] = xp.asarray(
        [1.0e5, 2.0e5, 3.0e5]
    )
    initial[:, species_to_index["H2CO"]] = xp.asarray(
        [2.0e7, 2.0e7, 2.0e7]
    )
    initial[:, species_to_index["H"]] = xp.asarray(
        [1.0e6, 1.0e6, 1.0e6]
    )

    temperature = Q_(xp.asarray([220.0, 260.0, 300.0]), "kelvin")

    def rhs(_ts, _time, state_vector, derivative_vector):
        with vector_array(state_vector, initial.shape, "r") as state_array:
            context = ReactionRateParameters(
                temperature=temperature,
                species_to_index=species_to_index,
                species_number_densities=Q_(
                    state_array,
                    "1 / centimeter**3",
                ),
            )
            rates = xp.stack(
                [reaction.rate(state_array, context) for reaction in reactions],
                axis=-1,
            )
            derivative = rates @ stoichiometry.T
            with vector_array(
                derivative_vector,
                initial.shape,
                "w",
            ) as derivative_array:
                derivative_array[:] = derivative

    if GPU_ACTIVE:
        state = PETSc.Vec().create(comm=PETSc.COMM_SELF)
        state.setSizes(initial.size)
        state.setType(PETSc.Vec.Type.CUDA)
        state.setUp()
        with vector_array(state, initial.shape, "w") as state_array:
            state_array[:] = initial
    else:
        state = PETSc.Vec().createSeq(initial.size, comm=PETSc.COMM_SELF)
        state.getArray()[:] = initial.ravel()

    ts = PETSc.TS().create(comm=PETSc.COMM_SELF)
    ts.setType(PETSc.TS.Type.BDF)
    ts.setRHSFunction(rhs)
    ts.setSolution(state)
    ts.setTime(0.0)
    ts.setTimeStep(1.0e-3)
    ts.setMaxTime(2.0e4)
    ts.setMaxSteps(10_000)
    ts.setExactFinalTime(PETSc.TS.ExactFinalTime.MATCHSTEP)
    ts.setTolerances(rtol=1.0e-7, atol=1.0e-3)
    ts.setFromOptions()
    ts.solve(state)

    with vector_array(state, initial.shape, "r") as final_array:
        final = to_host(final_array).copy()
    initial_host = to_host(initial)

    assert ts.getConvergedReason() > 0
    assert ts.getStepNumber() < 10_000
    assert np.all(final >= -1.0e-6)
    assert np.all(final[:, species_to_index["O1D"]] < initial_host[:, 1])
    assert np.all(final[:, species_to_index["OH"]] > initial_host[:, 2])
    assert np.all(final[:, species_to_index["H2CO"]] < initial_host[:, 3])
    assert np.all(final[:, species_to_index["H2"]] > initial_host[:, 5])

    oxygen_family_initial = (
        initial_host[:, 0] + initial_host[:, 1] + initial_host[:, 2]
    )
    oxygen_family_final = final[:, 0] + final[:, 1] + final[:, 2]
    np.testing.assert_allclose(
        oxygen_family_final,
        oxygen_family_initial,
        rtol=2.0e-7,
    )
