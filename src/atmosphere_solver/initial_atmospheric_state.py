from dataclasses import replace
from pathlib import Path
import json
import numpy as np
from src.chemical_types.define_molecule_set import MolecularSpeciesSet
from src.global_parameters.global_variables import PROJECT_DIR
from src.utilities.logging_module import log
from src.common.units import Q_
from src.exo_chem.easy_chem_driver import run_easy_chem_full_profile
from src.common.phys_constants import G

#
#  State builder class
#

class InitialAtmosphericStateBuilder:
    def __init__(self, *, chem_input, planet_data, atmosph_data, units):
        self.chem_input = chem_input
        self.planet_data = planet_data
        self.atmosph_data = atmosph_data
        # internal units
        self.units = units
    def build(self, *, network=None, pressure_guess=None):
        settings = self.chem_input.initial_chemistry
        if not settings:
            log.error("initial_chemistry is required")
        method = settings.get("method")
        # different methods
        if method == "profile":
            z, temperature, mole_frac = self._from_profile(settings)
        else:
            z = self._build_altitude_grid()
            temperature = self._initial_temperature_profile(z)
            if method == "prescribed":
                mole_frac = {
                    name: np.full_like(z, float(value))
                    for name, value in settings["mole_fractions"].items()
                }
            elif method == "easychem":
                mole_frac = self._from_easychem(
                    temperature, pressure_guess
                )
            else:
                log.error(
                    f"Unknown initialization method: {method!r}"
                )
        # input validation
        self._validate_profile(z, temperature, mole_frac)
        # update chemical species
        # if network is present
        # otherwise the species considered are ones already available
        mole_frac, species = self._update_mole_fractions(
            mole_frac,
            network,
            n_layers=len(z)
        )
        return {
            "altitude": z,
            "temperature": temperature,
            "gravity": self._gravity_profile(z),
            "species": species,
            "mole_frac": mole_frac
        }
    # build altitude grid
    def _build_altitude_grid(self):
        n_layers = int(self.atmosph_data["n_layers"])
        z_max = self.atmosph_data["z_max"].to(
            self.units["length"]
        ).magnitude
        # Internal solver variables are unitless floats
        return np.linspace(0.0, z_max, n_layers)
    # gravity profile
    def _gravity_profile(self, z):
        mass = self.planet_data.planet_mass.to(
            self.units["mass"]
        ).magnitude
        radius = self.planet_data.planet_radius.to(
            self.units["length"]
        ).magnitude
        grav_const = G.to(
            f"{self.units['length']}^3 / {self.units['mass']} / {self.units['time']}^2"
        ).magnitude
        return grav_const * mass / (radius + z)**2
    # set initial temperature profile
    def _initial_temperature_profile(self, z):
        temperature_guess = self.atmosph_data.get("initial_temperature")
        if temperature_guess is None:
            temperature_guess = Q_(288.0, "K")
        value = temperature_guess.to(self.units["temperature"]).magnitude
        return np.full_like(z, float(value))
    # set mole fractions
    def _update_mole_fractions(self, mole_frac, network, n_layers):
        """Map initial profiles into the reaction network's species order."""
        network_species = (
            () if network is None else network.dyn_species
        )
        network_indices = (
            {} if network is None else network.species_to_index
        )
        # non reacting chemicals
        non_react_spec = [
            name for name in mole_frac
            if name not in network_indices
        ]
        # final set chemical species
        species, extra_spec_to_index = self._update_chem_species_set(
            network_species, 
            non_react_spec
        )
        profiles = np.zeros(
            (n_layers, len(species)),
            dtype=float
        )
        initialized = set()
        # run over moles
        for name, values in mole_frac.items():
            index = network_indices.get(name)
            if index is None:
                index = extra_spec_to_index.get(name)
                if index is None:
                    log.error(
                        f"Cannot resolve initial species {name!r}"
                    )
                log.info(
                    f"Additional atmospheric species: {name}; "
                    f"profile index={index}; absent from reaction network"
                )
            if index in initialized:
                log.error(
                    f"Multiple input names initialize network index {index}"
                )
            profiles[:, index] = values
            initialized.add(index)
        for index, entity in enumerate(network_species):
            if index not in initialized:
                log.info(
                    f"Initial profile set to zero: index={index}, {entity}"
                )
        return profiles, species
    def _update_chem_species_set(self, network_species, non_react_spec):
        # generate non reactive set
        non_react_set = MolecularSpeciesSet.from_species_names(non_react_spec)
        # Continue ID after existing network species.
        next_id = max(
            (entity.ID for entity in network_species),
            default=0
        ) + 1
        species = list(network_species)
        extra_indices = {}
        # iterate extra species
        for entity in non_react_set.molecules:
            entity = replace(entity, ID=next_id)
            next_id += 1
            index = len(species)
            species.append(entity)
            for name in entity.aliases:
                extra_indices[name] = index
            log.info(
                f"Additional atmospheric species: index={index}, {entity}"
            )
        return tuple(species), extra_indices
    def _from_profile(self, settings):
        path = Path(settings["file"])
        if not path.is_absolute():
            path = Path(PROJECT_DIR) / path
        # open file
        with path.open("r", encoding="utf-8") as stream:
            data = json.load(stream)
        # data
        grid = data["grid"]
        T = data["reference_atmosphere"]["temperature"]
        composition = data["composition"]
        if grid["coordinate"] != "altitude":
            log.error("Expected an altitude grid")
        if composition["basis"] != "mole_fraction":
            log.error("Expected mole-fraction profiles")
        # set altitude
        z = np.asarray(
            Q_(grid["values"], grid["units"]).to(
                self.units["length"]
            ).magnitude,
            dtype=float
        )
        # temperature
        temperature = np.asarray(
            Q_(T["values"], T["units"]).to(
                self.units["temperature"]
            ).magnitude,
            dtype=float
        )
        # moles
        mole_frac = {
            name: np.asarray(values, dtype=float)
            for name, values in composition["profiles"].items()
        }
        # Keep the source grid here. Regridding is a separate operation.
        return z, temperature, mole_frac
    def _from_easychem(self, temperature, pressure_guess):
        if pressure_guess is None:
            log.error(
                "EasyChem initialization requires a pressure guess"
            )
        pressure = np.asarray(pressure_guess, dtype=float)
        if not self.chem_input.atomic_abundances:
            log.error("EasyChem requires atomic_abundances")
        if not self.chem_input.chemical_species:
            log.error(
                "EasyChem requires an explicit compatible species list"
            )
        result = run_easy_chem_full_profile(
            pressure=Q_(pressure, self.units["pressure"]),
            temperature=Q_(
                temperature, self.units["temperature"]
            ),
            atomic_abund=self.chem_input.atomic_abundances,
            chemical_species=self.chem_input.chemical_species
        )
        return {
            name: np.asarray(values, dtype=float)
            for name, values in result.mole_fraction_profiles.items()
        }
    # validation
    def _validate_profile(self, z, temperature, mole_frac):
        if z.ndim != 1 or z.size < 2:
            log.error("Altitude must be a 1D profile")
        if not bool(np.all(np.isfinite(z)).item()):
            log.error("Altitude must be finite")
        if not bool(np.all(np.diff(z) > 0).item()):
            log.error("Altitude must increase strictly")
        if temperature.shape != z.shape:
            log.error("Temperature must match altitude")
        if not bool(
            np.all(np.isfinite(temperature) & (temperature > 0)).item()
        ):
            log.error("Temperature must be finite and positive")
        if not mole_frac:
            log.error("Initial composition is empty")
        total = np.zeros_like(z)
        for name, values in mole_frac.items():
            if values.shape != z.shape:
                log.error(f"Invalid profile shape for {name}")
            if not bool(
                np.all(np.isfinite(values) & (values >= 0)).item()
            ):
                log.error(f"Invalid mole fractions for {name}")
            total += values
        if not bool(
            np.all(np.abs(total - 1.0) <= 1e-6).item()
        ):
            log.error("Mole fractions must sum to one in every layer")