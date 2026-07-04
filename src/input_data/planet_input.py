from src.planet_params.earth_params import get_earth_planetary_params, get_earth_stellar_params
from src.planet_params.europa_params import get_europa_planetary_params, get_europa_stellar_params
from src.planet_params.mars_params import get_mars_planetary_params, get_mars_stellar_params
from src.planet_params.titan_params import get_titan_planetary_params, get_titan_stellar_params
from src.planet_params.venus_params import get_venus_planetary_params, get_venus_stellar_params
from src.planet_params.planetary_params import PlanetaryEnvironmentParams
from src.stellar_params.stellar_data import StellarParams

#
#   Builder planetary input
#

class PlanetInputBuilder:
    #
    #    build planetary parameters from input data
    #
    def build_planetary_params(
        self,
        planet_model: str,
        planetary_data: dict | None,
    ) -> PlanetaryEnvironmentParams:
        if planet_model == "Earth":
            return self._apply_planetary_overrides(get_earth_planetary_params(), planetary_data)
        if planet_model == "Europa":
            return self._apply_planetary_overrides(get_europa_planetary_params(), planetary_data)
        if planet_model == "Mars":
            return self._apply_planetary_overrides(get_mars_planetary_params(), planetary_data)
        if planet_model == "Titan":
            return self._apply_planetary_overrides(get_titan_planetary_params(), planetary_data)
        if planet_model == "Venus":
            return self._apply_planetary_overrides(get_venus_planetary_params(), planetary_data)
        if planetary_data is None:
            return None
        # return planetary data
        return PlanetaryEnvironmentParams(
            name=planetary_data.get("name", planet_model or "custom"),
            planet_radius=self._parse_quantity_from_dict(
                input_dict=planetary_data.get("planet_radius"),
                required_keys=("units", "value"),
                desc="planet radius"
            ),
            planet_mass=self._parse_quantity_from_dict(
                input_dict=planetary_data.get("planet_mass"),
                required_keys=("units", "value"),
                desc="planet mass"
            ),
            orbital_distance=self._parse_quantity_from_dict(
                input_dict=planetary_data.get("orbital_distance"),
                required_keys=("units", "value"),
                desc="orbital distance"
            ),
            rotation_period=self._parse_quantity_from_dict(
                input_dict=planetary_data.get("rotation_period"),
                required_keys=("units", "value"),
                desc="planet rotation period",
            ),
            obliquity=float(planetary_data.get("obliquity", 0.0)),
            eccentricity=float(planetary_data.get("eccentricity", 0.0)),
            tidal_locked=bool(planetary_data.get("tidal_locked", False)),
            day_night_contrast=float(planetary_data.get("day_night_contrast", 0.0)),
            chemical_env=planetary_data.get("exo_chemistry"),
            atmosphere=self._build_atmosphere_params(planetary_data.get("atmosphere"))
        )
    # apply JSON overrides on top of preset planetary parameters
    def _apply_planetary_overrides(
        self,
        preset_data: PlanetaryEnvironmentParams,
        planetary_data: dict | None,
    ) -> PlanetaryEnvironmentParams:
        if planetary_data is None:
            return preset_data
        atmosphere_override = self._build_atmosphere_params(planetary_data.get("atmosphere"))
        if atmosphere_override is not None:
            preset_data.atmosphere = self._merge_atmosphere_data(
                base=preset_data.atmosphere or {},
                override=atmosphere_override,
            )
        return preset_data
    # merge atmosphere data, preserving nested hydrostatic solver defaults
    @staticmethod
    def _merge_atmosphere_data(base: dict, override: dict) -> dict:
        merged = {
            **base,
            **{key: value for key, value in override.items() if value is not None and key != "hydrostatic_solver"},
        }
        override_solver = override.get("hydrostatic_solver")
        if override_solver is not None:
            base_solver = base.get("hydrostatic_solver") or {}
            merged["hydrostatic_solver"] = {
                **base_solver,
                **{key: value for key, value in override_solver.items() if value is not None and key != "settings"},
                "settings": {
                    **(base_solver.get("settings") or {}),
                    **{
                        key: value
                        for key, value in (override_solver.get("settings") or {}).items()
                        if value is not None
                    },
                },
            }
        return merged
    # build atmospheric data
    def _build_atmosphere_params(self, atmosphere_data: dict | None) -> dict | None:
        if atmosphere_data is None:
            return None
        return {
            "n_layers": atmosphere_data.get("n_layers"),
            "z_max": self._parse_quantity_from_dict(
                input_dict=atmosphere_data.get("z_max"),
                required_keys=("units", "value"),
                desc="atmosphere z_max"
            ),
            "top_pressure": self._parse_quantity_from_dict(
                input_dict=atmosphere_data.get("top_pressure"),
                required_keys=("units", "value"),
                desc="atmosphere top_pressure"
            ),
            "atmosphere_mass_fraction": atmosphere_data.get("atmosphere_mass_fraction"),
            "hydrostatic_solver": self._build_hydrostatic_solver_params(atmosphere_data),
        }
    # build hydrostatic solver data
    def _build_hydrostatic_solver_params(self, atmosphere_data: dict) -> dict | None:
        hydrostatic_data = atmosphere_data.get("hydrostatic_solver")
        if hydrostatic_data is None:
            legacy_setting_keys = {
                "max_iter_loop",
                "rel_tol",
                "logp_tol",
                "abs_tol",
                "damping_loop",
                "anderson_depth",
            }
            if not any(key in atmosphere_data for key in legacy_setting_keys):
                return None
            hydrostatic_data = {"type": "standard"}
        if isinstance(hydrostatic_data, str):
            hydrostatic_data = {"type": hydrostatic_data}
        settings_data = hydrostatic_data.get("settings", {})
        return {
            "type": hydrostatic_data.get("type", "standard"),
            "settings": {
                "max_iter": settings_data.get("max_iter", atmosphere_data.get("max_iter_loop")),
                "rel_tol": settings_data.get("rel_tol", atmosphere_data.get("rel_tol")),
                "logp_tol": settings_data.get("logp_tol", atmosphere_data.get("logp_tol")),
                "abs_tol": self._parse_optional_quantity(
                    settings_data.get("abs_tol", atmosphere_data.get("abs_tol")),
                    desc="hydrostatic solver abs_tol",
                ),
                "damping": settings_data.get("damping", atmosphere_data.get("damping_loop")),
                "min_pressure": self._parse_optional_quantity(
                    settings_data.get("min_pressure"),
                    desc="hydrostatic solver min_pressure",
                ),
                "anderson_depth": settings_data.get(
                    "anderson_depth",
                    atmosphere_data.get("anderson_depth"),
                ),
            },
        }
    # parse an optional quantity, accepting already-built Quantity objects
    def _parse_optional_quantity(self, value, desc="quantity"):
        if value is None or hasattr(value, "to"):
            return value
        return self._parse_quantity_from_dict(
            input_dict=value,
            required_keys=("units", "value"),
            desc=desc,
        )
    #
    #    build stellar parameters from input data
    #
    def build_stellar_params(
        self,
        planet_model: str,
        stellar_data: dict | None,
    ) -> StellarParams | None:
        if planet_model == "Earth":
            return get_earth_stellar_params()
        if planet_model == "Europa":
            return get_europa_stellar_params()
        if planet_model == "Mars":
            return get_mars_stellar_params()
        if planet_model == "Titan":
            return get_titan_stellar_params()
        if planet_model == "Venus":
            return get_venus_stellar_params()
        if stellar_data is None:
            return None
        return StellarParams(
            name=stellar_data.get("name", "star"),
            spectral_class=stellar_data.get("spectral_class"),
            effective_temperature=self._parse_quantity_from_dict(
                input_dict=stellar_data.get("star_temperature"),
                required_keys=("units", "value"),
                desc="star temperature"
            ),
            radius=self._parse_quantity_from_dict(
                input_dict=stellar_data.get("star_radius"),
                required_keys=("units", "value"),
                desc="star radius"
            ),
            mass=self._parse_quantity_from_dict(
                input_dict=stellar_data.get("star_mass"),
                required_keys=("units", "value"),
                desc="star mass"
            ),
            luminosity=self._parse_quantity_from_dict(
                input_dict=stellar_data.get("luminosity"),
                required_keys=("units", "value"),
                desc="star luminosity"
            ),
        )
