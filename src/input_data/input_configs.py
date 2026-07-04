from dataclasses import dataclass, field
from typing import Any


@dataclass
class EnvironmentInputConfig:
    """Input contract for creating a local environment."""

    env_type: str | None = None
    source: str = "planetary_solver"
    data: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_raw(cls, raw_data: dict[str, Any]) -> "EnvironmentInputConfig":
        return cls(
            env_type=raw_data.get("environment"),
            source=raw_data.get("environment_source", "planetary_solver"),
            data=raw_data.get("environment_data") or {},
        )

    @property
    def is_explicit(self) -> bool:
        return self.source == "explicit"

    @property
    def uses_planetary_solver(self) -> bool:
        return self.source == "planetary_solver"


@dataclass
class ChemicalNetworkInputConfig:
    """Input contract for reaction-network setup."""

    network_type: str | None = None
    reaction_file: str | None = None
    data: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_raw(cls, raw_data: dict[str, Any]) -> "ChemicalNetworkInputConfig":
        data = raw_data.get("chemical_network") or {}
        return cls(
            network_type=data.get("type"),
            reaction_file=data.get("reaction_file"),
            data=data,
        )


@dataclass
class MoleculeInputConfig:
    """Input contract for molecular species / metabolite setup."""

    source: str | None = None
    molecule_type: str | None = None
    reaction_file: str | None = None
    data: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_raw(cls, raw_data: dict[str, Any]) -> "MoleculeInputConfig":
        data = raw_data.get("metabolites_data") or {}
        return cls(
            source=raw_data.get("molecule_source"),
            molecule_type=data.get("type"),
            reaction_file=data.get("reaction_file"),
            data=data,
        )


@dataclass
class EvolutionInputConfig:
    """Input contract for chemistry/evolution time integration."""

    data: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_raw(cls, raw_data: dict[str, Any]) -> "EvolutionInputConfig":
        return cls(data=raw_data.get("evol_params") or {})
