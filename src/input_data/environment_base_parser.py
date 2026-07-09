from abc import ABC, abstractmethod
from src.environment.external_drive_params import LiquidLevelParams
from src.environment.solvent import SolventData
from src.utilities.logging_module import log

#
#    Local environment abstract class
#

class LocalEnvironmentParser(ABC):
    def __init__(self, input_parser, env_data: dict | None):
        self.input_parser = input_parser
        self.env_data = env_data
        self.working_dir = getattr(input_parser, "working_dir", None)
    # parse Q from dict
    def _parse_quantity_from_dict(self, *args, **kwargs):
        return self.input_parser._parse_quantity_from_dict(*args, **kwargs)
    # abstract methods
    @abstractmethod
    def _build_local_env(self):
        raise NotImplementedError
    @abstractmethod
    def _get_solvent_data(self, env_data: dict) -> SolventData:
        raise NotImplementedError
    @abstractmethod
    def set_external_env_forces(self):
        raise NotImplementedError
    def _get_liquid_level_params(self, env_data: dict) -> LiquidLevelParams:
        solvent_data = env_data.get("solvent_data", {})
        if "liquid_level_params" in solvent_data:
            data = solvent_data.get("liquid_level_params")
        elif "liquid_level_params" in env_data:
            data = env_data.get("liquid_level_params")
        else:
            data = env_data.get("water_level_params")
        if data is None:
            log.error("Missing 'liquid_level_params' in input")
        try:
            model_type = data.get("model_type", data.get("type")).lower()
        except KeyError as e:
            log.error(f"Missing key in liquid_level_params: {e}")
        except AttributeError:
            log.error("Missing key in liquid_level_params: model_type")
        base_level = self._parse_quantity_from_dict(
            input_dict=data.get("base_level"),
            required_keys=("units", "value"),
            desc="liquid base level"
        )
        amplitude = self._parse_quantity_from_dict(
            input_dict=data.get("amplitude"),
            required_keys=("units", "value"),
            desc="liquid amplitude level"
        )
        period = self._parse_quantity_from_dict(
            input_dict=data.get("period"),
            required_keys=("units", "value"),
            desc="liquid oscillation period"
        )
        phase = float(data.get("phase", 0.0))
        switch_times = self._parse_quantity_from_dict(
            input_dict=data.get("switch_times"),
            required_keys=("units", "value"),
            desc="liquid level switch times"
        )
        levels = self._parse_quantity_from_dict(
            input_dict=data.get("levels"),
            required_keys=("units", "value"),
            desc="piecewise water levels"
        )
        return LiquidLevelParams(
            model_type=model_type,
            base_level=base_level,
            amplitude=amplitude,
            period=period,
            phase=phase,
            switch_times=switch_times,
            levels=levels,
        )
