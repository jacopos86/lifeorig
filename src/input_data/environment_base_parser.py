from abc import ABC, abstractmethod
from src.common.units import Q_
from src.environment.solvent import SolventData
from src.environment.external_drive_params import ExternalDriveParams
from src.utilities.logging_module import log

#
#    Local environment abstract class
#

class LocalEnvironmentParser(ABC):
    def __init__(self, env_data: dict | None, working_dir=None):
        self.env_data = env_data
        self.working_dir = working_dir
    # parse quantity
    def _quantity(self, input_dict: dict, required_keys=("units", "value"), desc="quantity"):
        if input_dict is None:
            return None
        missing = [key for key in required_keys if key not in input_dict]
        if missing:
            log.error(f"{desc} dictionary missing keys: {missing}")
        units = input_dict.get(required_keys[0])
        value = input_dict.get(required_keys[1])
        if units is None or value is None:
            log.error(f"{desc} dictionary keys cannot be None")
        return Q_(value, units)
    # abstract methods
    @abstractmethod
    def _build_local_env(self):
        raise NotImplementedError
    @abstractmethod
    def _get_solvent_data(self, env_data: dict) -> SolventData:
        raise NotImplementedError
    @abstractmethod
    def _parse_external_drive_params(self, data: dict, desc: str) -> ExternalDriveParams:
        raise NotImplementedError
    @abstractmethod
    def set_external_env_forces(self):
        raise NotImplementedError