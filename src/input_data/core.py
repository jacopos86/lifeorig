from abc import ABC, abstractmethod
import json
from src.common.units import Q_
from src.utilities.logging_module import log

#
#   pure abstract input class
#

class AbstractInput(ABC):
    _data = None
    """Pure abstract base for all input classes."""
    # load json file
    def _load_json(self, json_file: str) -> dict:
        try:
            f = open(json_file)
        except:
            msg = "file: " + json_file + " not found"
            raise Exception(msg)
        data = json.load(f)
        f.close()
        return data
    # read yaml
    def read_input_json(self, json_file: str):
        """Read JSON file and populate parameters fields."""
        # load JSON
        self._data = self._load_json(json_file)
        # parse common data
        self._parse_data()
        # validate directories
        self._validate()
    # parsing method
    @abstractmethod
    def _parse_data(self):
        """Read input data """
        raise NotImplementedError
    # validation method
    @abstractmethod
    def _validate(self):
        raise NotImplementedError
    @staticmethod
    def _parse_quantity_from_dict(input_dict: dict, required_keys=("units", "value"), desc="quantity") -> Q_:
        """
        Extract a Quantity from a dictionary, checking required keys.
        Works for both scalar and array values.
        Args:
            input_dict: dict containing quantity info (keys: units, value or values)
            required_keys: keys to check in the dict
            desc: description for error messages
        Returns:
            Q\_: instance with the proper unit
        """
        if input_dict is None:
            return None
        # check required keys
        missing = [k for k in required_keys if k not in input_dict]
        if missing:
            log.error(f"{desc} dictionary missing keys: {missing}")
        # extract value(s)
        val = input_dict.get(required_keys[1])
        units = input_dict.get(required_keys[0])
        if units is None or val is None:
            log.error(f"{desc} dictionary keys cannot be None")
        # convert to Quantity
        try:
            q = Q_(val, units)  # Q_ handles scalar or array
        except Exception as e:
            log.error(f"Failed to convert {desc} to Quantity: {e}")
        return q