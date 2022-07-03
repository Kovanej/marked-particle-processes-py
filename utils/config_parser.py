
from typing import Dict

from utils.const import CONFIG_VALID_KEYS, CONFIG_NON_NULLABLE, CONFIG_OPTIONAL_VALUES


class ConfigParser(object):

    def __init__(self, config: Dict):
        self.config = config
        self._validate_the_config()
        self._parse_the_config()

    def _parse_the_config(self):
        self._parse_the_non_nullable_values()
        self._parse_the_nullable_values()

    def _parse_the_non_nullable_values(self):
        for _ in CONFIG_NON_NULLABLE:
            setattr(self, _, self.config[_])

    def _parse_the_nullable_values(self):
        for attr_name, attr_value in CONFIG_OPTIONAL_VALUES.items():
            if attr_name in self.config.keys():
                setattr(self, attr_name, self.config[attr_name])
            else:
                setattr(self, attr_name, attr_value)

    def _validate_the_config(self):
        _invalid_keys = [
            _ for _ in self.config.keys() if _ not in CONFIG_VALID_KEYS
        ]
        if len(_invalid_keys) > 0:
            self.valid_config = False
            raise ValueError(
                f"Unknown key(s) in config.json: {_invalid_keys}"
            )
        _missing = [
            _ for _ in CONFIG_NON_NULLABLE if (_ not in self.config.keys()) or (self.config[_] is None)
        ]
        if len(_missing) > 0:
            self.valid_config = False
            raise ValueError(
                f"Following non-nullable parameters either missing in config or have assigned null value: {_missing}"
            )
        self.valid_config = True

