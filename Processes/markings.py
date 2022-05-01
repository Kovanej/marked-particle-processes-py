
from typing import Union, List, Optional
import logging
from datetime import datetime


class Mark(object):

    def __init__(
            self,
            mark_type: str,  # discrete or continuous
            mark_value: Optional[Union[float, int]],
            number_of_levels: Optional[int] = None  # only relevant for discrete mark_type
    ):
        logging.info(f"{datetime.now()}: Oh Hi Mark!")
        self.mark_type = mark_type
        self.number_of_levels = number_of_levels
        self.mark_value = mark_value
