# -*- coding: utf-8 -*-
# @Author: p-chambers
# @Date:   2016-11-17 17:32:25
# @Last Modified by:   p-chambers
# @Last Modified time: 2016-11-18 15:12:31
__all__ = ["flight_tools", "global_tools", "interpolate", "sim_manager",
            "simulator", "weather"]

from . import simulator
from . import weather


import logging
try:  # Python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass
logging.basicConfig(level=logging.DEBUG)
logging.getLogger(__name__).addHandler(NullHandler())