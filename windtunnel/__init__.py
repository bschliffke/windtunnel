# -*- coding: utf-8 -*-
"""Python package for basic boundary layer and concentration measurement analysis.
"""

import logging

from .utils import *
from .flow.stats import *
# from .concentration.stats import *
from .plots import *
from .timeseries import *
from .timeseries_nc import *
from .PuffConcentration import *
from .PointConcentration import *

# set standards for Logging module
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='windtunnel.log',
                    filemode='w')
# define a Handler which writes INFO messages or higher to the sys.stderr
console = logging.StreamHandler()
console.setLevel(logging.INFO)
# set a format which is simpler for console use
formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
# tell the handler to use this format
console.setFormatter(formatter)
# add the handler to the root logger
logger = logging.getLogger('')
logger.addHandler(console)


__all__ = [s for s in dir() if not s.startswith('_')]
