import os

if os.path.exists(__file__.replace('.py', '_local.py')):
    # noinspection PyUnresolvedReferences
    from Figures.figures_config_local import *

#TODO fill missing values here with some default
