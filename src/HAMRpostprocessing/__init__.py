# enables to directly use:
# from HAMRpostprocessing import Snapshot
from .snapshot import Snapshot
from .analysis import calc_Mtot

# lists everything that get's imported with:
# from HAMRpostprocessing import *
__all__ = ["Snapshot", "calc_Mtot"]