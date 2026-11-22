"""Initial and continued incommensurate simulations using the C++ library."""

from ._api import _Simulation


class Simulation(_Simulation):
    """Finite-magnetization model with centered cos(phi) measurements."""

    __slots__ = ()
    _commensurate = False
