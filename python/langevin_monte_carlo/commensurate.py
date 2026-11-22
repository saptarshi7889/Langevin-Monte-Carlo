"""Initial and continued commensurate simulations using the C++ library."""

from ._api import _Simulation


class Simulation(_Simulation):
    """Zero-magnetization model with unshifted cos(2*phi) measurements."""

    __slots__ = ()
    _commensurate = True
