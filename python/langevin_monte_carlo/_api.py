"""Shared, NumPy-independent Python interface to the C++ simulation libraries."""

from dataclasses import asdict, dataclass
import os

from . import _core


@dataclass(frozen=True)
class Parameters:
    """Physical parameters and the integer-step sampling schedule."""

    L: int = 16
    Beta: int = 16
    u: float = 1.0
    K: float = 1.0
    alpha: float = 1.0
    s: float = 0.5
    g: float = 0.0
    dt: float = 0.05
    seed: int = 1234
    burn_in: int = 1000
    sample_every: int = 20


@dataclass(frozen=True)
class Measurement:
    """One scheduled order-parameter observation."""

    step: int
    time: float
    order_parameter: float


@dataclass(frozen=True)
class State:
    """An owned, read-only snapshot; arrays support the Python buffer protocol."""

    parameters: Parameters
    commensurate: bool
    step: int
    samples: int
    order_sum: float
    phi: memoryview
    rough_x: memoryview
    rough_tau: memoryview

    @property
    def time(self) -> float:
        return self.step * self.parameters.dt

    @property
    def mean_order_parameter(self) -> float | None:
        return self.order_sum / self.samples if self.samples else None


def _view(data: bytes, shape=None) -> memoryview:
    return memoryview(data).cast("d", shape=shape) if shape else memoryview(data).cast("d")


class _Simulation:
    __slots__ = ("_handle",)
    _commensurate: bool

    def __init__(self, parameters: Parameters | None = None):
        if parameters is None:
            parameters = Parameters()
        if not isinstance(parameters, Parameters):
            raise TypeError("parameters must be a Parameters instance")
        self._handle = _core.create(self._commensurate, asdict(parameters))

    @classmethod
    def from_checkpoint(cls, path):
        """Restore the same C++ checkpoint format used by the CLI executables."""
        simulation = cls.__new__(cls)
        simulation._handle = _core.restore(cls._commensurate, os.fsencode(path))
        return simulation

    def advance_to(self, final_step: int) -> None:
        """Advance to an absolute step and accumulate scheduled measurements."""
        _core.advance_to(self._handle, final_step)

    def step(self) -> Measurement | None:
        """Advance once; return None unless this is a scheduled measurement."""
        observation = _core.step(self._handle)
        return Measurement(*observation) if observation is not None else None

    def save_checkpoint(self, path) -> None:
        """Write or replace a complete checkpoint; its parent must exist."""
        _core.save(self._handle, os.fsencode(path))

    @property
    def state(self) -> State:
        """Copy the complete visible state without exposing mutable C++ memory."""
        snapshot = _core.snapshot(self._handle)
        parameters = Parameters(**snapshot.pop("parameters"))
        snapshot["phi"] = _view(snapshot["phi"], (parameters.Beta, parameters.L))
        snapshot["rough_x"] = _view(snapshot["rough_x"])
        snapshot["rough_tau"] = _view(snapshot["rough_tau"])
        return State(parameters=parameters, **snapshot)

    @property
    def coefficients(self) -> memoryview:
        """An owned read-only snapshot, including the alpha/pi**2 prefactor."""
        return _view(_core.coefficients(self._handle))

    @property
    def current_step(self) -> int:
        return _core.progress(self._handle)[0]

    @property
    def samples(self) -> int:
        return _core.progress(self._handle)[1]

    @property
    def time(self) -> float:
        return _core.progress(self._handle)[2]
