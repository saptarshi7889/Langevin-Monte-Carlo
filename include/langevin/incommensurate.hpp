#pragma once

#include "simulation.hpp"

namespace langevin::incommensurate {
// A zero-field initial simulation and its checkpoint continuation share one API.
using Simulation = detail::Simulation<false>;
} // namespace langevin::incommensurate
