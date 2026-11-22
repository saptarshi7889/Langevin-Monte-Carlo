#pragma once

#include "simulation.hpp"

namespace langevin::commensurate {
// A zero-field initial simulation and its checkpoint continuation share one API.
using Simulation = detail::Simulation<true>;
} // namespace langevin::commensurate
