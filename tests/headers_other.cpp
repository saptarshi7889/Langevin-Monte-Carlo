// This translation unit includes only the commensurate public header.
#include <langevin/commensurate.hpp>

int other_translation_unit_samples()
{
    constexpr int SIZE = 2;
    constexpr int INTERVAL = 2;
    constexpr int FINAL_STEP = 4;
    langevin::Parameters parameters;
    parameters.L = parameters.Beta = SIZE;
    parameters.burn_in = 0;
    parameters.sample_every = INTERVAL;
    langevin::commensurate::Simulation simulation(parameters);
    simulation.advance_to(FINAL_STEP);
    return simulation.state().samples;
}
