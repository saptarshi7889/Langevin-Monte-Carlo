#include <langevin/incommensurate.hpp>
#include <langevin/commensurate.hpp>
#include <iostream>

namespace {
constexpr int SPACE_SIZE = 3;
constexpr int TIME_SIZE = 4;
constexpr int BURN_IN_STEPS = 2;
constexpr int SAMPLE_INTERVAL = 2;
constexpr int CHECKPOINT_STEP = 7;
constexpr int FINAL_STEP = 14;
constexpr double COMMENSURATE_G = 0.2;
}

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: library-example NEW_CHECKPOINT_PATH\n";
        return 1;
    }
    try {
        langevin::Parameters parameters;
        parameters.L = SPACE_SIZE;
        parameters.Beta = TIME_SIZE;
        parameters.burn_in = BURN_IN_STEPS;
        parameters.sample_every = SAMPLE_INTERVAL;

        langevin::incommensurate::Simulation initial(parameters);
        initial.advance_to(CHECKPOINT_STEP);
        initial.save_checkpoint(argv[1]);

        auto continued = langevin::incommensurate::Simulation::from_checkpoint(argv[1]);
        continued.advance_to(FINAL_STEP);
        std::cout << "Incommensurate samples: " << continued.state().samples << '\n';
        std::cout << "Field shape: " << continued.state().phi.rows() << " x "
                  << continued.state().phi.cols() << '\n';

        parameters.g = COMMENSURATE_G;
        langevin::commensurate::Simulation commensurate(parameters);
        while (commensurate.state().step < FINAL_STEP) {
            if (const auto sample = commensurate.step()) {
                std::cout << "Commensurate order at t=" << sample->time << ": "
                          << sample->order_parameter << '\n';
            }
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
