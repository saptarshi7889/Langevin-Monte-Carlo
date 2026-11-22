// Include each model alone in other translation units and both together here.
#include <langevin/incommensurate.hpp>
#include <langevin/commensurate.hpp>
#include <langevin/incommensurate.hpp>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <type_traits>

int other_translation_unit_samples();

namespace {
constexpr int SPACE_SIZE = 3;
constexpr int TIME_SIZE = 4;
constexpr int BURN_IN = 2;
constexpr int SAMPLE_INTERVAL = 2;
constexpr int FIRST_SAMPLE = 4;
constexpr int SPLIT_STEP = 7;
constexpr int FINAL_STEP = 14;
constexpr int EXPECTED_SAMPLES = 6;
constexpr double COMMENSURATE_G = 0.3;
constexpr unsigned SEED = 73;

void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

std::string contents(const std::filesystem::path& path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open saved API checkpoint");
    return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

template<class Model>
void check_model(const std::filesystem::path& path, bool commensurate)
{
    langevin::Parameters parameters;
    parameters.L = SPACE_SIZE;
    parameters.Beta = TIME_SIZE;
    parameters.burn_in = BURN_IN;
    parameters.sample_every = SAMPLE_INTERVAL;
    parameters.seed = SEED;
    if (commensurate) parameters.g = COMMENSURATE_G;
    Model initial(parameters);
    require(initial.state().phi.isZero() && initial.state().step == 0, "Initial state is not zero");
    require(initial.coefficients().size() == TIME_SIZE - 1, "Missing kernel coefficients");
    require(!initial.step(), "Premature measurement during burn-in");
    initial.advance_to(FIRST_SAMPLE - 1);
    const auto first = initial.step();
    require(first && first->step == FIRST_SAMPLE &&
            first->time == FIRST_SAMPLE * parameters.dt, "Incorrect measurement schedule");
    initial.advance_to(SPLIT_STEP);
    initial.save_checkpoint(path);
    auto resumed = Model::from_checkpoint(path);
    resumed.advance_to(FINAL_STEP);
    Model uninterrupted(parameters);
    uninterrupted.advance_to(FINAL_STEP);
    require(resumed.state().phi == uninterrupted.state().phi, "API continuation changes the trajectory");
    require(resumed.state().samples == EXPECTED_SAMPLES, "API continuation loses sample counts");
    resumed.save_checkpoint(path.string() + ".resumed");
    uninterrupted.save_checkpoint(path.string() + ".full");
    require(contents(path.string() + ".resumed") == contents(path.string() + ".full"),
            "API continuation does not preserve complete state");
    resumed.advance_to(FINAL_STEP); // Equal targets are explicitly a no-op.
    bool rejected = false;
    try { resumed.advance_to(SPLIT_STEP); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "Backward integration was accepted");
    rejected = false;
    try {
        if (commensurate) langevin::incommensurate::Simulation::from_checkpoint(path);
        else langevin::commensurate::Simulation::from_checkpoint(path);
    } catch (const std::runtime_error&) { rejected = true; }
    require(rejected, "API loaded the wrong model's checkpoint");
    parameters.Beta = SPACE_SIZE; // Odd imaginary-time size is unsupported.
    rejected = false;
    try { Model invalid(parameters); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "API constructor bypasses parameter validation");
}
} // namespace

static_assert(!std::is_same_v<langevin::incommensurate::Simulation,
                             langevin::commensurate::Simulation>);
static_assert(std::is_const_v<std::remove_reference_t<decltype(
    std::declval<langevin::incommensurate::Simulation&>().state())>>);

int main(int argc, char* argv[])
{
    try {
        require(argc == 2, "Expected a temporary directory argument");
        const std::filesystem::path directory(argv[1]);
        check_model<langevin::incommensurate::Simulation>(directory / "api-incommensurate.txt", false);
        check_model<langevin::commensurate::Simulation>(directory / "api-commensurate.txt", true);
        require(other_translation_unit_samples() == SAMPLE_INTERVAL, "Multi-translation-unit model failed");
        std::cout << "Public header API, model coexistence, multi-TU linking and exact continuation passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
