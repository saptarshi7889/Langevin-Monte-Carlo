#pragma once

#include "kernel.hpp"
#include <Eigen/Core>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace langevin {
inline constexpr int DEFAULT_SIZE = 16;
inline constexpr int DEFAULT_STEPS = 2000;
inline constexpr int DEFAULT_BURN_IN = 1000;
inline constexpr int DEFAULT_SAMPLE_EVERY = 20;
inline constexpr unsigned DEFAULT_SEED = 1234;
inline constexpr double DEFAULT_DT = 0.05;
inline constexpr double DEFAULT_U = 1.0;
inline constexpr double DEFAULT_K = 1.0;
inline constexpr double DEFAULT_ALPHA = 1.0;
inline constexpr double DEFAULT_S = 0.5;
inline constexpr double DEFAULT_G = 0.0;
inline constexpr double NOISE_VARIANCE = 2.0;
inline constexpr double LOCAL_FORCE_FACTOR = 2.0;
inline constexpr int CHECKPOINT_VERSION = 1;
inline constexpr int TEXT_PRECISION = std::numeric_limits<double>::max_digits10;

struct Parameters {
    int L = DEFAULT_SIZE, Beta = DEFAULT_SIZE;
    double u = DEFAULT_U, K = DEFAULT_K, alpha = DEFAULT_ALPHA, s = DEFAULT_S;
    double g = DEFAULT_G, dt = DEFAULT_DT;
    unsigned seed = DEFAULT_SEED;
    int burn_in = DEFAULT_BURN_IN, sample_every = DEFAULT_SAMPLE_EVERY;
};

inline void validate(const Parameters& p, bool commensurate) {
    if (p.L < 2 || p.Beta < 2 || p.Beta % 2 != 0 || p.burn_in < 0 ||
        p.sample_every <= 0 || !std::isfinite(p.u) || p.u <= 0 ||
        !std::isfinite(p.K) || p.K <= 0 || !std::isfinite(p.dt) || p.dt <= 0 ||
        !std::isfinite(p.alpha) || p.alpha < 0 || !std::isfinite(p.g) || p.g < 0 ||
        !std::isfinite(p.s) || p.s <= MIN_BATH_EXPONENT || p.s >= MAX_BATH_EXPONENT ||
        (!commensurate && p.g != 0)) {
        throw std::invalid_argument("Invalid parameters: require L >= 2, even Beta >= 2, positive u/K/dt, alpha/g >= 0, 0 < s < 2, burn_in >= 0, sample_every > 0; g requires the commensurate model");
    }
}

struct State {
    Parameters p;
    bool commensurate = false;
    int step = 0, samples = 0;
    double order_sum = 0.0;
    Eigen::MatrixXd phi;
    std::vector<double> rough_x, rough_tau;
    std::mt19937 generator;
    std::normal_distribution<double> normal{0.0, 1.0};
};

inline void write_parameters(std::ostream& out, const Parameters& p) {
    out << p.L << ' ' << p.Beta << ' ' << p.u << ' ' << p.K << ' ' << p.alpha << ' '
        << p.s << ' ' << p.g << ' ' << p.dt << ' ' << p.seed << ' ' << p.burn_in
        << ' ' << p.sample_every << '\n';
}

inline std::string parameter_key(const Parameters& p) {
    std::ostringstream out;
    out << std::setprecision(TEXT_PRECISION);
    write_parameters(out, p);
    return out.str();
}

inline std::ofstream output_file(const std::filesystem::path& path) {
    std::ofstream out(path);
    out.exceptions(std::ios::badbit | std::ios::failbit);
    out << std::setprecision(TEXT_PRECISION);
    return out;
}

inline void save_checkpoint(const std::filesystem::path& path, const State& state) {
    auto out = output_file(path);
    out << "LANGEVIN_CHECKPOINT " << CHECKPOINT_VERSION << '\n' << state.commensurate << '\n';
    write_parameters(out, state.p);
    out << state.step << ' ' << state.samples << ' ' << state.order_sum << '\n';
    for (int row = 0; row < state.p.Beta; ++row) {
        for (int col = 0; col < state.p.L; ++col) out << state.phi(row, col) << ' ';
        out << '\n';
    }
    for (double value : state.rough_x) out << value << ' ';
    out << '\n';
    for (double value : state.rough_tau) out << value << ' ';
    out << '\n' << state.generator << '\n' << state.normal << '\n';
    out.close();
}

inline State load_checkpoint(const std::string& path, bool commensurate) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open checkpoint: " + path);
    State state;
    std::string magic;
    int version = 0;
    input >> magic >> version >> state.commensurate;
    if (!input || magic != "LANGEVIN_CHECKPOINT" || version != CHECKPOINT_VERSION ||
        state.commensurate != commensurate) {
        throw std::runtime_error("Invalid checkpoint version or model");
    }
    auto& p = state.p;
    input >> p.L >> p.Beta >> p.u >> p.K >> p.alpha >> p.s >> p.g >> p.dt
          >> p.seed >> p.burn_in >> p.sample_every;
    if (!input) throw std::runtime_error("Incomplete checkpoint parameters");
    validate(p, commensurate);
    input >> state.step >> state.samples >> state.order_sum;
    if (!input || state.step < 0 || state.samples < 0 ||
        state.samples != (state.step > p.burn_in ? (state.step - p.burn_in) / p.sample_every : 0)) {
        throw std::runtime_error("Invalid checkpoint counters");
    }
    state.phi.resize(p.Beta, p.L);
    for (int row = 0; row < p.Beta; ++row)
        for (int col = 0; col < p.L; ++col) input >> state.phi(row, col);
    state.rough_x.resize(p.L - 1);
    state.rough_tau.resize(p.Beta - 1);
    for (double& value : state.rough_x) input >> value;
    for (double& value : state.rough_tau) input >> value;
    input >> state.generator >> state.normal;
    if (!input || !state.phi.allFinite() || !std::isfinite(state.order_sum))
        throw std::runtime_error("Incomplete or nonfinite checkpoint");
    for (const auto* values : {&state.rough_x, &state.rough_tau})
        for (double value : *values)
            if (!std::isfinite(value) || value < 0) throw std::runtime_error("Invalid checkpoint roughness");
    input >> std::ws;
    if (!input.eof()) throw std::runtime_error("Unexpected data after checkpoint");
    return state;
}

inline int integer(const std::string& text) {
    std::size_t end = 0;
    const int value = std::stoi(text, &end);
    if (end != text.size()) throw std::invalid_argument("Invalid integer: " + text);
    return value;
}

inline double real(const std::string& text) {
    std::size_t end = 0;
    const double value = std::stod(text, &end);
    if (end != text.size() || !std::isfinite(value)) throw std::invalid_argument("Invalid number: " + text);
    return value;
}

inline std::map<std::string, std::string> arguments(int argc, char** argv) {
    std::map<std::string, std::string> result;
    const std::string allowed = " L Beta u K alpha s g dt seed burn_in sample_every steps output checkpoint id ";
    for (int index = 1; index < argc; ++index) {
        const std::string arg(argv[index]);
        const auto equals = arg.find('=');
        if (equals == std::string::npos || equals == 0 || equals + 1 == arg.size())
            throw std::invalid_argument("Arguments must be name=value; use --help");
        const std::string name = arg.substr(0, equals);
        if (allowed.find(" " + name + " ") == std::string::npos ||
            !result.emplace(name, arg.substr(equals + 1)).second)
            throw std::invalid_argument("Unknown or duplicate argument: " + name);
    }
    return result;
}

inline void apply_parameters(Parameters& p, const std::map<std::string, std::string>& args) {
    for (const auto& entry : args) {
        const auto& name = entry.first;
        const auto& value = entry.second;
        if (name == "L") p.L = integer(value);
        else if (name == "Beta") p.Beta = integer(value);
        else if (name == "burn_in") p.burn_in = integer(value);
        else if (name == "sample_every") p.sample_every = integer(value);
        else if (name == "seed") {
            const int seed = integer(value);
            if (seed < 0) throw std::invalid_argument("seed must be nonnegative");
            p.seed = static_cast<unsigned>(seed);
        }
        else if (name == "u") p.u = real(value);
        else if (name == "K") p.K = real(value);
        else if (name == "alpha") p.alpha = real(value);
        else if (name == "s") p.s = real(value);
        else if (name == "g") p.g = real(value);
        else if (name == "dt") p.dt = real(value);
    }
}

template<class Force>
inline void advance(State& state, const std::vector<double>& coefficients, Force force) {
    const auto& p = state.p;
    Eigen::MatrixXd noise(p.Beta, p.L), spatial(p.Beta, p.L), temporal(p.Beta, p.L);
    Eigen::MatrixXd interaction(p.Beta, p.L);
    for (int row = 0; row < p.Beta; ++row)
        for (int col = 0; col < p.L; ++col) noise(row, col) = state.normal(state.generator);
    const auto drift = [&](const Eigen::MatrixXd& phi) -> Eigen::MatrixXd {
        Laplace_X(phi, spatial, p.Beta, p.L);
        Laplace_T(phi, temporal, p.Beta, p.L);
        force(phi, interaction, coefficients, p);
        return p.u / (PI * p.K) * spatial + 1.0 / (PI * p.u * p.K) * temporal + interaction;
    };
    const Eigen::MatrixXd first = drift(state.phi);
    const Eigen::MatrixXd increment = std::sqrt(NOISE_VARIANCE * p.dt) * noise;
    const Eigen::MatrixXd predicted = state.phi + p.dt * first + increment;
    const Eigen::MatrixXd second = drift(predicted);
    state.phi += HALF * p.dt * (first + second) + increment;
    if (!state.phi.allFinite()) throw std::runtime_error("Nonfinite field; check parameters and time step");
    ++state.step;
}

inline double measure(State& state) {
    // Only the incommensurate action permits arbitrary uniform field shifts.
    Eigen::MatrixXd observed = state.phi;
    if (!state.commensurate) observed.array() -= observed.mean();
    RoughnessX(observed, state.rough_x.data(), state.p.Beta, state.p.L);
    RoughnessT(observed, state.rough_tau.data(), state.p.Beta, state.p.L);
    double order = 0;
    OrderParameter(observed, order, state.p.Beta, state.p.L);
    order /= static_cast<double>(state.p.L) * state.p.Beta;
    state.order_sum += order;
    ++state.samples;
    return order;
}

template<class Force>
inline int run(int argc, char** argv, bool commensurate, bool continued, Force force) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--help") {
            std::cout << "Arguments (name=value): L Beta u K alpha s g dt seed burn_in sample_every steps output id checkpoint\n"
                      << "steps is the FINAL absolute step, time = step*dt; burn_in and sample_every are step counts.\n"
                      << "Initial defaults: L=16 Beta=16 u=1 K=1 alpha=1 s=0.5 g=0 dt=0.05 seed=1234\n"
                      << "steps=2000 burn_in=1000 sample_every=20 id=0 output=run-0\n"
                      << "Continuation requires checkpoint=PATH and restores all parameters and RNG state.\n";
            return 0;
        }
        const auto args = arguments(argc, argv);
        if (continued != (args.count("checkpoint") != 0))
            throw std::invalid_argument(continued ? "Continuation requires checkpoint=PATH" : "Use the continuation executable to load a checkpoint");
        State state;
        if (continued) state = load_checkpoint(args.at("checkpoint"), commensurate);
        else state.commensurate = commensurate;
        const auto original = parameter_key(state.p);
        apply_parameters(state.p, args);
        validate(state.p, commensurate);
        if (continued && original != parameter_key(state.p))
            throw std::invalid_argument("Cannot change physical, sampling or seed parameters of a checkpoint");
        const int steps = args.count("steps") ? integer(args.at("steps")) : DEFAULT_STEPS;
        const int id = args.count("id") ? integer(args.at("id")) : 0;
        if (id < 0 || steps <= state.step) throw std::invalid_argument("Require id >= 0 and steps > checkpoint step (or zero for an initial run)");
        if (!continued) {
            state.phi = Eigen::MatrixXd::Zero(state.p.Beta, state.p.L);
            state.rough_x.assign(state.p.L - 1, 0.0);
            state.rough_tau.assign(state.p.Beta - 1, 0.0);
            state.generator.seed(state.p.seed);
        }
        const std::filesystem::path output = args.count("output") ? args.at("output") : "run-" + std::to_string(id);
        if (!std::filesystem::create_directories(output))
            throw std::runtime_error("Output directory already exists; choose a new output=PATH");
        std::vector<double> coefficients(state.p.Beta - 1);
        Kernel(coefficients.data(), state.p.Beta, state.p.alpha, state.p.s);
        auto kernel_file = output_file(output / "kernel.csv");
        kernel_file << "separation,coefficient\n";
        for (int r = 1; r < state.p.Beta; ++r) kernel_file << r << ',' << coefficients[r - 1] << '\n';
        auto measurements = output_file(output / "measurements.csv");
        measurements << "step,time,order_parameter\n";
        while (state.step < steps) {
            advance(state, coefficients, force);
            if (state.step > state.p.burn_in && (state.step - state.p.burn_in) % state.p.sample_every == 0)
                measurements << state.step << ',' << state.step * state.p.dt << ',' << measure(state) << '\n';
        }
        save_checkpoint(output / "checkpoint.txt", state);
        auto configuration = output_file(output / "configuration.dat");
        configuration << state.phi << '\n';
        auto summary = output_file(output / "summary.txt");
        summary << "model " << (commensurate ? "commensurate" : "incommensurate") << '\n'
                << "id " << id << "\nstep " << state.step << "\ntime " << state.step * state.p.dt
                << "\nsamples " << state.samples << '\n';
        if (state.samples) {
            summary << "order_parameter " << state.order_sum / state.samples << '\n';
            auto roughness = output_file(output / "roughness.csv");
            roughness << "axis,separation,mean_squared_difference\n";
            for (int r = 1; r < state.p.L; ++r) roughness << "space," << r << ',' << state.rough_x[r - 1] / state.samples << '\n';
            for (int r = 1; r < state.p.Beta; ++r) roughness << "imaginary_time," << r << ',' << state.rough_tau[r - 1] / state.samples << '\n';
        } else summary << "order_parameter unavailable (no samples)\n";
        kernel_file.close();
        measurements.close();
        configuration.close();
        summary.close();
        std::cout << "Saved " << state.samples << " cumulative samples at step " << state.step << " to " << output << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
} // namespace langevin
