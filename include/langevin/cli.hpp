#pragma once

#include "simulation.hpp"
#include <iostream>
#include <map>

namespace langevin::cli {
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

template<bool Commensurate>
inline int run(int argc, char** argv, bool continued) {
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
        using Simulation = detail::Simulation<Commensurate>;
        std::optional<Simulation> restored;
        Parameters parameters;
        if (continued) {
            restored.emplace(Simulation::from_checkpoint(args.at("checkpoint")));
            parameters = restored->state().p;
        }
        const auto original = parameter_key(parameters);
        apply_parameters(parameters, args);
        validate(parameters, Commensurate);
        if (continued && original != parameter_key(parameters))
            throw std::invalid_argument("Cannot change physical, sampling or seed parameters of a checkpoint");
        const int steps = args.count("steps") ? integer(args.at("steps")) : DEFAULT_STEPS;
        const int id = args.count("id") ? integer(args.at("id")) : 0;
        if (id < 0 || steps <= (continued ? restored->state().step : 0)) throw std::invalid_argument("Require id >= 0 and steps > checkpoint step (or zero for an initial run)");
        Simulation simulation = continued ? std::move(*restored) : Simulation(parameters);
        const auto& state = simulation.state();
        const std::filesystem::path output = args.count("output") ? args.at("output") : "run-" + std::to_string(id);
        if (!std::filesystem::create_directories(output))
            throw std::runtime_error("Output directory already exists; choose a new output=PATH");
        const auto& coefficients = simulation.coefficients();
        auto kernel_file = output_file(output / "kernel.csv");
        kernel_file << "separation,coefficient\n";
        for (int r = 1; r < state.p.Beta; ++r) kernel_file << r << ',' << coefficients[r - 1] << '\n';
        auto measurements = output_file(output / "measurements.csv");
        measurements << "step,time,order_parameter\n";
        while (state.step < steps) {
            if (const auto observation = simulation.step()) {
                measurements << observation->step << ',' << observation->time << ','
                             << observation->order_parameter << '\n';
            }
        }
        simulation.save_checkpoint(output / "checkpoint.txt");
        auto configuration = output_file(output / "configuration.dat");
        configuration << state.phi << '\n';
        auto summary = output_file(output / "summary.txt");
        summary << "model " << (Commensurate ? "commensurate" : "incommensurate") << '\n'
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
} // namespace langevin::cli
