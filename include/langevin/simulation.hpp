#pragma once

#include "kernel.hpp"
#include "numerics.hpp"
#include <Eigen/Core>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <optional>
#include <utility>

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

namespace detail {
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
    OrderParameter(observed, order, state.p.Beta, state.p.L, state.commensurate);
    order /= static_cast<double>(state.p.L) * state.p.Beta;
    state.order_sum += order;
    ++state.samples;
    return order;
}

template<bool Commensurate>
struct Interaction {
    void operator()(const Eigen::MatrixXd& phi, Eigen::MatrixXd& result,
                    const std::vector<double>& coefficients, const Parameters& p) const {
        if constexpr (Commensurate) {
            Diss_comm(phi, result, coefficients.data(), p.Beta, p.L,
                      LOCAL_FORCE_FACTOR * p.g / (PI * PI));
        } else {
            Diss(phi, result, coefficients.data(), p.Beta, p.L);
        }
    }
};
} // namespace detail

struct Measurement {
    int step;
    double time;
    double order_parameter;
};

namespace detail {
// The public model headers give each specialization a model-specific name.
template<bool Commensurate>
class Simulation {
public:
    explicit Simulation(const Parameters& parameters = Parameters{}) {
        validate(parameters, Commensurate);
        state_.p = parameters;
        state_.commensurate = Commensurate;
        state_.phi = Eigen::MatrixXd::Zero(parameters.Beta, parameters.L);
        state_.rough_x.assign(parameters.L - 1, 0.0);
        state_.rough_tau.assign(parameters.Beta - 1, 0.0);
        state_.generator.seed(parameters.seed);
        initialize_kernel();
    }

    static Simulation from_checkpoint(const std::filesystem::path& path) {
        return Simulation(load_checkpoint(path.string(), Commensurate));
    }

    // State is read-only so callers cannot accidentally invalidate the kernel,
    // sampling counters or reproducible random state.
    const State& state() const noexcept { return state_; }
    const std::vector<double>& coefficients() const noexcept { return coefficients_; }

    // Advance exactly one step; return a measurement only at scheduled steps.
    std::optional<Measurement> step() {
        if (state_.step == std::numeric_limits<int>::max()) {
            throw std::overflow_error("Langevin step counter exhausted");
        }
        advance(state_, coefficients_, Interaction<Commensurate>{});
        const auto& p = state_.p;
        if (state_.step > p.burn_in && (state_.step - p.burn_in) % p.sample_every == 0) {
            const double order = measure(state_);
            return Measurement{state_.step, state_.step * p.dt, order};
        }
        return std::nullopt;
    }

    // The target is an absolute step, not an additional number of steps.
    // Equal targets are a no-op; running backward is an error.
    void advance_to(int final_step) {
        if (final_step < state_.step) {
            throw std::invalid_argument("Cannot advance to a step earlier than the current state");
        }
        while (state_.step < final_step) step();
    }

    void save_checkpoint(const std::filesystem::path& path) const {
        langevin::save_checkpoint(path, state_);
    }

private:
    State state_;
    std::vector<double> coefficients_;

    explicit Simulation(State state) : state_(std::move(state)) { initialize_kernel(); }

    void initialize_kernel() {
        coefficients_.resize(state_.p.Beta - 1);
        kernel(coefficients_.data(), state_.p.Beta, state_.p.alpha, state_.p.s);
    }
};
} // namespace detail
} // namespace langevin
