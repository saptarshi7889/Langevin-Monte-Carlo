#define main simulation_entry_point
#ifdef COMMENSURATE
#include "../src/commensurate/Langevin_simulation_initial_commensurate.cpp"
#else
#include "../src/incommensurate/Langevin_simulation_initial_simulation.cpp"
#endif
#undef main

#include <algorithm>
#include <stdexcept>

namespace {
constexpr int ROWS = 4;
constexpr int COLS = 3;
constexpr double TOLERANCE = 1e-10;
constexpr double DIFFERENCE_STEP = 1e-6;
constexpr double FORCE_TOLERANCE = 1e-8;
constexpr double TEST_ALPHA = 0.7;
constexpr double TEST_LOCAL_FORCE = 0.3;
constexpr int HARMONIC_STEPS = 100000;
constexpr int HARMONIC_BURN_IN = 10000;
constexpr int HARMONIC_STRIDE = 10;
constexpr double HARMONIC_DT = 0.05;
constexpr double VARIANCE_TOLERANCE = 0.02;

void require(bool condition, const char* description) {
    if (!condition) throw std::runtime_error(description);
}

void near(double actual, double expected, double tolerance, const char* description) {
    require(std::isfinite(actual) && std::abs(actual - expected) < tolerance, description);
}

void force(const MatrixXd& phi, MatrixXd& result, const vector<double>& coefficients,
           const langevin::Parameters& p) {
#ifdef COMMENSURATE
    Diss_comm(phi, result, coefficients.data(), p.Beta, p.L,
              langevin::LOCAL_FORCE_FACTOR * p.g / (Pi * Pi));
#else
    Diss(phi, result, coefficients.data(), p.Beta, p.L);
#endif
}

void test_kernel() {
    constexpr int EXTENT = 16;
    constexpr double HALF = 0.5;
    vector<double> coefficients(EXTENT - 1);
    Kernel(coefficients.data(), EXTENT, TEST_ALPHA, 1.0);
    for (int r = 1; r < EXTENT; ++r) {
        double expected = 0;
        for (int image = -EXTENT / 2; image < EXTENT / 2; ++image) {
            const double distance = r + image * EXTENT;
            expected += TEST_ALPHA / (Pi * Pi * (distance * distance - HALF * HALF));
        }
        near(coefficients[r - 1], expected, TOLERANCE, "Ohmic kernel normalization or image range");
    }
    for (double exponent : {0.5, 1.5}) {
        Kernel(coefficients.data(), EXTENT, TEST_ALPHA, exponent);
        double beta_value = std::beta(1.0 - exponent / 2, exponent + 1);
        // Independent beta recurrence for positive integer distances.
        vector<double> beta_by_distance(EXTENT * EXTENT + EXTENT, 0);
        for (std::size_t distance = 1; distance < beta_by_distance.size(); ++distance) {
            beta_by_distance[distance] = beta_value;
            beta_value *= (distance - exponent / 2) / (distance + exponent / 2 + 1);
        }
        for (int r = 1; r < EXTENT; ++r) {
            double expected = 0;
            for (int image = -EXTENT / 2; image < EXTENT / 2; ++image)
                expected += beta_by_distance[std::abs(r + image * EXTENT)];
            near(coefficients[r - 1], TEST_ALPHA * expected / (Pi * Pi), TOLERANCE,
                 "Generalized kernel differs from beta recurrence");
            require(coefficients[r - 1] > 0, "Kernel must be positive");
        }
    }
    Kernel(coefficients.data(), EXTENT, 0, 0.5);
    require(std::all_of(coefficients.begin(), coefficients.end(), [](double x) { return x == 0; }),
            "Zero bath coupling");
}

void test_laplacians_and_roughness() {
    MatrixXd phi(ROWS, COLS), result(ROWS, COLS);
    for (int row = 0; row < ROWS; ++row)
        for (int col = 0; col < COLS; ++col)
            phi(row, col) = cos(2 * Pi * col / COLS) + sin(2 * Pi * row / ROWS);
    Laplace_X(phi, result, ROWS, COLS);
    for (int row = 0; row < ROWS; ++row)
        for (int col = 0; col < COLS; ++col)
            near(result(row, col), -4 * pow(sin(Pi / COLS), 2) * cos(2 * Pi * col / COLS),
                 TOLERANCE, "Spatial periodic Laplacian Fourier eigenmode");
    Laplace_T(phi, result, ROWS, COLS);
    for (int row = 0; row < ROWS; ++row)
        for (int col = 0; col < COLS; ++col)
            near(result(row, col), -4 * pow(sin(Pi / ROWS), 2) * sin(2 * Pi * row / ROWS),
                 TOLERANCE, "Temporal periodic Laplacian Fourier eigenmode");
    vector<double> space(COLS - 1, 0), time(ROWS - 1, 0);
    RoughnessX(phi, space.data(), ROWS, COLS);
    RoughnessT(phi, time.data(), ROWS, COLS);
    for (int r = 1; r < COLS; ++r)
        near(space[r - 1], 1 - cos(2 * Pi * r / COLS), TOLERANCE, "Spatial roughness");
    for (int r = 1; r < ROWS; ++r)
        near(time[r - 1], 1 - cos(2 * Pi * r / ROWS), TOLERANCE, "Temporal roughness");
}

void test_action_gradient() {
    MatrixXd phi(ROWS, COLS), result(ROWS, COLS);
    for (int row = 0; row < ROWS; ++row)
        for (int col = 0; col < COLS; ++col) phi(row, col) = (row - col) * TEST_LOCAL_FORCE;
    // Use a symmetric periodic kernel to test differentiation of a scalar action.
    vector<double> coefficients{TEST_ALPHA, TEST_LOCAL_FORCE, TEST_ALPHA};
    const auto action = [&](const MatrixXd& field) {
        double value = 0;
        for (int col = 0; col < COLS; ++col)
            for (int row = 0; row < ROWS; ++row) {
#ifdef COMMENSURATE
                value -= TEST_LOCAL_FORCE / 4 * cos(4 * field(row, col));
#endif
                for (int other = row + 1; other < ROWS; ++other) {
                    const double weight = coefficients[other - row - 1] / 2;
                    value -= weight * cos(2 * (field(other, col) - field(row, col)));
#ifdef COMMENSURATE
                    value -= weight * cos(2 * (field(other, col) + field(row, col)));
#endif
                }
            }
        return value;
    };
#ifdef COMMENSURATE
    Diss_comm(phi, result, coefficients.data(), ROWS, COLS, TEST_LOCAL_FORCE);
#else
    Diss(phi, result, coefficients.data(), ROWS, COLS);
#endif
    for (int row = 0; row < ROWS; ++row)
        for (int col = 0; col < COLS; ++col) {
            MatrixXd plus = phi, minus = phi;
            plus(row, col) += DIFFERENCE_STEP;
            minus(row, col) -= DIFFERENCE_STEP;
            near(result(row, col), -(action(plus) - action(minus)) / (2 * DIFFERENCE_STEP),
                 FORCE_TOLERANCE, "Interaction drift must equal negative action gradient");
        }
}

void test_order_parameter_phase() {
    langevin::State state;
    state.p.L = COLS;
    state.p.Beta = ROWS;
    state.phi = MatrixXd::Constant(ROWS, COLS, Pi / 4);
    const MatrixXd original = state.phi;
    state.rough_x.assign(COLS - 1, 0);
    state.rough_tau.assign(ROWS - 1, 0);
#ifdef COMMENSURATE
    state.commensurate = true;
    near(langevin::measure(state), 0, TOLERANCE, "Commensurate measurements must preserve absolute phase");
#else
    near(langevin::measure(state), 1, TOLERANCE, "Incommensurate observable uses mean-subtracted field");
#endif
    require(state.phi == original, "Measurements must not shift the simulation field");
#ifdef COMMENSURATE
    langevin::Parameters p;
    p.L = COLS;
    p.Beta = ROWS;
    p.g = TEST_LOCAL_FORCE;
    vector<double> coefficients(ROWS - 1, 0);
    MatrixXd result(ROWS, COLS);
    state.phi.setConstant(Pi / 8);
    force(state.phi, result, coefficients, p);
    near(result(0, 0), -2 * p.g / (Pi * Pi), TOLERANCE, "Physical g to local force conversion");
#endif
}

void test_harmonic_equilibrium() {
    langevin::State state;
    state.p.L = state.p.Beta = 2;
    state.p.alpha = 0;
    state.p.dt = HARMONIC_DT;
    state.generator.seed(langevin::DEFAULT_SEED);
    state.phi = MatrixXd::Zero(state.p.Beta, state.p.L);
    vector<double> coefficients(state.p.Beta - 1, 0);
    double sum = 0, squares = 0;
    int samples = 0;
    for (int step = 0; step < HARMONIC_STEPS; ++step) {
        langevin::advance(state, coefficients, force);
        if (step >= HARMONIC_BURN_IN && step % HARMONIC_STRIDE == 0) {
            const double mode = (state.phi(0, 0) - state.phi(0, 1) + state.phi(1, 0) - state.phi(1, 1)) / 4;
            sum += mode;
            squares += mode * mode;
            ++samples;
        }
    }
    const double variance = squares / samples - pow(sum / samples, 2);
    const double equilibrium = Pi / 16; // Noise covariance 2 and spatial decay rate 4/pi.
    near(variance, equilibrium, VARIANCE_TOLERANCE, "Harmonic equilibrium variance / noise normalization");
}
} // namespace

int main() {
    try {
        test_kernel();
        test_laplacians_and_roughness();
        test_action_gradient();
        test_order_parameter_phase();
        test_harmonic_equilibrium();
        std::cout << "Passed kernel, Fourier Laplacians, rectangular roughness, action-gradient, field-phase and harmonic equilibrium tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
