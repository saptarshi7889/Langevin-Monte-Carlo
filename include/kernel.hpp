#pragma once

#include <cmath>
#include <stdexcept>

namespace langevin {
inline constexpr double PI = 3.14159265358979323846;
inline constexpr double MIN_BATH_EXPONENT = 0.0;
inline constexpr double MAX_BATH_EXPONENT = 2.0;
inline constexpr double HALF = 0.5;

// Reflect negative integer images analytically before evaluating the beta function.
// For integer n != 0 and 0 < s < 2, B(n-s/2,s+1) is even in n.
inline void kernel(double* values, int extent, double alpha, double exponent) {
    if (extent < 2 || extent % 2 != 0 || !std::isfinite(alpha) || alpha < 0 ||
        !std::isfinite(exponent) || exponent <= MIN_BATH_EXPONENT ||
        exponent >= MAX_BATH_EXPONENT) {
        throw std::invalid_argument("Kernel requires an even Beta >= 2, alpha >= 0 and 0 < s < 2");
    }
    for (int r = 1; r < extent; ++r) {
        double sum = 0.0;
        // Preserve the original half-open image range.
        for (int image = -extent / 2; image < extent / 2; ++image) {
            const double distance = std::abs(static_cast<double>(r) +
                                             static_cast<double>(image) * extent);
            const double argument = distance - HALF * exponent;
            sum += std::exp(std::lgamma(argument) + std::lgamma(exponent + 1.0) -
                            std::lgamma(argument + exponent + 1.0));
        }
        values[r - 1] = alpha * sum / (PI * PI);
    }
}
} // namespace langevin
