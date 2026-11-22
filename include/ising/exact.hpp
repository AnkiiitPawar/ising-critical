// Exact results to test the simulations against.
//
// 1. Small lattices: enumerate all 2^N states (Gray code, one spin flip per
//    step) into the density of states g(E, M). Every equilibrium average at
//    every temperature then follows exactly.
// 2. The infinite lattice: Onsager's energy (1944), the Onsager-Yang
//    spontaneous magnetization (Yang 1952) and T_c = 2 / ln(1 + sqrt 2).
#pragma once

#include <bit>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "lattice.hpp"

namespace ising {

inline constexpr double kPi = 3.14159265358979323846;

/// T_c = 2 / ln(1 + sqrt 2) = 2.269185...
inline double critical_temperature() { return 2.0 / std::log(1.0 + std::sqrt(2.0)); }

/// Complete elliptic integral of the first kind K(k) (modulus k), from the
/// arithmetic-geometric mean: K(k) = pi / (2 AGM(1, sqrt(1 - k^2))).
inline double elliptic_k(double k) {
    double a = 1.0, b = std::sqrt(1.0 - k * k);
    for (int i = 0; i < 60 && std::fabs(a - b) > 1e-16 * a; ++i) {
        const double an = 0.5 * (a + b);
        b = std::sqrt(a * b);
        a = an;
    }
    return kPi / (2.0 * a);
}

/// Onsager's internal energy per spin of the infinite lattice:
/// u = -coth(2b) [1 + (2/pi)(2 tanh^2(2b) - 1) K(kappa)],  kappa = 2 sinh(2b) / cosh^2(2b).
inline double onsager_energy(double T) {
    const double b = 1.0 / T;
    const double t2 = std::tanh(2 * b) * std::tanh(2 * b);
    const double kappa = 2 * std::sinh(2 * b) / (std::cosh(2 * b) * std::cosh(2 * b));
    const double coth = 1.0 / std::tanh(2 * b);
    if (std::fabs(2 * t2 - 1) < 1e-15) return -coth;  // at T_c the K term vanishes
    return -coth * (1 + (2 / kPi) * (2 * t2 - 1) * elliptic_k(kappa));
}

/// Spontaneous magnetization per spin: (1 - sinh^-4(2b))^(1/8) below T_c, 0 above.
inline double onsager_magnetization(double T) {
    const double s = std::sinh(2.0 / T);
    const double x = 1.0 - 1.0 / (s * s * s * s);
    return x > 0 ? std::pow(x, 0.125) : 0.0;
}

/// Equilibrium averages per spin: e = E/N and m = M/N.
struct Averages {
    double e, e2, abs_m, m2, m4;
    double specific_heat(double T, int N) const { return N * (e2 - e * e) / (T * T); }
    double binder() const { return 1.0 - m4 / (3.0 * m2 * m2); }
};

/// g(E, M) for an L x L periodic lattice, by visiting all 2^N states.
struct DensityOfStates {
    int L, N;
    std::vector<std::uint64_t> count;  // index ((E + 2N) / 4) * (N + 1) + (M + N) / 2

    int energy_levels() const { return N + 1; }
    std::uint64_t at(int E, int M) const { return count[((E + 2 * N) / 4) * (N + 1) + (M + N) / 2]; }

    Averages averages(double T) const {
        long double z = 0, e = 0, e2 = 0, am = 0, m2 = 0, m4 = 0;
        for (int ei = 0; ei < energy_levels(); ++ei) {
            const int E = 4 * ei - 2 * N;
            const long double w = std::exp(-static_cast<long double>(E + 2 * N) / T);  // shifted by E_min
            for (int mi = 0; mi <= N; ++mi) {
                const std::uint64_t g = count[ei * (N + 1) + mi];
                if (g == 0) continue;
                const long double p = w * g;
                const long double en = static_cast<long double>(E) / N;
                const long double mn = static_cast<long double>(2 * mi - N) / N;
                z += p;
                e += p * en;
                e2 += p * en * en;
                am += p * std::fabs(mn);
                m2 += p * mn * mn;
                m4 += p * mn * mn * mn * mn;
            }
        }
        return {static_cast<double>(e / z), static_cast<double>(e2 / z), static_cast<double>(am / z),
                static_cast<double>(m2 / z), static_cast<double>(m4 / z)};
    }
};

/// Enumerate all states of an L x L lattice in Gray-code order: state t
/// differs from state t-1 in the spin at the lowest set bit of t, so E and M
/// change by one local update per step. Practical up to L = 5 (2^25 states).
inline DensityOfStates enumerate(int L) {
    if (L < 2 || L * L > 30) throw std::invalid_argument("enumerate: need 2 <= L and L*L <= 30");
    Lattice lat(L);
    lat.fill(-1);
    const int N = lat.N;
    DensityOfStates dos{L, N, std::vector<std::uint64_t>(static_cast<std::size_t>(N + 1) * (N + 1), 0)};
    int E = lat.energy(), M = lat.magnetization();
    const std::uint64_t states = std::uint64_t{1} << N;
    for (std::uint64_t t = 0; t < states; ++t) {
        if (t > 0) {
            const int i = std::countr_zero(t);
            const int s = lat.spin[i];
            E += 2 * s * lat.neighbour_sum(i);
            M -= 2 * s;
            lat.spin[i] = static_cast<std::int8_t>(-s);
        }
        ++dos.count[((E + 2 * N) / 4) * (N + 1) + (M + N) / 2];
    }
    return dos;
}

}  // namespace ising
