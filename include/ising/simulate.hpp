// Running a chain and turning its samples into thermodynamic observables.
#pragma once

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "lattice.hpp"
#include "stats.hpp"
#include "updates.hpp"

namespace ising {

enum class Algorithm { Metropolis, Wolff };

inline const char* name(Algorithm a) { return a == Algorithm::Metropolis ? "metropolis" : "wolff"; }

/// Equilibrate for `equilibration` sweeps, then record `sweeps` samples.
///
/// Metropolis: one sample after every sweep. Wolff: one sample after every k
/// cluster flips, with k = round(N / <cluster size>) fixed from the
/// equilibration phase, so one sample costs about one sweep-equivalent.
/// Sampling after a fixed number of steps keeps the samples unbiased; stopping
/// once N spins have flipped would not (see Wolff::sweep). `sweeps_per_sample`,
/// if given, receives the measured average cost of one sample in sweeps.
inline std::vector<Sample> run(Lattice& lat, double T, Algorithm algo, int equilibration, int sweeps, Rng& rng,
                               double* sweeps_per_sample = nullptr) {
    std::vector<Sample> out;
    out.reserve(sweeps);
    const double inv_n = 1.0 / lat.N;
    if (algo == Algorithm::Metropolis) {
        const Metropolis mc(1.0 / T);
        int E = lat.energy(), M = lat.magnetization();
        for (int t = 0; t < equilibration; ++t) mc.sweep(lat, rng, E, M);
        for (int t = 0; t < sweeps; ++t) {
            mc.sweep(lat, rng, E, M);
            out.push_back({E * inv_n, M * inv_n});
        }
        if (sweeps_per_sample) *sweeps_per_sample = 1.0;
    } else {
        Wolff mc(1.0 / T);
        for (int t = 0; t < equilibration; ++t) mc.sweep(lat, rng);
        // Calibrate k: clusters per sweep-equivalent, averaged over further equilibrated sweeps.
        const int calibration = std::max(100, equilibration / 2);
        long long clusters = 0;
        for (int t = 0; t < calibration; ++t) clusters += mc.sweep(lat, rng);
        const int k = std::max(1, static_cast<int>(std::lround(static_cast<double>(clusters) / calibration)));
        long long flipped = 0;
        for (int t = 0; t < sweeps; ++t) {
            flipped += mc.steps(lat, rng, k);
            out.push_back({lat.energy() * inv_n, lat.magnetization() * inv_n});
        }
        if (sweeps_per_sample) *sweeps_per_sample = static_cast<double>(flipped) / sweeps / lat.N;
    }
    return out;
}

/// The standard observables with jackknife errors.
struct Observables {
    Estimate e, abs_m, specific_heat, susceptibility, chi_m2, binder;
};

inline Observables observables(const std::vector<Sample>& xs, double T, int N, int bins = 50) {
    return {
        jackknife(xs, bins, [](const Moments& m) { return m[0]; }),
        jackknife(xs, bins, [](const Moments& m) { return m[2]; }),
        jackknife(xs, bins, [&](const Moments& m) { return N * (m[1] - m[0] * m[0]) / (T * T); }),
        // Finite-lattice susceptibility built from |m|, which stays finite at every T.
        jackknife(xs, bins, [&](const Moments& m) { return N * (m[3] - m[2] * m[2]) / T; }),
        // N <m^2> / T: grows like L^(gamma/nu) at T_c.
        jackknife(xs, bins, [&](const Moments& m) { return N * m[3] / T; }),
        jackknife(xs, bins, [](const Moments& m) { return 1.0 - m[4] / (3.0 * m[3] * m[3]); }),
    };
}

/// Single-histogram reweighting (Ferrenberg and Swendsen 1988): samples drawn
/// at temperature T0 also estimate averages at a nearby T, with weights
/// exp(-(1/T - 1/T0) E). Reliable while T stays within the range of energies
/// the run actually visited.
inline Moments reweight(const std::vector<Sample>& xs, int N, double T0, double T) {
    const double db = 1.0 / T - 1.0 / T0;
    double max_exp = -1e300;
    for (const auto& s : xs) max_exp = std::max(max_exp, -db * s.e * N);
    Moments mo{};
    double z = 0;
    for (const auto& s : xs) {
        const double w = std::exp(-db * s.e * N - max_exp);
        const double m2 = s.m * s.m;
        z += w;
        mo[0] += w * s.e;
        mo[1] += w * s.e * s.e;
        mo[2] += w * std::fabs(s.m);
        mo[3] += w * m2;
        mo[4] += w * m2 * m2;
    }
    for (auto& x : mo) x /= z;
    return mo;
}

}  // namespace ising
