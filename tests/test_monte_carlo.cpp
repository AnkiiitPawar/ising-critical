// Both Markov chains must reproduce exact equilibrium averages on a lattice
// small enough to enumerate. Seeds are fixed, so the outcome is deterministic.
#include <cstdio>

#include "check.hpp"
#include "ising/exact.hpp"
#include "ising/simulate.hpp"

using namespace ising;

// |estimate - exact| must be within 5 standard errors, and the error must be
// small enough for the comparison to mean something.
static void agree(const char* what, Estimate got, double exact, double max_error) {
    CHECK(got.error > 0);
    CHECK(got.error < max_error);
    CHECK_NEAR(got.value, exact, 5 * got.error);
    std::printf("  %-26s %+.5f +- %.5f   exact %+.5f\n", what, got.value, got.error, exact);
}

static void matches_exact_enumeration() {
    const auto dos = enumerate(4);
    for (const Algorithm algo : {Algorithm::Metropolis, Algorithm::Wolff}) {
        for (const double T : {1.8, critical_temperature(), 3.0}) {
            Lattice lat(4);
            Rng rng(12345);
            const auto xs = run(lat, T, algo, 2000, 200000, rng);
            const auto obs = observables(xs, T, lat.N);
            const auto ex = dos.averages(T);
            std::printf("%s, L = 4, T = %.4f\n", name(algo), T);
            agree("energy per spin", obs.e, ex.e, 0.005);
            agree("|m|", obs.abs_m, ex.abs_m, 0.005);
            agree("specific heat", obs.specific_heat, ex.specific_heat(T, lat.N), 0.02);
            agree("N <m^2> / T", obs.chi_m2, lat.N * ex.m2 / T, 0.05);
            agree("Binder cumulant", obs.binder, ex.binder(), 0.01);
        }
    }
}

static void incremental_bookkeeping() {
    Lattice lat(8);
    Rng rng(3);
    lat.randomize(rng);
    int E = lat.energy(), M = lat.magnetization();
    const Metropolis mc(1 / 2.3);
    for (int t = 0; t < 500; ++t) mc.sweep(lat, rng, E, M);
    CHECK(E == lat.energy());
    CHECK(M == lat.magnetization());
}

static void wolff_limits() {
    Rng rng(1);
    Lattice lat(8);
    Wolff cold(50.0);  // p_add = 1 - e^-100 rounds to 1: the cluster is the whole ordered lattice
    CHECK(cold.step(lat, rng) == lat.N);
    CHECK(lat.magnetization() == -lat.N);
    Wolff hot(0.0);  // p_add = 0: only the seed
    CHECK(hot.step(lat, rng) == 1);
}

static void reweighting() {
    const auto dos = enumerate(4);
    const double t0 = critical_temperature();
    Lattice lat(4);
    Rng rng(99);
    const auto xs = run(lat, t0, Algorithm::Wolff, 2000, 400000, rng);
    // At the simulation temperature, reweighting changes nothing.
    const auto same = reweight(xs, lat.N, t0, t0);
    const auto plain = moments_of(xs.data(), xs.data() + xs.size());
    for (int k = 0; k < 5; ++k) CHECK_NEAR(same[k], plain[k], 1e-12);
    // Nearby temperatures, from the same samples.
    for (const double T : {2.1, 2.45}) {
        const auto rw = reweight(xs, lat.N, t0, T);
        const auto ex = dos.averages(T);
        CHECK_NEAR(rw[0], ex.e, 0.01);
        CHECK_NEAR(1 - rw[4] / (3 * rw[3] * rw[3]), ex.binder(), 0.01);
    }
}

int main() {
    matches_exact_enumeration();
    incremental_bookkeeping();
    wolff_limits();
    reweighting();
    return finish("test_monte_carlo");
}
