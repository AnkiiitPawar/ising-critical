// Critical slowing down: the integrated autocorrelation time at T_c grows
// like L^z. Local Metropolis updates have z near 2.17; Wolff clusters make z
// small. Times are converted to sweeps: for Wolff, one sample costs k cluster
// flips, and the measured spins flipped per sample / N gives its cost in sweeps.
#include <cmath>
#include <vector>

#include "common.hpp"
#include "ising/exact.hpp"
#include "ising/simulate.hpp"

using namespace ising;

static void measure(std::FILE* f, Algorithm algo, int L, int sweeps) {
    const double tc = critical_temperature();
    Lattice lat(L);
    Rng rng(31337u + 7u * static_cast<unsigned>(L) + (algo == Algorithm::Wolff ? 1u : 0u));
    double spp = 1.0;  // sweeps per sample
    const auto xs = run(lat, tc, algo, 20000, sweeps, rng, &spp);
    std::vector<double> m(xs.size()), e(xs.size());
    for (std::size_t i = 0; i < xs.size(); ++i) {
        m[i] = std::fabs(xs[i].m);
        e[i] = xs[i].e;
    }
    for (const auto& [obs, series] : {std::pair{"abs_m", &m}, std::pair{"e", &e}}) {
        const auto t = integrated_autocorrelation_time(*series);
        std::fprintf(f, "%s,%d,%s,%.4f,%.4f,%zu,%d,%.4f\n", name(algo), L, obs, t.tau * spp, t.error * spp, t.window,
                     sweeps, spp);
        std::printf("%-10s L = %3d  tau(%s) = %.2f +- %.2f sweeps (%.3f sweeps per sample)\n", name(algo), L, obs,
                    t.tau * spp, t.error * spp, spp);
    }
    std::fflush(f);
}

int main() {
    std::FILE* f = open_csv("dynamics.csv", "algorithm,L,observable,tau,tau_err,window,samples,sweeps_per_sample");
    for (const int L : {8, 12, 16, 24, 32}) measure(f, Algorithm::Metropolis, L, L <= 16 ? 1000000 : 2000000);
    for (const int L : {8, 16, 32, 64, 128}) measure(f, Algorithm::Wolff, L, L <= 64 ? 100000 : 50000);
    std::fclose(f);
}
