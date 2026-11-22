// Run one simulation and print the observables with error bars.
//
//   ising-run --L 64 --T 2.2692 --algo wolff --sweeps 20000
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "ising/exact.hpp"
#include "ising/simulate.hpp"

using namespace ising;

int main(int argc, char** argv) {
    int L = 32, sweeps = 20000, equilibration = 2000;
    double T = critical_temperature();
    unsigned long long seed = 1;
    Algorithm algo = Algorithm::Wolff;
    for (int i = 1; i + 1 < argc; i += 2) {
        const std::string k = argv[i], v = argv[i + 1];
        if (k == "--L") L = std::atoi(v.c_str());
        else if (k == "--T") T = std::atof(v.c_str());
        else if (k == "--sweeps") sweeps = std::atoi(v.c_str());
        else if (k == "--equilibration") equilibration = std::atoi(v.c_str());
        else if (k == "--seed") seed = std::strtoull(v.c_str(), nullptr, 10);
        else if (k == "--algo") algo = (v == "metropolis") ? Algorithm::Metropolis : Algorithm::Wolff;
        else {
            std::fprintf(stderr,
                         "usage: ising-run [--L 32] [--T 2.2692] [--algo wolff|metropolis] [--sweeps 20000]\n"
                         "                 [--equilibration 2000] [--seed 1]\n");
            return 2;
        }
    }
    if (L < 2 || T <= 0 || sweeps < 100) {
        std::fprintf(stderr, "need L >= 2, T > 0 and at least 100 sweeps\n");
        return 2;
    }
    Lattice lat(L);
    Rng rng(seed);
    double spp = 1.0;  // sweeps per sample (Wolff: k clusters per sample)
    const auto xs = run(lat, T, algo, equilibration, sweeps, rng, &spp);
    const auto o = observables(xs, T, lat.N);
    std::vector<double> m(xs.size()), e(xs.size());
    for (std::size_t i = 0; i < xs.size(); ++i) {
        m[i] = std::fabs(xs[i].m);
        e[i] = xs[i].e;
    }
    const auto tm = integrated_autocorrelation_time(m), te = integrated_autocorrelation_time(e);
    std::printf("%s, L = %d, T = %.5f, %d samples (%.3f sweeps each)\n", name(algo), L, T, sweeps, spp);
    auto row = [](const char* label, Estimate x) { std::printf("  %-22s %.6f +- %.6f\n", label, x.value, x.error); };
    row("energy per spin", o.e);
    row("|m|", o.abs_m);
    row("specific heat", o.specific_heat);
    row("susceptibility", o.susceptibility);
    row("Binder cumulant", o.binder);
    std::printf("  tau_int(|m|)          %.2f +- %.2f sweeps\n", tm.tau * spp, tm.error * spp);
    std::printf("  tau_int(e)            %.2f +- %.2f sweeps\n", te.tau * spp, te.error * spp);
    if (T != critical_temperature())
        std::printf("  Onsager (infinite L): e = %.6f, m = %.6f\n", onsager_energy(T), onsager_magnetization(T));
    return 0;
}
