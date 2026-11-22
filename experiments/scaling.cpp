// Finite-size scaling at T_c: <|m|> ~ L^(-beta/nu) and N <m^2> / T ~ L^(gamma/nu),
// with exact 2D Ising exponents beta/nu = 1/8 and gamma/nu = 7/4.
#include "common.hpp"
#include "ising/exact.hpp"
#include "ising/simulate.hpp"

using namespace ising;

int main() {
    const double tc = critical_temperature();
    std::FILE* f = open_csv("scaling.csv", "L,abs_m,abs_m_err,chi_m2,chi_m2_err,binder,binder_err,e,e_err,sweeps");
    for (const int L : {8, 16, 32, 64, 128, 256}) {
        const int sweeps = L <= 64 ? 50000 : (L == 128 ? 20000 : 10000);
        Lattice lat(L);
        Rng rng(5150u + static_cast<unsigned>(L));
        const auto xs = run(lat, tc, Algorithm::Wolff, 2000, sweeps, rng);
        const auto o = observables(xs, tc, lat.N);
        std::fprintf(f, "%d,%.7f,%.7f,%.5f,%.5f,%.6f,%.6f,%.7f,%.7f,%d\n", L, o.abs_m.value, o.abs_m.error,
                     o.chi_m2.value, o.chi_m2.error, o.binder.value, o.binder.error, o.e.value, o.e.error, sweeps);
        std::fflush(f);
        std::printf("L = %d done\n", L);
    }
    std::fclose(f);
}
