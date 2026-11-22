// Energy, magnetization, specific heat and susceptibility across T_c for
// L = 16 and 64, next to Onsager's exact infinite-lattice results.
#include <vector>

#include "common.hpp"
#include "ising/exact.hpp"
#include "ising/simulate.hpp"

using namespace ising;

int main() {
    std::vector<double> temps;
    for (int k = 0; k <= 20; ++k) temps.push_back(1.5 + 0.1 * k);
    for (int k = 0; k <= 7; ++k) temps.push_back(2.20 + 0.02 * k);
    std::FILE* f = open_csv("thermodynamics.csv",
                            "L,T,e,e_err,abs_m,abs_m_err,c,c_err,chi,chi_err,binder,binder_err,onsager_e,onsager_m");
    for (const int L : {16, 64}) {
        const int sweeps = L == 16 ? 50000 : 10000;
        for (const double T : temps) {
            Lattice lat(L);
            Rng rng(1000u * static_cast<unsigned>(L) + static_cast<unsigned>(T * 1000));
            const auto xs = run(lat, T, Algorithm::Wolff, 2000, sweeps, rng);
            const auto o = observables(xs, T, lat.N);
            std::fprintf(f, "%d,%.4f,%.7f,%.7f,%.7f,%.7f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.7f,%.7f\n", L, T, o.e.value,
                         o.e.error, o.abs_m.value, o.abs_m.error, o.specific_heat.value, o.specific_heat.error,
                         o.susceptibility.value, o.susceptibility.error, o.binder.value, o.binder.error,
                         onsager_energy(T), onsager_magnetization(T));
            std::fflush(f);
        }
        std::printf("L = %d done\n", L);
    }
    std::fclose(f);
}
