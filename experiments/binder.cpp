// Binder cumulant U4 = 1 - <m^4> / (3 <m^2>^2) near T_c for L = 8 ... 64.
//
// U4 tends to 2/3 in the ordered phase and to 0 in the disordered phase. At
// T_c it takes a universal, size-independent value, so the curves for
// different L cross there. Smooth curves come from one long run per L at T_c
// with histogram reweighting; independent runs at five temperatures check
// the reweighting.
#include <vector>

#include "common.hpp"
#include "ising/exact.hpp"
#include "ising/simulate.hpp"

using namespace ising;

int main() {
    const double tc = critical_temperature();
    std::FILE* curves = open_csv("binder_curves.csv", "L,T,binder,chi_m2");
    std::FILE* points = open_csv("binder_points.csv", "L,T,binder,binder_err");
    for (const int L : {8, 16, 32, 64}) {
        const int sweeps = L <= 32 ? 200000 : 100000;
        Lattice lat(L);
        Rng rng(77u + static_cast<unsigned>(L));
        const auto xs = run(lat, tc, Algorithm::Wolff, 5000, sweeps, rng);
        for (int k = 0; k <= 100; ++k) {
            const double T = 2.22 + 0.001 * k;
            const auto m = reweight(xs, lat.N, tc, T);
            std::fprintf(curves, "%d,%.4f,%.6f,%.6f\n", L, T, 1 - m[4] / (3 * m[3] * m[3]), lat.N * m[3] / T);
        }
        for (const double T : {2.24, 2.255, 2.27, 2.285, 2.30}) {
            Lattice l2(L);
            Rng r2(911u * static_cast<unsigned>(L) + static_cast<unsigned>(T * 1e3));
            const auto ys = run(l2, T, Algorithm::Wolff, 2000, 40000, r2);
            const auto o = observables(ys, T, l2.N);
            std::fprintf(points, "%d,%.4f,%.6f,%.6f\n", L, T, o.binder.value, o.binder.error);
        }
        std::fflush(curves);
        std::fflush(points);
        std::printf("L = %d done\n", L);
    }
    std::fclose(curves);
    std::fclose(points);
}
