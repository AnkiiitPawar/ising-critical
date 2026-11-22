// Equilibrium configurations of a 256 x 256 lattice below, at and above T_c,
// saved as binary PBM images (1 bit per spin).
#include <cstdio>
#include <string>

#include "common.hpp"
#include "ising/exact.hpp"
#include "ising/simulate.hpp"

using namespace ising;

int main() {
    const int L = 256;
    for (const double T : {2.0, critical_temperature(), 2.6}) {
        Lattice lat(L);
        Rng rng(2024u + static_cast<unsigned>(T * 100));
        lat.randomize(rng);
        Wolff mc(1.0 / T);
        for (int t = 0; t < 3000; ++t) mc.sweep(lat, rng);
        char name[64];
        std::snprintf(name, sizeof name, "docs/data/snapshot_T%.3f.pbm", T);
        std::FILE* f = std::fopen(name, "wb");
        std::fprintf(f, "P4\n%d %d\n", L, L);
        for (int y = 0; y < L; ++y) {
            for (int x = 0; x < L; x += 8) {
                unsigned char byte = 0;
                for (int b = 0; b < 8; ++b)
                    if (lat.spin[y * L + x + b] > 0) byte = static_cast<unsigned char>(byte | (0x80 >> b));
                std::fputc(byte, f);
            }
        }
        std::fclose(f);
        std::printf("wrote %s (m = %+.3f)\n", name, static_cast<double>(lat.magnetization()) / lat.N);
    }
}
