// The two-dimensional Ising model on an L x L square lattice with periodic
// boundaries: spins s_i = +-1, energy E = -J sum_<ij> s_i s_j with J = 1 and
// no external field. Temperatures are in units of J / k_B.
#pragma once

#include <cstdint>
#include <random>
#include <vector>

namespace ising {

using Rng = std::mt19937_64;

/// Uniform double in [0, 1) from the top 53 bits of one draw.
inline double uniform(Rng& rng) { return static_cast<double>(rng() >> 11) * 0x1.0p-53; }

/// Uniform integer in [0, n), by multiply-and-shift (bias below 2^-32 for n < 2^32).
inline int below(Rng& rng, int n) {
    return static_cast<int>(((rng() >> 32) * static_cast<std::uint64_t>(n)) >> 32);
}

struct Lattice {
    int L;
    int N;
    std::vector<std::int8_t> spin;
    std::vector<std::int32_t> right, left, up, down;

    explicit Lattice(int side) : L(side), N(side * side), spin(N, 1), right(N), left(N), up(N), down(N) {
        for (int y = 0; y < L; ++y) {
            for (int x = 0; x < L; ++x) {
                const int i = y * L + x;
                right[i] = y * L + (x + 1) % L;
                left[i] = y * L + (x + L - 1) % L;
                down[i] = ((y + 1) % L) * L + x;
                up[i] = ((y + L - 1) % L) * L + x;
            }
        }
    }

    int neighbour_sum(int i) const { return spin[right[i]] + spin[left[i]] + spin[up[i]] + spin[down[i]]; }

    /// E = -sum over bonds of s_i s_j, each bond counted once (right and down).
    int energy() const {
        int e = 0;
        for (int i = 0; i < N; ++i) e -= spin[i] * (spin[right[i]] + spin[down[i]]);
        return e;
    }

    int magnetization() const {
        int m = 0;
        for (auto s : spin) m += s;
        return m;
    }

    void fill(std::int8_t s) { spin.assign(N, s); }

    void randomize(Rng& rng) {
        for (auto& s : spin) s = (rng() >> 63) ? 1 : -1;
    }
};

}  // namespace ising
