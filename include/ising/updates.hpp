// Two Markov chains with the Boltzmann distribution as their stationary law.
#pragma once

#include <cmath>
#include <vector>

#include "lattice.hpp"

namespace ising {

/// Single-spin-flip Metropolis. A sweep visits every site once in typewriter
/// order and flips s_i with probability min(1, exp(-beta * dE)).
struct Metropolis {
    double beta;
    double accept[2];  // exp(-4 beta), exp(-8 beta): the only positive dE values

    explicit Metropolis(double b) : beta(b), accept{std::exp(-4 * b), std::exp(-8 * b)} {}

    /// One sweep. E and M are kept up to date incrementally.
    void sweep(Lattice& lat, Rng& rng, int& E, int& M) const {
        for (int i = 0; i < lat.N; ++i) {
            const int s = lat.spin[i];
            const int dE = 2 * s * lat.neighbour_sum(i);  // -8, -4, 0, 4 or 8
            if (dE <= 0 || uniform(rng) < accept[(dE >> 2) - 1]) {
                lat.spin[i] = static_cast<std::int8_t>(-s);
                E += dE;
                M -= 2 * s;
            }
        }
    }
};

/// Wolff (1989) single-cluster update: grow a cluster of like spins from a
/// random seed, adding each aligned neighbour with probability
/// p = 1 - exp(-2 beta), then flip the whole cluster. Every cluster move is
/// accepted, and near T_c the clusters are as large as the correlated regions.
struct Wolff {
    double p_add;
    std::vector<std::int32_t> stack;

    explicit Wolff(double beta) : p_add(1 - std::exp(-2 * beta)) {}

    /// Grow and flip one cluster. Returns its size.
    int step(Lattice& lat, Rng& rng) {
        const int seed = below(rng, lat.N);
        const std::int8_t s0 = lat.spin[seed];
        lat.spin[seed] = static_cast<std::int8_t>(-s0);  // flipped spins count as visited
        stack.clear();
        stack.push_back(seed);
        int size = 1;
        while (!stack.empty()) {
            const int i = stack.back();
            stack.pop_back();
            for (const int j : {lat.right[i], lat.left[i], lat.up[i], lat.down[i]}) {
                if (lat.spin[j] == s0 && uniform(rng) < p_add) {
                    lat.spin[j] = static_cast<std::int8_t>(-s0);
                    stack.push_back(j);
                    ++size;
                }
            }
        }
        return size;
    }

    /// Flip `k` clusters. Returns the total number of spins flipped.
    long long steps(Lattice& lat, Rng& rng, int k) {
        long long flipped = 0;
        for (int i = 0; i < k; ++i) flipped += step(lat, rng);
        return flipped;
    }

    /// Flip clusters until at least N spins have flipped: one sweep-equivalent.
    /// The stopping time depends on the clusters, so the state reached this way
    /// is NOT a fair equilibrium sample (large clusters, typical of ordered
    /// states, end a sweep early). Use it only to equilibrate; measure after a
    /// fixed number of steps instead. Returns the number of clusters flipped.
    int sweep(Lattice& lat, Rng& rng) {
        int flipped = 0, clusters = 0;
        while (flipped < lat.N) {
            flipped += step(lat, rng);
            ++clusters;
        }
        return clusters;
    }
};

}  // namespace ising
