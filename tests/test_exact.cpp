// The exact results the simulations are tested against must themselves be right.
#include <cstdint>

#include "check.hpp"
#include "ising/exact.hpp"

using namespace ising;

// Enumerate by setting every spin from the bits of the state index and
// recomputing E and M from scratch: slow, but independent of the Gray code.
static DensityOfStates naive(int L) {
    Lattice lat(L);
    const int N = lat.N;
    DensityOfStates dos{L, N, std::vector<std::uint64_t>(static_cast<std::size_t>(N + 1) * (N + 1), 0)};
    for (std::uint64_t t = 0; t < (std::uint64_t{1} << N); ++t) {
        for (int i = 0; i < N; ++i) lat.spin[i] = ((t >> i) & 1) ? 1 : -1;
        const int E = lat.energy(), M = lat.magnetization();
        ++dos.count[((E + 2 * N) / 4) * (N + 1) + (M + N) / 2];
    }
    return dos;
}

static void lattice_basics() {
    Lattice lat(3);
    CHECK(lat.energy() == -2 * lat.N);  // all spins up: every one of the 2N bonds satisfied
    CHECK(lat.magnetization() == lat.N);
    CHECK(lat.right[2] == 0 && lat.left[0] == 2 && lat.up[0] == 6 && lat.down[6] == 0);
    // Flipping one spin changes E by 2 s h, the rule every update relies on.
    Rng rng(7);
    Lattice big(6);
    big.randomize(rng);
    for (int k = 0; k < 200; ++k) {
        const int i = below(rng, big.N);
        const int before = big.energy();
        const int dE = 2 * big.spin[i] * big.neighbour_sum(i);
        big.spin[i] = static_cast<std::int8_t>(-big.spin[i]);
        CHECK(big.energy() - before == dE);
    }
}

static void density_of_states() {
    for (int L : {3, 4}) {
        const auto gray = enumerate(L), slow = naive(L);
        CHECK(gray.count == slow.count);
    }
    for (int L : {3, 4, 5}) {
        const auto dos = enumerate(L);
        const int N = dos.N;
        std::uint64_t total = 0;
        for (auto c : dos.count) total += c;
        CHECK(total == (std::uint64_t{1} << N));
        CHECK(dos.at(-2 * N, N) == 1 && dos.at(-2 * N, -N) == 1);  // the two ground states
        std::uint64_t first = 0;  // first excitation: one flipped spin, 4 broken bonds
        for (int M = -N; M <= N; M += 2) first += dos.at(-2 * N + 8, M);
        CHECK(first == static_cast<std::uint64_t>(2 * N));
        for (int E = -2 * N; E <= 2 * N; E += 4)
            for (int M = -N; M <= N; M += 2) CHECK(dos.at(E, M) == dos.at(E, -M));
    }
    // Infinite temperature: independent spins, so <e> = 0 and <m^2> = 1/N.
    const auto a = enumerate(4).averages(1e9);
    CHECK_NEAR(a.e, 0.0, 1e-8);
    CHECK_NEAR(a.m2, 1.0 / 16, 1e-8);
    // Very low temperature: ordered, Binder cumulant -> 2/3.
    const auto b = enumerate(4).averages(0.2);
    CHECK_NEAR(b.abs_m, 1.0, 1e-6);
    CHECK_NEAR(b.binder(), 2.0 / 3.0, 1e-6);
}

static void onsager() {
    CHECK_NEAR(critical_temperature(), 2.269185314213022, 1e-14);
    CHECK_NEAR(elliptic_k(0.0), kPi / 2, 1e-15);
    CHECK_NEAR(elliptic_k(1 / std::sqrt(2.0)), 1.8540746773013719, 1e-14);  // Gamma(1/4)^2 / (4 sqrt pi)
    const double tc = critical_temperature();
    CHECK_NEAR(onsager_energy(tc), -std::sqrt(2.0), 1e-12);  // u(T_c) = -sqrt 2
    // Continuous through T_c, where K diverges but its coefficient vanishes.
    CHECK_NEAR(onsager_energy(tc * (1 - 1e-9)), -std::sqrt(2.0), 1e-6);
    CHECK_NEAR(onsager_energy(tc * (1 + 1e-9)), -std::sqrt(2.0), 1e-6);
    CHECK_NEAR(onsager_energy(0.3), -2.0, 1e-9);  // ground state
    CHECK_NEAR(onsager_energy(50.0), -2 * std::tanh(1 / 50.0), 2e-4);  // high-T series, leading term
    CHECK(onsager_magnetization(tc + 0.01) == 0.0);
    CHECK_NEAR(onsager_magnetization(1.0), std::pow(1 - std::pow(std::sinh(2.0), -4), 0.125), 1e-15);
    CHECK(onsager_magnetization(tc - 1e-3) > 0.4);  // the 1/8 power rises steeply below T_c
}

int main() {
    lattice_basics();
    density_of_states();
    onsager();
    return finish("test_exact");
}
