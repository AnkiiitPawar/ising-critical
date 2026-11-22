// The error-analysis tools, checked on processes with known answers.
#include <random>

#include "check.hpp"
#include "ising/stats.hpp"

using namespace ising;

static void fft_round_trip() {
    std::vector<std::complex<double>> a(64), b;
    for (int i = 0; i < 64; ++i) a[i] = {std::sin(0.3 * i), std::cos(1.7 * i)};
    b = a;
    fft(b, false);
    fft(b, true);
    for (int i = 0; i < 64; ++i) CHECK_NEAR(std::abs(a[i] - b[i]), 0.0, 1e-12);
}

// AR(1): x_{t+1} = rho x_t + sqrt(1 - rho^2) noise has autocorrelation rho^t,
// so tau_int = 1/2 + sum_{t>=1} rho^t = (1 + rho) / (2 (1 - rho)).
static void ar1_autocorrelation_time() {
    for (const double rho : {0.0, 0.5, 0.8, 0.95}) {
        std::mt19937_64 rng(42);
        std::normal_distribution<double> noise;
        std::vector<double> x(1 << 20);
        double v = 0;
        for (auto& xi : x) {
            v = rho * v + std::sqrt(1 - rho * rho) * noise(rng);
            xi = v;
        }
        const auto r = autocorrelation(x, 10);
        CHECK_NEAR(r[0], 1.0, 1e-12);
        CHECK_NEAR(r[1], rho, 0.01);
        const auto tau = integrated_autocorrelation_time(x);
        const double exact = (1 + rho) / (2 * (1 - rho));
        CHECK_NEAR(tau.tau, exact, 5 * tau.error + 1e-3);
        CHECK(tau.error < 0.03 * exact + 1e-3);
    }
}

static void jackknife_of_independent_samples() {
    std::mt19937_64 rng(5);
    std::normal_distribution<double> noise;
    std::vector<Sample> xs(100000);
    for (auto& s : xs) s = {noise(rng), 0.0};
    const auto est = jackknife(xs, 100, [](const Moments& m) { return m[0]; });
    CHECK_NEAR(est.value, 0.0, 5 * est.error);
    CHECK_NEAR(est.error, 1 / std::sqrt(100000.0), 0.25 / std::sqrt(100000.0));
}

int main() {
    fft_round_trip();
    ar1_autocorrelation_time();
    jackknife_of_independent_samples();
    return finish("test_stats");
}
