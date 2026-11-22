// Error analysis for correlated Monte Carlo time series.
//
// Successive measurements from a Markov chain are correlated, so the naive
// standard error is too small. Two tools:
//  - the jackknife over contiguous bins, for error bars on any function of
//    averages (specific heat, susceptibility, Binder cumulant);
//  - the integrated autocorrelation time tau_int, which measures how many
//    steps one independent sample costs.
#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <vector>

namespace ising {

/// One measurement: e = E/N and m = M/N.
struct Sample {
    double e, m;
};

/// The five moments the observables are built from.
using Moments = std::array<double, 5>;  // <e>, <e^2>, <|m|>, <m^2>, <m^4>

inline Moments moments_of(const Sample* begin, const Sample* end) {
    Moments mo{};
    for (auto p = begin; p != end; ++p) {
        const double m2 = p->m * p->m;
        mo[0] += p->e;
        mo[1] += p->e * p->e;
        mo[2] += std::fabs(p->m);
        mo[3] += m2;
        mo[4] += m2 * m2;
    }
    const double n = static_cast<double>(end - begin);
    for (auto& x : mo) x /= n;
    return mo;
}

struct Estimate {
    double value, error;
};

/// Jackknife over `bins` contiguous blocks. `f` maps the five moments to an
/// observable. Bins should be much longer than the autocorrelation time.
inline Estimate jackknife(const std::vector<Sample>& xs, int bins, const std::function<double(const Moments&)>& f) {
    const std::size_t per = xs.size() / bins;
    std::vector<Moments> bin(bins);
    for (int b = 0; b < bins; ++b) bin[b] = moments_of(xs.data() + b * per, xs.data() + (b + 1) * per);
    Moments total{};
    for (const auto& m : bin)
        for (int k = 0; k < 5; ++k) total[k] += m[k] / bins;
    std::vector<double> loo(bins);
    double mean_loo = 0;
    for (int b = 0; b < bins; ++b) {
        Moments m{};
        for (int k = 0; k < 5; ++k) m[k] = (total[k] * bins - bin[b][k]) / (bins - 1);
        loo[b] = f(m);
        mean_loo += loo[b] / bins;
    }
    double var = 0;
    for (double v : loo) var += (v - mean_loo) * (v - mean_loo);
    var *= static_cast<double>(bins - 1) / bins;
    return {f(total), std::sqrt(var)};
}

/// In-place radix-2 FFT; the length must be a power of two.
inline void fft(std::vector<std::complex<double>>& a, bool inverse) {
    const std::size_t n = a.size();
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const double ang = 2 * 3.14159265358979323846 / static_cast<double>(len) * (inverse ? 1 : -1);
        const std::complex<double> wl(std::cos(ang), std::sin(ang));
        for (std::size_t i = 0; i < n; i += len) {
            std::complex<double> w(1);
            for (std::size_t k = 0; k < len / 2; ++k) {
                const auto u = a[i + k], v = a[i + k + len / 2] * w;
                a[i + k] = u + v;
                a[i + k + len / 2] = u - v;
                w *= wl;
            }
        }
    }
    if (inverse)
        for (auto& x : a) x /= static_cast<double>(n);
}

/// Normalised autocorrelation function rho(t) for t = 0 .. max_lag, via FFT.
inline std::vector<double> autocorrelation(const std::vector<double>& x, std::size_t max_lag) {
    const std::size_t n = x.size();
    double mean = 0;
    for (double v : x) mean += v / static_cast<double>(n);
    std::size_t size = 1;
    while (size < 2 * n) size <<= 1;  // zero padding avoids wrap-around
    std::vector<std::complex<double>> a(size);
    for (std::size_t i = 0; i < n; ++i) a[i] = x[i] - mean;
    fft(a, false);
    for (auto& v : a) v = std::norm(v);
    fft(a, true);
    std::vector<double> rho(std::min(max_lag, n - 1) + 1);
    const double c0 = a[0].real() / static_cast<double>(n);
    for (std::size_t t = 0; t < rho.size(); ++t) rho[t] = (a[t].real() / static_cast<double>(n - t)) / c0;
    return rho;
}

/// Integrated autocorrelation time tau_int = 1/2 + sum_{t>=1} rho(t), summed
/// up to the first window W with W >= c * tau_int(W) (Sokal's automatic
/// windowing). The variance of a mean over n samples is then
/// 2 tau_int sigma^2 / n. Returns {tau_int, statistical error, window}.
struct Tau {
    double tau, error;
    std::size_t window;
};

inline Tau integrated_autocorrelation_time(const std::vector<double>& x, double c = 6.0) {
    const auto rho = autocorrelation(x, x.size() / 2);
    double tau = 0.5;
    std::size_t w = 1;
    for (; w < rho.size(); ++w) {
        tau += rho[w];
        if (static_cast<double>(w) >= c * tau) break;
    }
    // Madras-Sokal estimate of the error of tau_int.
    const double err = tau * std::sqrt(2.0 * (2.0 * static_cast<double>(w) + 1.0) / static_cast<double>(x.size()));
    return {tau, err, w};
}

}  // namespace ising
