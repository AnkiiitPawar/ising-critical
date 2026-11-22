// A minimal test harness: no dependencies.
#pragma once

#include <cmath>
#include <cstdio>

inline int g_failures = 0;
inline int g_checks = 0;

#define CHECK(cond)                                                              \
    do {                                                                         \
        ++g_checks;                                                              \
        if (!(cond)) {                                                           \
            std::fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
            ++g_failures;                                                        \
        }                                                                        \
    } while (0)

#define CHECK_NEAR(a, b, tol)                                                    \
    do {                                                                         \
        ++g_checks;                                                              \
        const double va_ = (a), vb_ = (b), vt_ = (tol);                          \
        if (!(std::fabs(va_ - vb_) <= vt_)) {                                    \
            std::fprintf(stderr, "%s:%d: |%s - %s| = |%.10g - %.10g| > %.3g\n", __FILE__, __LINE__, #a, #b, va_, vb_, vt_); \
            ++g_failures;                                                        \
        }                                                                        \
    } while (0)

inline int finish(const char* name) {
    if (g_failures) {
        std::fprintf(stderr, "%s: %d of %d checks failed\n", name, g_failures, g_checks);
        return 1;
    }
    std::printf("%s: all %d checks passed\n", name, g_checks);
    return 0;
}
