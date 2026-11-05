// test_solar_disk_chord -- the area integrals behind illumination.solar_disk_scheme = chord
// (illumination/solar_disk_chord.h). Header-only under test.
//   1. straight horizontal edge vs the analytic segment area (acos p - p sqrt(1 - p^2)) / pi;
//   2. a fully visible disk above the face plane returns cos(i) exactly, uniform and limb-darkened
//      (limb darkening redistributes brightness; it must not change the disk total);
//   3. a disk cut only by the face plane vs ptm::disk::diskMeanCosAboveHorizon;
//   4. random piecewise-linear edges with a tilted, clipped cosine field vs a 3000^2 midpoint grid.
#include "illumination/solar_disk_chord.h"
#include "illumination/primary_disk.h"
#include <cstdio>
#include <random>

using namespace ptm::sdisk;

static int fails = 0;
static void check(bool ok, const char* what, double got, double want, double tol) {
    std::printf("%s %-58s got %.12f want %.12f (tol %.1e)\n", ok ? "PASS" : "FAIL", what, got, want, tol);
    if (!ok) ++fails;
}

static double brute(double A, double Bx, double By, const double* xs, const double* ys, int nc, int side,
                    const LimbLaw& L, int N) {
    double acc = 0.0; const double dA = 4.0 / N / N;
    for (int i = 0; i < N; ++i) {
        const double x = -1 + (i + 0.5) * 2.0 / N;
        for (int j = 0; j < N; ++j) {
            const double y = -1 + (j + 0.5) * 2.0 / N, r2 = x * x + y * y;
            if (r2 >= 1) continue;
            if (nc > 0) {
                double e = ys[0];
                if (nc > 1) {
                    int k = 0; while (k < nc - 2 && x > xs[k + 1]) ++k;
                    const double t = (x - xs[k]) / (xs[k + 1] - xs[k]); e = ys[k] + t * (ys[k + 1] - ys[k]);
                }
                if ((side > 0 && y < e) || (side < 0 && y > e)) continue;
            }
            const double w = A + Bx * x + By * y;
            if (w > 0) acc += w * L(std::sqrt(1 - r2)) * dA;
        }
    }
    return acc / (M_PI * (L.uniform ? 1.0 : L.mean()));
}

int main() {
    LimbLaw lin; lin.uniform = false; lin.a = {{0.4, 0.6, 0, 0, 0, 0}};   // linear u = 0.6
    for (double p : {-0.9, -0.3, 0.0, 0.5, 0.97, 0.995}) {
        const double xs[1] = {0}, ys[1] = {p};
        const double want = (std::acos(p) - p * std::sqrt(1 - p * p)) / M_PI;
        char what[96]; std::snprintf(what, sizeof what, "1. horizontal edge p=%.3f", p);
        const double got = chordFraction(1, 0, 0, xs, ys, 1, +1);
        check(std::fabs(got - want) < 1e-12, what, got, want, 1e-12);
    }
    check(std::fabs(fullDiskFraction(0.7, 0.001, 0.002) - 0.7) < 1e-12, "2. full disk above plane, uniform",
          fullDiskFraction(0.7, 0.001, 0.002), 0.7, 1e-12);
    check(std::fabs(fullDiskFraction(0.7, 0.001, 0.002, lin) - 0.7) < 1e-12, "2. full disk above plane, linear LD",
          fullDiskFraction(0.7, 0.001, 0.002, lin), 0.7, 1e-12);
    const double alpha = 0.00465;
    for (double h_over_a : {-0.8, -0.3, 0.0, 0.4, 0.9}) {
        const double c0 = std::sin(h_over_a * alpha), B = std::tan(alpha) * std::cos(h_over_a * alpha);
        const double got = fullDiskFraction(c0, 0.0, B), want = ptm::disk::diskMeanCosAboveHorizon(c0, alpha);
        char what[96]; std::snprintf(what, sizeof what, "3. plane-only clip h/alpha=%+.1f (rel)", h_over_a);
        check(std::fabs(got - want) <= 2e-3 * want + 1e-12, what, got, want, 2e-3);
    }
    std::mt19937 g(1); std::uniform_real_distribution<double> U(-1, 1);
    double mx = 0, mxl = 0;
    for (int it = 0; it < 20; ++it) {
        const double xs[3] = {-0.5, 0, 0.5}, ys[3] = {U(g) * 0.8, U(g) * 0.8, U(g) * 0.8};
        const double A = U(g) * 0.5 + 0.2, Bx = U(g) * 0.01, By = U(g) * 0.6; const int side = U(g) > 0 ? 1 : -1;
        mx = std::max(mx, std::fabs(chordFraction(A, Bx, By, xs, ys, 3, side) - brute(A, Bx, By, xs, ys, 3, side, uniformLaw(), 3000)));
        mxl = std::max(mxl, std::fabs(chordFraction(A, Bx, By, xs, ys, 3, side, lin) - brute(A, Bx, By, xs, ys, 3, side, lin, 3000)));
    }
    check(mx < 2e-5, "4. random PL edges vs 3000^2 grid, uniform (max abs)", mx, 0.0, 2e-5);
    check(mxl < 2e-5, "4. random PL edges vs 3000^2 grid, linear LD (max abs)", mxl, 0.0, 2e-5);
    std::printf("%s: %d failure(s)\n", fails ? "FAILED" : "PASSED", fails);
    return fails ? 1 : 0;
}
