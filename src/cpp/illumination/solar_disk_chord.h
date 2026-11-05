// solar_disk_chord.h -- area integrals for the probe-then-chord solar-disk rule.
//
// Header-only geometry, no rays (illumination_model.cpp traces them). The disk is the unit disk in a
// limb-aligned frame: x along u, y along v, where v is the local vertical projected onto the plane of
// the sky (y = +1 is the top limb, toward the zenith). A direction (x, y) on the disk is
// normalize(s + tan(alpha) * (x u + y v)). The face's incidence cosine across the disk is, to O(alpha^2),
// the linear field w(x, y) = A + Bx x + By y with A = n.s, Bx = tan(alpha) n.u, By = tan(alpha) n.v.
//
// chordFraction() integrates max(0, w) * L(mu) over the part of the disk on one side of a
// piecewise-linear edge, divided by pi * <L>. So a fully visible disk entirely above the face plane
// returns A = cos(i), and flux = S * chordFraction(...). L is the limb-darkening law (uniform: L = 1).
//
// The integral is evaluated piece by piece in x: the x range is cut at every point where the active
// lower or upper bound can change (PL knots, edge/circle, weight-line/circle and edge/weight-line
// crossings), so on each piece the bounds are fixed smooth functions. With x = sin(theta) the disk
// boundary sqrt(1 - x^2) = cos(theta) is smooth too, and 16-point Gauss-Legendre in theta (and, for a
// limb-darkened disk, in phi with y = sqrt(1 - x^2) sin(phi)) is exact to ~1e-12 for these integrands.
// No sampling: the result does not depend on a sample count.
//
// Fermi solar-disk spec 2026-11-05, section 2 (~/projects/solar_disk_spec_2026-11/SPEC.md).
#ifndef PTM_SOLAR_DISK_CHORD_H
#define PTM_SOLAR_DISK_CHORD_H

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace ptm {
namespace sdisk {

// Limb-darkening law I(mu)/I(1) = sum_k a[k] mu^k, k = 0..5 (Neckel & Labs 1994 form, eq. 5).
// mean = <I>/I(1) over the disk = 2 sum_k a[k]/(k+2) (their eq. 6, F/I).
struct LimbLaw {
    std::array<double, 6> a{{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    bool uniform = true;
    double mean() const {
        double m = 0.0;
        for (int k = 0; k < 6; ++k) m += 2.0 * a[k] / (k + 2.0);
        return m;
    }
    double operator()(double mu) const {
        double v = a[5];
        for (int k = 4; k >= 0; --k) v = v * mu + a[k];
        return v;
    }
};

inline const LimbLaw& uniformLaw() {
    static const LimbLaw u{};
    return u;
}

namespace detail {
// 16-point Gauss-Legendre on [-1, 1].
inline const std::array<double, 16>& glX() {
    static const std::array<double, 16> x{{
        -0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
        -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.0950125098376374,
         0.0950125098376374,  0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
         0.7554044083550030,  0.8656312023878318,  0.9445750230732326,  0.9894009349916499}};
    return x;
}
inline const std::array<double, 16>& glW() {
    static const std::array<double, 16> w{{
        0.0271524594117541, 0.0622535239386479, 0.0951585116824928, 0.1246289712555339,
        0.1495959888165767, 0.1691565193950025, 0.1826034150449236, 0.1894506104550685,
        0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
        0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541}};
    return w;
}

struct Line { double p, q; };   // y = p + q x

// x in [-1, 1] where the line meets the unit circle (0 or 2 roots appended).
inline void lineCircle(const Line& L, std::vector<double>& out) {
    // x^2 + (p + q x)^2 = 1  ->  (1 + q^2) x^2 + 2 p q x + p^2 - 1 = 0
    const double a = 1.0 + L.q * L.q, b = 2.0 * L.p * L.q, c = L.p * L.p - 1.0;
    const double disc = b * b - 4.0 * a * c;
    if (disc <= 0.0) return;
    const double sq = std::sqrt(disc);
    out.push_back((-b - sq) / (2.0 * a));
    out.push_back((-b + sq) / (2.0 * a));
}
inline void lineLine(const Line& L1, const Line& L2, std::vector<double>& out) {
    const double dq = L1.q - L2.q;
    if (std::fabs(dq) > 1e-300) out.push_back((L2.p - L1.p) / dq);
}
}  // namespace detail

// Integral over the unit disk of max(0, A + Bx x + By y) * L(mu) / (pi <L>), restricted to y >= e(x)
// (side = +1) or y <= e(x) (side = -1). e is piecewise linear through (xs[k], ys[k]), k < nc, xs
// ascending, extended linearly beyond the end knots. nc = 0 means no edge (the whole disk).
inline double chordFraction(double A, double Bx, double By, const double* xs, const double* ys, int nc,
                            int side, const LimbLaw& law = uniformLaw()) {
    using detail::Line;
    // Edge segments: segment k covers (knot k, knot k+1], the first and last extend to -/+ infinity.
    std::vector<Line> seg;
    std::vector<double> knots;   // interior knots, the x where the active segment changes
    if (nc == 1) {
        seg.push_back(Line{ys[0], 0.0});
    } else if (nc >= 2) {
        for (int k = 0; k + 1 < nc; ++k) {
            const double q = (ys[k + 1] - ys[k]) / (xs[k + 1] - xs[k]);
            seg.push_back(Line{ys[k] - q * xs[k], q});
            if (k > 0) knots.push_back(xs[k]);
        }
    }
    auto edgeAt = [&](double x) {
        size_t k = 0;
        while (k < knots.size() && x > knots[k]) ++k;
        return seg[k].p + seg[k].q * x;
    };
    // Weight zero line A + Bx x + By y = 0.
    const bool have_wline = std::fabs(By) > 1e-300;
    const Line wline = have_wline ? Line{-A / By, -Bx / By} : Line{0.0, 0.0};

    std::vector<double> cut{-1.0, 1.0};
    for (double k : knots) cut.push_back(k);
    for (const Line& s : seg) detail::lineCircle(s, cut);
    if (have_wline) {
        detail::lineCircle(wline, cut);
        for (const Line& s : seg) detail::lineLine(s, wline, cut);
    } else if (std::fabs(Bx) > 1e-300) {
        cut.push_back(-A / Bx);
    }
    for (double& c : cut) c = std::max(-1.0, std::min(1.0, c));
    std::sort(cut.begin(), cut.end());

    const auto& gx = detail::glX();
    const auto& gw = detail::glW();
    double acc = 0.0;
    for (size_t j = 0; j + 1 < cut.size(); ++j) {
        const double ta = std::asin(cut[j]), tb = std::asin(cut[j + 1]);
        if (!(tb - ta > 1e-15)) continue;
        const double hm = 0.5 * (tb - ta), hc = 0.5 * (tb + ta);
        for (int g = 0; g < 16; ++g) {
            const double th = hc + hm * gx[g];
            const double x = std::sin(th), cth = std::cos(th);   // dx = cos(theta) dtheta
            double lo = -cth, hi = cth;
            if (!seg.empty()) {
                const double e = edgeAt(x);
                if (side > 0) lo = std::max(lo, e); else hi = std::min(hi, e);
            }
            const double a = A + Bx * x;
            if (have_wline) {
                const double y0 = -a / By;
                if (By > 0.0) lo = std::max(lo, y0); else hi = std::min(hi, y0);
            } else if (a <= 0.0) {
                continue;
            }
            if (!(hi > lo)) continue;
            double inner;
            if (law.uniform) {
                inner = a * (hi - lo) + 0.5 * By * (hi * hi - lo * lo);
            } else {
                // y = c sin(phi), mu = c cos(phi), dy = c cos(phi) dphi, c = cos(theta)
                const double c = cth;
                if (!(c > 0.0)) continue;
                const double pa = std::asin(std::max(-1.0, std::min(1.0, lo / c)));
                const double pb = std::asin(std::max(-1.0, std::min(1.0, hi / c)));
                const double pm = 0.5 * (pb - pa), pc = 0.5 * (pb + pa);
                inner = 0.0;
                for (int h = 0; h < 16; ++h) {
                    const double ph = pc + pm * gx[h];
                    const double y = c * std::sin(ph), cp = c * std::cos(ph);
                    inner += gw[h] * std::max(0.0, a + By * y) * law(cp) * cp;
                }
                inner *= pm;
            }
            acc += gw[g] * hm * inner * cth;
        }
    }
    const double norm = law.uniform ? M_PI : M_PI * law.mean();
    return acc / norm;
}

// The whole disk, clipped by the face plane only (no terrain edge).
inline double fullDiskFraction(double A, double Bx, double By, const LimbLaw& law = uniformLaw()) {
    return chordFraction(A, Bx, By, nullptr, nullptr, 0, +1, law);
}

}  // namespace sdisk
}  // namespace ptm

#endif  // PTM_SOLAR_DISK_CHORD_H
