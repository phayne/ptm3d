#!/usr/bin/env python3
"""The depth-weighted emission model reproduces its reference values exactly.

``scripts/emission/`` was carved out of an older, instrument-specific analysis
directory. The move was meant to be a pure relocation plus a docstring rewrite:
the physics — Maxwell-Garnett/Bruggeman mixing, the Hufford/Maetzler ice
permittivity, Fresnel emissivity and the plane-parallel emission integral — was
not supposed to change by one bit. A rename that quietly perturbs a number is
exactly the failure this file exists to catch, so the comparison below is
``==`` and not ``math.isclose``.

The expected values were generated from the pre-move module and are pasted in
literally. ``repr`` of a Python float round-trips through the parser exactly,
so these are the original doubles, not a rounded copy of them.

    python3 src/python/tests/test_emission_model.py

Runs on the repository's declared numpy floor (>=1.20) and upward.

Deliberate physics changes since the move, each with its own commit, and
each re-pinning only the values it moves:

* 2026-11-12: ``absorption_coefficient`` returned k0 Im(n), the attenuation
  rate of the FIELD amplitude, while the transfer integral needs the
  INTENSITY rate 2 k0 Im(n) (the function's own docstring gave the low-loss
  power form). Every ``kappa`` literal below is exactly twice its previous
  value (doubling is exact in binary), and every ``Tb`` row was regenerated.
  ``physics_checks`` tests the coefficient against the plane-wave intensity
  decay directly, so this file no longer only compares the module with
  itself.
* 2026-11-12: ``radiative_transfer_emission_angle`` stretched the optical
  depth by 1/cos of the OUTSIDE emission angle; the ray inside the regolith
  travels at the refracted angle, sin(theta_t) = sin(theta)/n'. The angled
  ``Tb`` rows were regenerated; ``physics_checks`` brackets the angled result
  between nadir results on depth grids stretched by the extreme refracted
  path factors of the profile.
* 2026-11-12: the depth integral used the trapezoid rule on the
  exponential kernel (weights summing to about 1.014 on the reference grid)
  and dropped the half-space below the column once its optical depth
  reached 5 (an isothermal 100 K column jumped from 100.00 to 99.39 K as
  tau_max crossed 5). Each cell is now integrated exactly and the half-space
  is always included, so an isothermal column returns t T to 1e-12 at any
  depth. Every ``Tb`` row was regenerated; ``physics_checks`` adds the
  isothermal, linear-profile, periodic-profile and thin-skin limits.

The ``==`` above is exact, and exactness is not portable across math-library
backends. Measured: 306/306 reproduce on an MKL-backed numpy 1.26.4 and on
numpy 2.2.6/2.4.3, while an OpenBLAS-backed numpy 1.26.4 and 1.24.4 move 27 of
them by 1-3 ULP. The movement is in the ``dielectric_mixing`` permittivity,
upstream of the emission integral, so a handful of last-bit deltas here is a
backend difference; anything larger, or anything confined to the ``Tb`` rows,
is this module changing.
"""

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "scripts" / "emission"))

import dielectric_mixing as dielectric              # noqa: E402
import depth_weighted_emission as emission          # noqa: E402


DIELECTRIC = [
    (97.5, 80.0, 0.0, 'maxwell_garnett', 3.0126334999999997, 0.0040555958654052, 4.774690262448833, 0.92767978782399),
    (97.5, 80.0, 0.0, 'bruggeman', 3.0126335, 0.004055595865405199, 4.774690262448832, 0.92767978782399),
    (97.5, 80.0, 0.3, 'maxwell_garnett', 2.2973427602771395, 0.0025258829007025244, 3.4053678404589083, 0.9579781447465723),
    (97.5, 80.0, 0.3, 'bruggeman', 2.2704057076945636, 0.0024173297650368773, 3.278293889849833, 0.9591289550910814),
    (97.5, 80.0, 0.5, 'maxwell_garnett', 1.880230431392494, 0.0016788224059870764, 2.5018597251521872, 0.9754919175975098),
    (97.5, 80.0, 0.5, 'bruggeman', 1.827437056660556, 0.0014634428643470496, 2.212169010142664, 0.9776205994375803),
    (97.5, 80.0, 0.7, 'maxwell_garnett', 1.4108496411387301, 0.0005633378079529836, 0.9691529141620197, 0.9926321314473004),
    (97.5, 80.0, 0.7, 'bruggeman', 1.4431435076497312, 0.0006977138405570995, 1.1868240726747172, 0.9916369475793617),
    (97.5, 100.0, 0.0, 'maxwell_garnett', 3.0308335, 0.004161038000038464, 4.884097457598793, 0.9269268386734373),
    (97.5, 100.0, 0.0, 'bruggeman', 3.0308335, 0.004161038000038464, 4.884097457598793, 0.9269268386734373),
    (97.5, 100.0, 0.3, 'maxwell_garnett', 2.3086765802966083, 0.002590914337670014, 3.484457873857908, 0.9574936574056541),
    (97.5, 100.0, 0.3, 'bruggeman', 2.281250701395648, 0.0024787697976014697, 3.3536165580699286, 0.9586657566437566),
    (97.5, 100.0, 0.5, 'maxwell_garnett', 1.8877629319941138, 0.0017218181615273262, 2.5608096241834857, 0.9751857276306481),
    (97.5, 100.0, 0.5, 'bruggeman', 1.8339994729253533, 0.0014992175736153004, 2.2621885813351774, 0.977357742497008),
    (97.5, 100.0, 0.7, 'maxwell_garnett', 1.4133704236255487, 0.0005746581965960823, 0.9877462243725957, 0.9925558886016266),
    (97.5, 100.0, 0.7, 'bruggeman', 1.446269525430922, 0.0007135374003275263, 1.2124278173844487, 0.9915385537515081),
    (97.5, 132.0, 0.0, 'maxwell_garnett', 3.0599535, 0.004345834286688125, 5.076675685921039, 0.9257247376097689),
    (97.5, 132.0, 0.0, 'bruggeman', 3.0599535, 0.0043458342866881244, 5.076675685921038, 0.9257247376097689),
    (97.5, 132.0, 0.3, 'maxwell_garnett', 2.3268049417129055, 0.0027049336386109285, 3.6236006913440324, 0.9567184443579613),
    (97.5, 132.0, 0.3, 'bruggeman', 2.2985899842354813, 0.002586536989893402, 3.486194801892195, 0.9579248356810278),
    (97.5, 132.0, 0.5, 'maxwell_garnett', 1.8998090692930767, 0.0017972197169783728, 2.664464441607973, 0.974694875867538),
    (97.5, 132.0, 0.5, 'bruggeman', 1.8444787959524647, 0.001562047504162584, 2.350288339880908, 0.976936927553129),
    (97.5, 132.0, 0.7, 'maxwell_garnett', 1.4173736030140285, 0.0005946836233533108, 1.0207222466478838, 0.9924342972589911),
    (97.5, 132.0, 0.7, 'bruggeman', 1.4512502085125358, 0.000741398488197141, 1.2576051741235716, 0.9913810583402245),
    (233.0, 80.0, 0.0, 'maxwell_garnett', 3.0126335, 0.009691834221942681, 27.267627407061536, 0.9276795215954695),
    (233.0, 80.0, 0.0, 'bruggeman', 3.0126335, 0.009691834221942681, 27.267627407061536, 0.9276795215954695),
    (233.0, 80.0, 0.3, 'maxwell_garnett', 2.297343090080996, 0.006036211506748229, 19.447609762529115, 0.9579779908599411),
    (233.0, 80.0, 0.3, 'bruggeman', 2.270406432565838, 0.0057767967727520115, 18.7219042630574, 0.9591287946999701),
    (233.0, 80.0, 0.5, 'maxwell_garnett', 1.8802307678925707, 0.004011954067876127, 14.287795117290777, 0.9754918322391761),
    (233.0, 80.0, 0.5, 'bruggeman', 1.827438225408763, 0.003497251189113669, 12.633403333790078, 0.9776204974406351),
    (233.0, 80.0, 0.7, 'maxwell_garnett', 1.4108513499631132, 0.0013462302858940753, 5.534695909060971, 0.9926320718485103),
    (233.0, 80.0, 0.7, 'bruggeman', 1.443144697786011, 0.0016673553501867667, 6.777792085779146, 0.9916368977042535),
    (233.0, 100.0, 0.0, 'maxwell_garnett', 3.0308335000000004, 0.009943813887271405, 27.892436300266837, 0.9269265605654641),
    (233.0, 100.0, 0.0, 'bruggeman', 3.0308335, 0.009943813887271405, 27.892436300266837, 0.9269265605654641),
    (233.0, 100.0, 0.3, 'maxwell_garnett', 2.308676921819879, 0.006191619911773039, 19.899281853425627, 0.957493496400409),
    (233.0, 100.0, 0.3, 'bruggeman', 2.2812514554075607, 0.005923622601564315, 19.15206127254529, 0.9586655889272175),
    (233.0, 100.0, 0.5, 'maxwell_garnett', 1.88776328026632, 0.004114702842306727, 14.624449866302127, 0.9751856382065459),
    (233.0, 100.0, 0.5, 'bruggeman', 1.8340006913778106, 0.003582743442882536, 12.919057809389646, 0.9773576357736834),
    (233.0, 100.0, 0.7, 'maxwell_garnett', 1.4133722069565235, 0.0013732829342035496, 5.640879011389446, 0.9925558262386058),
    (233.0, 100.0, 0.7, 'bruggeman', 1.4462707674947228, 0.0017051694825070672, 6.924011146695767, 0.9915385015135607),
    (233.0, 132.0, 0.0, 'maxwell_garnett', 3.0599535000000007, 0.010385429628700692, 28.992222901716904, 0.9257244379486886),
    (233.0, 132.0, 0.0, 'bruggeman', 3.0599535, 0.010385429628700692, 28.992222901716907, 0.9257244379486886),
    (233.0, 132.0, 0.3, 'maxwell_garnett', 2.326805304654401, 0.006464096703169958, 20.693906844715514, 0.9567182704406261),
    (233.0, 132.0, 0.3, 'bruggeman', 2.298590791270048, 0.006181158380613449, 19.90919656182693, 0.9579246546764332),
    (233.0, 132.0, 0.5, 'maxwell_garnett', 1.8998094390965181, 0.004294893119904759, 15.216408161561723, 0.9746947790687831),
    (233.0, 132.0, 0.5, 'bruggeman', 1.8444801047408985, 0.0037328906265820905, 13.422182570612488, 0.9769368122581751),
    (233.0, 132.0, 0.7, 'maxwell_garnett', 1.417375521604015, 0.0014211382368364953, 5.829199144005945, 0.9924342298862531),
    (233.0, 132.0, 0.7, 'bruggeman', 1.4512515449112924, 0.0017717501175872573, 7.182012030401747, 0.991381001818287),
    (343.5, 80.0, 0.0, 'maxwell_garnett', 3.0126335, 0.01428817620273524, 59.263641362641806, 0.9276791428888831),
    (343.5, 80.0, 0.0, 'bruggeman', 3.0126335, 0.01428817620273524, 59.263641362641806, 0.9276791428888831),
    (343.5, 80.0, 0.3, 'maxwell_garnett', 2.297343559217419, 0.008898876316357166, 42.26757085343754, 0.9579777719581597),
    (343.5, 80.0, 0.3, 'bruggeman', 2.270407463681909, 0.008516433219588247, 40.69030876781773, 0.959128566545732),
    (343.5, 80.0, 0.5, 'maxwell_garnett', 1.8802312465536324, 0.005914616851051246, 31.053199474036866, 0.9754917108180264),
    (343.5, 80.0, 0.5, 'bruggeman', 1.8274398879359086, 0.005155814587669897, 27.457514119607012, 0.9776203523512769),
    (343.5, 80.0, 0.7, 'maxwell_garnett', 1.4108537807462118, 0.001984673040150346, 12.029117484281091, 0.9926319870697414),
    (343.5, 80.0, 0.7, 'bruggeman', 1.4431463907392963, 0.002458093003018815, 14.730885916800895, 0.9916368267574446),
    (343.5, 100.0, 0.0, 'maxwell_garnett', 3.0308334999999995, 0.014659656953981663, 60.62160144730701, 0.9269261649606109),
    (343.5, 100.0, 0.0, 'bruggeman', 3.0308335, 0.014659656953981665, 60.62160144730701, 0.9269261649606109),
    (343.5, 100.0, 0.3, 'maxwell_garnett', 2.3086774076264844, 0.009127986876218774, 43.24923586646469, 0.9574932673725416),
    (343.5, 100.0, 0.3, 'bruggeman', 2.281252527975225, 0.008732890840613442, 41.625212281966846, 0.9586653503527145),
    (343.5, 100.0, 0.5, 'maxwell_garnett', 1.887763775672508, 0.006066094002486244, 31.784886006809334, 0.9751855110019182),
    (343.5, 100.0, 0.5, 'bruggeman', 1.8340024246080264, 0.005281851181835888, 28.078355248834754, 0.9773574839611474),
    (343.5, 100.0, 0.7, 'maxwell_garnett', 1.4133747437238349, 0.0020245550048549805, 12.25989392997592, 0.9925557375277569),
    (343.5, 100.0, 0.7, 'bruggeman', 1.4462725343138059, 0.0025138401991268737, 15.048677676389905, 0.9915384272056527),
    (343.5, 132.0, 0.0, 'maxwell_garnett', 3.0599535, 0.015310708486947133, 63.0118777826606, 0.9257240116850248),
    (343.5, 132.0, 0.0, 'bruggeman', 3.0599535, 0.015310708486947131, 63.0118777826606, 0.9257240116850248),
    (343.5, 132.0, 0.3, 'maxwell_garnett', 2.32680582092713, 0.009529685279399838, 44.97627499053147, 0.9567180230455792),
    (343.5, 132.0, 0.3, 'bruggeman', 2.298591939260597, 0.009112562313450043, 43.270770155544476, 0.9579243971999682),
    (343.5, 132.0, 0.5, 'maxwell_garnett', 1.8998099651294955, 0.0063317391031969235, 33.07144985002486, 0.9746946413738594),
    (343.5, 132.0, 0.5, 'bruggeman', 1.844481966471709, 0.005503205002857689, 29.171845721448815, 0.9769366482526803),
    (343.5, 132.0, 0.7, 'maxwell_garnett', 1.4173782507745856, 0.0020951049359380835, 12.669185237512721, 0.9924341340491535),
    (343.5, 132.0, 0.7, 'bruggeman', 1.451253445919756, 0.0026119961653755857, 15.60941524692299, 0.9913809214164779),
]
PROFILES = [
    (0, 97.5, 0.0, 'fresnel', 'maxwell_garnett', 0.0, 88.39473061520377, 88.39473061520377),
    (1, 233.0, 0.35, 0.75, 'bruggeman', 30.0, 73.74232108880177, 73.83191679968131),
    (2, 343.5, 0.6, 'fresnel', 'maxwell_garnett', 60.0, 99.78332061609036, 97.46451337201222),
    (3, 97.5, 0.0, 0.75, 'bruggeman', 0.0, 74.65977525363212, 74.65977525363212),
    (4, 233.0, 0.35, 'fresnel', 'maxwell_garnett', 30.0, 103.95467373829656, 104.05247551028448),
    (5, 343.5, 0.6, 0.75, 'bruggeman', 60.0, 84.27679403561595, 85.30321937169715),
]

def profile(index):
    """The depth grid and T(z) the reference values were generated from."""
    z = np.concatenate([[0.0], np.geomspace(2e-3, 1.5, 23)])
    return z, 95.0 + (index + 1) * 4.0 * np.exp(-z / (0.02 * (index + 1)))


def physics_checks(failures):
    """Independent of the module's own numbers: the absorption coefficient is
    the e-folding rate of a plane wave's intensity |exp(i k0 n z)|^2, and it
    matches the low-loss form 2 pi nu sqrt(eps') tan(delta) / c."""
    for freq, temp, por, mixing, *_ in DIELECTRIC:
        eps = complex(dielectric.porous_ice_permittivity(freq, temp, por, mixing))
        k0 = 2.0 * np.pi * freq * 1e9 / dielectric.C_LIGHT
        n = np.sqrt(eps)
        dz = 0.01
        alpha = np.log(abs(np.exp(1j * k0 * n * 0.0)) ** 2
                       / abs(np.exp(1j * k0 * n * dz)) ** 2) / dz
        low_loss = k0 * np.sqrt(eps.real) * eps.imag / eps.real
        kappa = float(dielectric.absorption_coefficient(freq, temp, por, mixing))
        tag = f"{freq:g} GHz, {temp:g} K, porosity {por:g}, {mixing}"
        if abs(kappa / alpha - 1.0) > 1e-9:
            failures.append(f"kappa vs plane-wave intensity decay [{tag}]: "
                            f"{kappa!r} vs {alpha!r}")
        if abs(kappa / low_loss - 1.0) > 1e-4:
            failures.append(f"kappa vs low-loss formula [{tag}]: "
                            f"{kappa!r} vs {low_loss!r}")

    # Refraction. With a fixed emissivity and depth-independent properties,
    # substituting u = z / cos(theta_t) turns the slant integral into a nadir
    # integral of the same profile on depths stretched by 1/cos(theta_t). Here
    # n' varies slightly with temperature, so the angled result must lie
    # between the nadir results on grids stretched by the smallest and largest
    # 1/cos(theta_t) in the column. An outside-angle law (factor
    # 1/cos(theta) = 2 at 60 deg) lands far outside that bracket.
    n_bracket = 0
    for index, freq, por, emis, mixing, angle, *_ in PROFILES:
        if angle == 0.0:
            continue
        z, temperature = profile(index)
        n_re = np.real(np.sqrt(np.array([
            dielectric.porous_ice_permittivity(freq, t, por, mixing)
            for t in temperature])))
        s2 = np.sin(np.radians(angle)) ** 2
        stretch = 1.0 / np.sqrt(1.0 - s2 / n_re ** 2)
        a = emission.radiative_transfer_1d(z * stretch.max(), temperature,
                                           freq, por, 0.9, mixing)
        b = emission.radiative_transfer_1d(z * stretch.min(), temperature,
                                           freq, por, 0.9, mixing)
        lo, hi = min(a, b), max(a, b)
        tb = emission.radiative_transfer_emission_angle(
            z, temperature, freq, por, angle, 0.9, mixing)
        if not (lo - 1e-9 <= tb <= hi + 1e-9):
            failures.append(f"refracted path [profile {index} at {angle:g} deg]:"
                            f" {tb!r} outside [{lo!r}, {hi!r}]")
        n_bracket += 1

    # The Planck integral (brightness_temperature). Its refracted path obeys
    # the same bracket, and its Rayleigh-Jeans brightness temperature sits
    # below the temperature-weighted result by t (h nu/2k - (h nu/k)^2/12T),
    # the Planck curvature, which a temperature-weighted integral cannot show.
    n_planck = 0
    for index, freq, por, emis, mixing, angle, *_ in PROFILES:
        z, temperature = profile(index)
        bt = emission.brightness_temperature(z, temperature, freq, por, angle,
                                             0.9, mixing)
        tw = emission.radiative_transfer_emission_angle(
            z, temperature, freq, por, angle, 0.9, mixing)
        x = dielectric.H_PLANCK * freq * 1e9 / dielectric.K_BOLTZ
        expected = 0.9 * (x / 2.0 - x ** 2 / (12.0 * tw / 0.9))
        # A temperature-weighted integral would give an offset of zero.
        if abs((tw - bt.rayleigh_jeans) / expected - 1.0) > 1e-3:
            failures.append(f"Planck integral [profile {index}]: T-weighted - RJ "
                            f"= {tw - bt.rayleigh_jeans!r}, expected "
                            f"{expected!r}")
        if not bt.planck > bt.rayleigh_jeans:
            failures.append(f"Planck integral [profile {index}]: Planck Tb "
                            f"{bt.planck!r} not above RJ {bt.rayleigh_jeans!r}")
        n_planck += 2

    # The depth integral itself, against closed forms.
    import warnings
    n_col = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")   # shallow columns are tested on purpose
        # Isothermal: t T exactly for any column optical depth, both forms.
        kap = float(dielectric.absorption_coefficient(233.0, 100.0, 0.5))
        for tau_max in (0.01, 0.5, 4.9, 5.1, 7.0, 50.0):
            z = np.geomspace(1e-4, tau_max / kap, 40)
            z = np.concatenate([[0.0], z])
            T = np.full_like(z, 100.0)
            tb = emission.radiative_transfer_1d(z, T, 233.0, 0.5, 0.9)
            bt = emission.brightness_temperature(z, T, 233.0, 0.5, 45.0, 0.9)
            rj_exact = 0.9 * float(
                dielectric.C_LIGHT ** 2 * emission.planck_nu(100.0, 233e9)
                / (2 * dielectric.K_BOLTZ * 233e9 ** 2))
            if abs(tb - 90.0) > 1e-12 * 100 or abs(bt.rayleigh_jeans - rj_exact) > 1e-12 * 100:
                failures.append(f"isothermal column, tau_max {tau_max}: "
                                f"{tb!r} and RJ {bt.rayleigh_jeans!r}, expected "
                                f"90.0 and {rj_exact!r}")
            n_col += 1
        # Linear profile, constant kappa, column to depth Z with the half-space
        # below at T(Z): T0 + g (1 - exp(-kappa Z)) / kappa.
        k, T0, g = 14.6, 90.0, 40.0
        for Z in (0.05, 0.3, 2.0):
            z = np.concatenate([[0.0], np.geomspace(1e-4, Z, 30)])
            w, wb, _ = emission.column_weights(z, np.full_like(z, k))
            got = np.dot(w, T0 + g * z) + wb * (T0 + g * Z)
            want = T0 + g * -np.expm1(-k * Z) / k
            if abs(got - want) > 1e-12 * want:
                failures.append(f"linear profile, Z {Z}: {got!r} vs {want!r}")
            n_col += 1
        # Periodic conduction profile T0 + A e^{-z/l} cos(wt - z/l): response
        # amplitude 1/sqrt((1+b)^2 + b^2), lag atan(b/(1+b)), b = L_e/l.
        l, A = 0.05, 20.0
        for L_e in (0.01, 0.0684, 0.2):
            k = 1.0 / L_e
            z = np.linspace(0.0, 20 * max(l, L_e), 40001)
            phase = np.linspace(0, 2 * np.pi, 73)[:-1]
            w, wb, _ = emission.column_weights(z, np.full_like(z, k))
            resp = np.array([np.dot(w, A * np.exp(-z / l) * np.cos(p - z / l))
                             + wb * A * np.exp(-z[-1] / l) * np.cos(p - z[-1] / l)
                             for p in phase])
            c = np.sum(resp * np.exp(-1j * phase)) * 2 / len(phase)
            b = L_e / l
            amp, lag = abs(c) / A, -np.angle(c)
            if (abs(amp - 1 / np.hypot(1 + b, b)) > 1e-6
                    or abs(lag - np.arctan(b / (1 + b))) > 1e-6):
                failures.append(f"periodic profile, L_e/l {b:.3f}: amplitude "
                                f"{amp!r}, lag {lag!r}")
            n_col += 1
        # Thin electrical skin (L_e = 1 um) recovers the surface: for
        # 100 + 30 exp(-z/d) the exact result is 100 + 30 k d / (k d + 1).
        z = np.concatenate([[0.0], np.geomspace(1e-9, 1.0, 400)])
        k, d = 1e6, 0.05
        w, wb, _ = emission.column_weights(z, np.full_like(z, k))
        f = 100.0 + 30.0 * np.exp(-z / d)
        got = np.dot(w, f) + wb * f[-1]
        if abs(got - (100.0 + 30.0 * k * d / (k * d + 1))) > 1e-6:
            failures.append(f"thin skin: {got!r} vs "
                            f"{100.0 + 30.0 * k * d / (k * d + 1)!r}")
        n_col += 1
    return 2 * len(DIELECTRIC) + n_bracket + n_planck + n_col


def main():
    failures = []

    def check(label, got, expected):
        if got != expected:
            failures.append(f"{label}: got {got!r}, expected {expected!r} "
                            f"(delta {got - expected!r})")

    for freq, temp, por, mixing, eps_re, eps_im, kappa, emis in DIELECTRIC:
        eps = dielectric.porous_ice_permittivity(freq, temp, por, mixing)
        tag = f"{freq:g} GHz, {temp:g} K, porosity {por:g}, {mixing}"
        check(f"Re(eps) [{tag}]", float(np.real(eps)), eps_re)
        check(f"Im(eps) [{tag}]", float(np.imag(eps)), eps_im)
        check(f"kappa [{tag}]",
              float(dielectric.absorption_coefficient(freq, temp, por, mixing)),
              kappa)
        check(f"emissivity [{tag}]",
              float(dielectric.fresnel_emissivity(freq, temp, por, 0.0, mixing)),
              emis)

    for index, freq, por, emis, mixing, angle, tb_nadir, tb_angle in PROFILES:
        z, temperature = profile(index)
        tag = f"profile {index} at {freq:g} GHz"
        check(f"Tb nadir [{tag}]",
              emission.radiative_transfer_1d(z, temperature, freq, por, emis,
                                             mixing),
              tb_nadir)
        check(f"Tb at {angle:g} deg [{tag}]",
              emission.radiative_transfer_emission_angle(
                  z, temperature, freq, por, angle, emis, mixing),
              tb_angle)

    # radiative_transfer_batch must agree with the scalar path face by face,
    # which is the only property of it the callers rely on.
    z, temperature = profile(2)
    stack = np.vstack([profile(i)[1] for i in range(6)])
    batch = emission.radiative_transfer_batch(z, stack, 233.0, 0.4, 0.9,
                                              "maxwell_garnett")
    for i in range(stack.shape[0]):
        check(f"batch face {i} vs scalar",
              float(batch[i]),
              emission.radiative_transfer_1d(z, stack[i], 233.0, 0.4, 0.9,
                                             "maxwell_garnett"))

    n_checks = 4 * len(DIELECTRIC) + 2 * len(PROFILES) + stack.shape[0]
    n_checks += physics_checks(failures)
    if failures:
        print(f"[test_emission_model] FAIL -- {len(failures)} of {n_checks} "
              f"values moved:")
        for line in failures:
            print(f"  {line}")
        return 1

    print(f"[test_emission_model] PASS -- {n_checks} checks: values reproduced "
          f"exactly ({len(DIELECTRIC)} dielectric states, {len(PROFILES)} "
          f"profiles) and the absorption coefficient matches the plane-wave "
          f"intensity decay")
    return 0


if __name__ == "__main__":
    sys.exit(main())
