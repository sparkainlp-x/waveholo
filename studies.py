#!/usr/bin/env python3
"""v0.1.3 studies: every table and figure in RESULTS_v0.1.3.md.

Run:  python3 studies.py [--quick] [--outdir results]
All randomness is seeded; the output text is written to
results/studies_output.txt and figures to results/*.png.

S1  Reference piston: imag_energy depends on it, realness does not.
S2  The v0.1.2 "comb" is the lambda/2 piston oscillation, aliased.
S3  Sub-bin carrier estimation vs integer rounding.
S4  Residual realness bias = twin/DC cross-talk (window, bandwidth).
S5  Off-grid focus: coarse sweep + golden-section refinement.
S6  Noise with a random reference piston (honest robustness).
S7  Phase and mixed objects: ac_realness vs the alternatives.
S8  Two-plane sample: per-tile focus map recovers both planes.
S9  Spectral identity: realness depends only on F(k) F(-k).
S10 Real fly hologram (pyDHM): metrics side by side.
"""

import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wave_hologram as W  # noqa: E402
from noise_robustness import add_noise  # noqa: E402

N = W.N
OUT = None
LOG = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def fig_ax(nrows=1, ncols=1, size=(8, 4)):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure
    fig = Figure(figsize=size, dpi=110)
    axes = fig.subplots(nrows, ncols)
    return fig, axes


def save(fig, name):
    fig.tight_layout()
    path = os.path.join(OUT, name)
    fig.savefig(path)
    say(f"  figure: {path}")


def holo_field(obj, z, car=W.CAR_DEFAULT, ref_phase=0.0, win=W.WIN_DEFAULT,
               noise=None, rng=None, carrier_y=0.0, estimate=True):
    sensor = W.angular_spectrum(np.asarray(obj, complex), z)
    h = W.record_hologram(sensor, car, carrier_y=carrier_y, ref_phase=ref_phase)
    if noise is not None:
        h = add_noise(h, noise[0], noise[1], rng)
    carrier = (W.estimate_carrier_2d(h, dc_radius=30.0 / N) if estimate
               else (-float(carrier_y), -float(car)))
    return W.reconstruct_field(h, win=win, carrier=carrier), carrier


COARSE = np.linspace(-10e-3, 0.0, 101)


# ---------------------------------------------------------------- S1
def s1_piston():
    say("\n== S1. REFERENCE PISTON ==")
    say("Truth -5.00 mm, clean. The reference phase is a free constant in")
    say("any real setup; v0.1.2 always simulated it as 0.")
    obj = W.make_object()
    say(f"{'ref phase':>10} {'imag@truth':>10} {'real@truth':>10} "
        f"{'imag focus':>10} {'real focus':>10}")
    phases = np.linspace(0, np.pi, 9)
    rows = []
    for ph in phases:
        f0, _ = holo_field(obj, 5e-3, ref_phase=ph)
        at = W.angular_spectrum(f0, -5e-3)
        ri = W.find_focus(f0, COARSE, metric="imag_energy")
        rr = W.find_focus(f0, COARSE, metric="realness")
        rows.append((ph, W.imag_energy(at), W.realness(at), ri["z"], rr["z"]))
        say(f"{ph:>10.3f} {W.imag_energy(at):>10.4f} {W.realness(at):>10.4f} "
            f"{ri['z'] * 1e3:>+10.3f} {rr['z'] * 1e3:>+10.3f}")
    rows = np.array(rows)
    say(f"imag_energy at truth spans {rows[:, 1].min():.4f}..{rows[:, 1].max():.4f}"
        f" (sin^2 of the piston); realness spans "
        f"{rows[:, 2].min():.6f}..{rows[:, 2].max():.6f}.")
    say(f"worst imag focus error {np.max(np.abs(rows[:, 3] + 5e-3)) * 1e3:.3f} mm;"
        f" worst realness error {np.max(np.abs(rows[:, 4] + 5e-3)) * 1e3:.3f} mm")
    fig, ax = fig_ax()
    ax.plot(rows[:, 0], rows[:, 1], "o-", label="imag_energy (v0.1.2)")
    ax.plot(rows[:, 0], rows[:, 2], "s-", label="realness (v0.1.3)")
    ax.plot(phases, np.sin(phases) ** 2 * (1 - rows[0, 1]) + rows[0, 1] *
            np.cos(phases) ** 2, "k:", lw=1, label="sin^2 prediction")
    ax.set_xlabel("reference piston [rad]")
    ax.set_ylabel("metric value AT TRUE FOCUS")
    ax.set_title("At true focus, imag_energy is whatever the piston makes it")
    ax.legend()
    save(fig, "s1_piston.png")


# ---------------------------------------------------------------- S2
def s2_comb(quick):
    say("\n== S2. THE 'COMB' IS THE lambda/2 PISTON OSCILLATION, ALIASED ==")
    obj = W.make_object()
    f0, _ = holo_field(obj, 5.03e-3)
    prop = W.Propagator(f0)
    fine = np.linspace(-5.04e-3, -5.02e-3, 801 if not quick else 201)
    ie = np.array([W.imag_energy(prop(z)) for z in fine])
    rl = np.array([W.realness(prop(z)) for z in fine])
    mins = [i for i in range(1, len(ie) - 1) if ie[i] < ie[i - 1] and ie[i] < ie[i + 1]]
    period = float(np.median(np.diff(fine[mins]))) if len(mins) > 2 else np.nan
    say(f"imag_energy minima spacing on a {np.diff(fine)[0] * 1e9:.0f} nm grid: "
        f"{period * 1e6:.4f} um (lambda/2 = {W.WL / 2 * 1e6:.4f} um)")
    alias = 1e-6 - 4 * (W.WL / 2)
    say(f"sampling at 1 um aliases a {W.WL / 2 * 1e6:.3f} um period to "
        f"1/|1/0.266 - 4|: apparent period "
        f"{1 / abs(1 / (W.WL / 2 * 1e6) - round(1 / (W.WL / 2 * 1e6))):.2f} um"
        f" (v0.1.2 reported teeth every 4 um); residual {alias * 1e9:.1f} nm")
    say(f"realness over the same 20 um: {rl.min():.5f}..{rl.max():.5f} "
        f"(variation {100 * (rl.max() - rl.min()) / rl.min():.2f}%) - smooth")
    coarse1 = np.arange(-5.04e-3, -5.02e-3 + 1e-12, 1e-6)
    ie1 = np.array([W.imag_energy(prop(z)) for z in coarse1])
    fig, ax = fig_ax()
    ax.plot(fine * 1e3, ie, lw=0.8, label="imag_energy, 25 nm steps")
    ax.plot(coarse1 * 1e3, ie1, "o-", ms=3, label="imag_energy, 1 um steps (v0.1.2 'comb')")
    ax.plot(fine * 1e3, rl, lw=2, label="realness (piston-invariant)")
    ax.set_xlabel("reconstruction distance [mm]")
    ax.set_ylabel("metric")
    ax.set_title("No comb: a lambda/2 oscillation sampled at 1 um")
    ax.legend(fontsize=8)
    save(fig, "s2_no_comb.png")


# ---------------------------------------------------------------- S3
def s3_carrier():
    say("\n== S3. SUB-BIN CARRIER ==")
    obj = W.make_object()
    sensor = W.angular_spectrum(obj.astype(complex), 5e-3)
    say(f"{'true car':>8} {'int est':>7} {'subbin est':>10} {'err':>7} "
        f"{'imag int':>8} {'imag sub':>8} {'real int':>8} {'real sub':>8}  (focus err, mm)")
    errs = []
    for car in np.arange(80.0, 81.01, 0.1):
        h = W.record_hologram(sensor, car)
        ci = W.estimate_carrier(h, subbin=False)
        cs = W.estimate_carrier_2d(h, dc_radius=30.0 / N)
        out = []
        for carrier in ((0.0, -float(ci)), cs):
            f0 = W.reconstruct_field(h, carrier=carrier)
            for m in ("imag_energy", "realness"):
                out.append(abs(W.find_focus(f0, COARSE, metric=m)["z"] + 5e-3) * 1e3)
        errs.append((abs(-cs[1] - car), *out))
        say(f"{car:>8.2f} {ci:>7d} {-cs[1]:>10.4f} {abs(-cs[1] - car):>7.4f} "
            f"{out[0]:>8.3f} {out[2]:>8.3f} {out[1]:>8.3f} {out[3]:>8.3f}")
    e = np.array(errs)
    say(f"max sub-bin carrier error {e[:, 0].max():.4f} bins; worst focus error: "
        f"imag+int {e[:, 1].max():.3f}, imag+sub {e[:, 3].max():.3f}, "
        f"real+int {e[:, 2].max():.3f}, real+sub {e[:, 4].max():.3f} mm")
    rng = np.random.default_rng(5)
    errs2 = []
    for _ in range(20):
        cy, cx = rng.uniform(-50, 50), rng.uniform(55, 100)
        h = W.record_hologram(sensor, cx, carrier_y=cy,
                              ref_phase=rng.uniform(0, 2 * np.pi))
        ky, kx = W.estimate_carrier_2d(h, dc_radius=30.0 / N)
        errs2.append(np.hypot(ky + cy, kx + cx))
    say(f"20 random 2-D carriers (|ky|<50, kx 55..100): median error "
        f"{np.median(errs2):.4f} bins, max {np.max(errs2):.4f}")


# ---------------------------------------------------------------- S4
def lowpass_object(obj, cutoff):
    F = np.fft.fftshift(np.fft.fft2(obj))
    V, U = np.mgrid[0:N, 0:N] - N // 2
    F[np.hypot(V, U) > cutoff] = 0
    return np.real(np.fft.ifft2(np.fft.ifftshift(F)))


def s4_bias():
    say("\n== S4. RESIDUAL BIAS = SIDEBAND CROSS-TALK ==")
    obj = W.make_object()
    sensor = W.angular_spectrum(obj.astype(complex), 5e-3)
    say("Controls (truth -5.00 mm): true sensor field and an ideal")
    say("low-passed sideband give exactly 0; the hologram adds twin and")
    say("|U|^2 terms whose broadband tails leak into the sideband window.")
    r = W.find_focus(sensor, COARSE, tol=1e-10)
    say(f"  true sensor field: err {abs(r['z'] + 5e-3) * 1e6:.3f} um")
    rows = []
    for name, o in (("sharp object", obj),
                    ("band-limited (r<40 bins)", lowpass_object(obj, 40)),
                    ("band-limited (r<24 bins)", lowpass_object(obj, 24))):
        line = f"  {name:>26}:"
        for car, win in ((80, 56), (80, 48), (80, 40), (80, 32), (80, 24),
                         (100, 27)):
            f0, _ = holo_field(o, 5e-3, car=car, win=win, ref_phase=0.9)
            rr = W.find_focus(f0, COARSE, tol=1e-10)
            err = (rr["z"] + 5e-3) * 1e6
            rows.append((name, car, win, err))
            line += f" c{car}/w{win} {err:+6.1f}"
        say(line + "  (signed err, um)")
    fig, ax = fig_ax()
    for name in ("sharp object", "band-limited (r<40 bins)", "band-limited (r<24 bins)"):
        sel = [(w, e) for (n, c, w, e) in rows if n == name and c == 80]
        ax.plot([s[0] for s in sel], [s[1] for s in sel], "o-", label=name)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("sideband half-width [bins] (carrier 80)")
    ax.set_ylabel("focus bias [um] (truth 5 mm)")
    ax.set_title("Bias tracks how much twin/DC leaks into the window")
    ax.legend()
    save(fig, "s4_bias_vs_window.png")


# ---------------------------------------------------------------- S5
def s5_refine():
    say("\n== S5. OFF-GRID FOCUS: COARSE 0.1 mm SWEEP + GOLDEN SECTION ==")
    obj = W.make_object()
    rng = np.random.default_rng(23)
    say(f"{'truth':>8} {'noise':>6} {'coarse':>8} {'refined':>9} {'err um':>7} "
        f"{'imag 0.1mm':>10} {'imag err':>8} {'evals':>5}")
    for z in (3.3e-3, 5.03e-3, 6.137e-3, 7.77e-3, 9.1e-3):
        for noise in (None, ("gaussian", 0.2)):
            f0, _ = holo_field(obj, z, ref_phase=rng.uniform(0, 2 * np.pi),
                               noise=noise, rng=rng)
            calls = [0]

            def counted(f):
                calls[0] += 1
                return W.realness(f)
            r = W.find_focus(f0, COARSE, metric=counted, tol=1e-8)
            ri = W.find_focus(f0, COARSE, metric="imag_energy", refine=False)
            say(f"{-z * 1e3:>+8.3f} {('0.20' if noise else '-'):>6} "
                f"{r['coarse_z'] * 1e3:>+8.2f} {r['z'] * 1e3:>+9.4f} "
                f"{abs(r['z'] + z) * 1e6:>7.1f} {ri['z'] * 1e3:>+10.2f} "
                f"{abs(ri['z'] + z) * 1e3:>8.3f} {calls[0]:>5}")


# ---------------------------------------------------------------- S6
def s6_noise(quick):
    say("\n== S6. NOISE WITH A RANDOM REFERENCE PISTON ==")
    say("Each trial draws a new piston and a new noise field; carrier is")
    say("estimated. Entry: median |error| mm / trials within 0.1 mm.")
    obj = W.make_object()
    seeds = 6 if quick else 12
    cases = [("gaussian", g) for g in (0.1, 0.2, 0.35, 0.5, 0.8, 1.0, 1.5, 2.0)]
    cases += [("poisson", p) for p in (1e4, 1e3, 3e2, 1e2, 3e1)]
    say(f"{'noise':>16} {'imag_energy':>14} {'realness':>14} {'ac_realness':>14}")
    table = []
    for kind, lv in cases:
        res = {m: [] for m in ("imag_energy", "realness", "ac_realness")}
        for seed in range(seeds):
            rng = np.random.default_rng(1000 + seed)
            f0, _ = holo_field(obj, 5e-3, ref_phase=rng.uniform(0, 2 * np.pi),
                               noise=(kind, lv), rng=rng)
            for m in res:
                res[m].append(abs(W.find_focus(f0, COARSE, metric=m)["z"] + 5e-3) * 1e3)
        line = f"{kind + ' ' + format(lv, 'g'):>16}"
        for m, v in res.items():
            v = np.array(v)
            line += f" {np.median(v):>8.3f}/{(v < 0.1).sum():>2d}/{seeds}"
            table.append((kind, lv, m, np.median(v), (v < 0.1).mean()))
        say(line)
    fig, ax = fig_ax()
    for m, st in (("imag_energy", "o-"), ("realness", "s-"), ("ac_realness", "^-")):
        xs = [lv for (k, lv, mm, _, _) in table if k == "gaussian" and mm == m]
        ys = [ok for (k, lv, mm, _, ok) in table if k == "gaussian" and mm == m]
        ax.plot(xs, ys, st, label=m)
    ax.set_xlabel("gaussian noise sigma / mean hologram")
    ax.set_ylabel("fraction of trials within 0.1 mm")
    ax.set_title("Random piston: imag_energy fails even at low noise")
    ax.legend()
    save(fig, "s6_noise_random_piston.png")


# ---------------------------------------------------------------- S7
def s7_phase_objects():
    say("\n== S7. PHASE AND MIXED OBJECTS (truth -5.00, piston 0.7) ==")
    say("Entry: |error| mm / contrast (median/min). ok = error < 0.1 mm.")
    amp = W.make_object()
    y, x = np.mgrid[0:N, 0:N]
    disk = (np.hypot(x - 128, y - 128) < 100).astype(float)
    small = (np.hypot(x - 60, y - 200) < 30).astype(float)
    rng = np.random.default_rng(3)
    smooth = np.zeros((N, N))
    for _ in range(12):
        cx, cy = rng.integers(40, 216, 2)
        smooth += np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * 12 ** 2))
    smooth /= smooth.max()
    objs = [
        ("amplitude, dark bg", amp),
        ("amplitude, bright bg", 0.3 + 0.7 * amp),
        ("binary phase pi/8", np.exp(1j * np.pi / 8 * amp)),
        ("binary phase pi/2", np.exp(1j * np.pi / 2 * amp)),
        ("binary phase 3pi/4", np.exp(1j * 3 * np.pi / 4 * amp)),
        ("binary phase pi", np.exp(1j * np.pi * amp)),
        ("phase disk pi/2", np.exp(1j * np.pi / 2 * disk)),
        ("1 material amp+phase", (1 - 0.4 * amp) * np.exp(1j * 0.8 * amp)),
        ("3-level phase", np.exp(1j * (np.pi / 2 * amp + np.pi / 3 * small))),
        ("2 materials, own supp.", (0.3 + 0.7 * amp) *
         np.exp(1j * 1.2 * np.roll(amp, (30, -40), axis=(0, 1)))),
        ("smooth phase 0.5 rad", np.exp(1j * 0.5 * smooth)),
        ("smooth phase 2 rad", np.exp(1j * 2 * smooth)),
    ]
    mets = ["imag_energy", "realness", "ac_realness", "amplitude_variance",
            "tamura", "phase_variance"]
    say(f"{'object':>24}" + "".join(f"{m[:13]:>14}" for m in mets))
    score = {m: 0 for m in mets}
    for name, o in objs:
        f0, _ = holo_field(o, 5e-3, ref_phase=0.7)
        line = f"{name:>24}"
        for m in mets:
            r = W.find_focus(f0, COARSE, metric=m)
            e = abs(r["z"] + 5e-3) * 1e3
            score[m] += e < 0.1
            line += f"  {e:>6.3f}/{min(r['contrast'], 99):>4.1f}"
        say(line)
    say("objects found within 0.1 mm: " +
        ", ".join(f"{m} {score[m]}/{len(objs)}" for m in mets))


# ---------------------------------------------------------------- S8
def s8_two_plane():
    say("\n== S8. TWO-PLANE SAMPLE: PER-TILE FOCUS MAP ==")
    amp = W.make_object()
    y, x = np.mgrid[0:N, 0:N]
    dots = np.zeros((N, N))
    for dx0, dy0 in ((-70, -70), (70, -70), (0, 70)):
        dots[np.hypot(x - (128 + dx0), y - (128 + dy0)) < 10] = 1
    p1 = np.clip(amp - dots, 0, 1)
    sensor = (W.angular_spectrum(p1.astype(complex), 5e-3)
              + W.angular_spectrum(dots.astype(complex), 7.5e-3))
    h = W.record_hologram(sensor, ref_phase=1.3)
    f0 = W.reconstruct_field(h, carrier=W.estimate_carrier_2d(h, dc_radius=30.0 / N))
    sweep = np.linspace(-10e-3, 0, 101)
    g = W.find_focus(f0, sweep)
    say(f"global realness focus: {g['z'] * 1e3:+.3f} mm (one number for two planes)")
    tile = 32
    zmap, _ = W.local_focus_map(f0, sweep, tile=tile)
    nt = N // tile
    truth = np.full((nt, nt), np.nan)
    for i in range(nt):
        for j in range(nt):
            a = p1[i * tile:(i + 1) * tile, j * tile:(j + 1) * tile].sum()
            b = dots[i * tile:(i + 1) * tile, j * tile:(j + 1) * tile].sum()
            if a + b > 0:
                truth[i, j] = -5.0 if a > b else -7.5
    both = ~np.isnan(zmap) & ~np.isnan(truth)
    err = np.abs(zmap[both] * 1e3 - truth[both])
    say(f"tiles with signal: {both.sum()}; |error| median {np.median(err):.3f} mm,"
        f" within 0.25 mm: {(err < 0.25).sum()}/{both.sum()}, "
        f"assigned to the right plane (nearer truth): "
        f"{(err < 1.25).sum()}/{both.sum()}")
    fig, axes = fig_ax(1, 3, size=(13, 4))
    axes[0].imshow(p1 + 0.5 * dots, origin="lower", cmap="gray")
    axes[0].set_title("sample: ring+cross @5 mm, dots @7.5 mm")
    im = axes[1].imshow(-zmap * 1e3, origin="lower", cmap="viridis", vmin=4.5, vmax=8)
    axes[1].set_title("per-tile realness focus [mm]")
    fig.colorbar(im, ax=axes[1], fraction=0.046)
    im2 = axes[2].imshow(-truth, origin="lower", cmap="viridis", vmin=4.5, vmax=8)
    axes[2].set_title("truth per tile [mm]")
    fig.colorbar(im2, ax=axes[2], fraction=0.046)
    for a in axes:
        a.set_xticks([])
        a.set_yticks([])
    save(fig, "s8_two_plane_depth_map.png")


# ---------------------------------------------------------------- S9
def pair_spectrum_curve(field0, sweep, dx=W.DX, wl=W.WL):
    """realness(z) from the pair spectrum P(k)=F(k)F(-k): no FFT per z."""
    return W.pair_spectrum_focus_curve(field0, sweep, dx=dx, wl=wl)


def s9_identity():
    say("\n== S9. SPECTRAL IDENTITY ==")
    say("sum U_z^2 = sum_k F(k) F(-k) exp(2 i z kz(k)), and sum |U_z|^2 is")
    say("constant (Parseval). So realness(z) depends ONLY on the pair")
    say("spectrum P(k) = F(k)F(-k); paraxially kz ~ k - pi*lambda*|f|^2, so")
    say("the focus curve is |1-D Fourier transform of P over |f|^2| at")
    say("'frequency' lambda*z: autofocus is frequency estimation.")
    obj = W.make_object()
    f0, _ = holo_field(obj, 5e-3, ref_phase=0.4)
    t = time.time()
    direct = np.array([W.realness(W.angular_spectrum(f0, z)) for z in COARSE])
    t1 = time.time() - t
    t = time.time()
    spec = pair_spectrum_curve(f0, COARSE)
    t2 = time.time() - t
    say(f"max |direct - pair-spectrum| over 101 planes: "
        f"{np.max(np.abs(direct - spec)):.2e}  (direct {t1:.2f} s, "
        f"pair-spectrum {t2:.2f} s)")
    # Consequence: noise enters through P(k) only; tested weightings of P
    # (|P|^g, |f|^a, high-pass) did not beat the flat sum (see RESULTS).


# ---------------------------------------------------------------- S10
def s10_fly():
    say("\n== S10. REAL FLY HOLOGRAM (pyDHM fly-20x.bmp, 532 nm, 2.4 um) ==")
    path = os.path.join(HERE, "real_data", "fly-20x.bmp")
    if not os.path.exists(path):
        say("  fly-20x.bmp missing; skipped")
        return
    import matplotlib.image as mpimg
    h = mpimg.imread(path).astype(float)
    if h.ndim == 3:
        h = h[..., :3].mean(-1)
    car = W.estimate_carrier_2d(h, dc_radius=0.1, half="any")
    win = 200
    f0 = W.reconstruct_field(h, carrier=car, win=win, crop=True)
    pitch = W.cropped_pitch(h.shape, win, 2.4e-6)
    say(f"carrier (ky, kx) = ({car[0]:+.3f}, {car[1]:+.3f}) bins; cropped "
        f"field {f0.shape}, pitch {pitch[0] * 1e6:.2f} um")
    sweep = np.arange(-10, 10.001, 0.25) * 1e-3
    curves = {}
    for m in ("imag_energy", "realness", "ac_realness", "amplitude_variance",
              "tamura"):
        r = W.find_focus(f0, sweep, metric=m, dx=pitch, wl=0.532e-6)
        curves[m] = (r["curve"], r["z"])
        say(f"  {m:>20}: focus {r['z'] * 1e3:+.3f} mm, contrast {r['contrast']:.3f}")
    fig, axes = fig_ax(1, 2, size=(13, 4))
    for m, (c, z) in curves.items():
        if m == "imag_energy":
            continue
        cn = (c - c.min()) / (c.max() - c.min())
        axes[0].plot(sweep * 1e3, cn, label=f"{m} ({z * 1e3:+.2f})")
    axes[0].set_xlabel("z [mm]")
    axes[0].set_ylabel("normalised metric (min = focus)")
    axes[0].set_title("Fly: ac_realness agrees with amp-variance and Tamura")
    axes[0].legend(fontsize=8)
    c, z = curves["imag_energy"]
    axes[1].plot(sweep * 1e3, c, "k.-", ms=3)
    axes[1].set_title("imag_energy on the fly: lambda/2 aliasing, no basin")
    axes[1].set_xlabel("z [mm]")
    save(fig, "s10_fly_metrics.png")
    zf = curves["ac_realness"][1]
    best = W.angular_spectrum(f0, zf, dx=pitch, wl=0.532e-6)
    fig, axes = fig_ax(1, 2, size=(10, 5))
    axes[0].imshow(np.abs(f0), cmap="gray")
    axes[0].set_title("|U| at hologram plane (z = 0)")
    axes[1].imshow(np.abs(best), cmap="gray")
    axes[1].set_title(f"|U| at ac_realness focus {zf * 1e3:+.2f} mm")
    for a in axes:
        a.set_xticks([])
        a.set_yticks([])
    save(fig, "s10_fly_refocus.png")


def main(argv=None):
    global OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--outdir", default=os.path.join(HERE, "results"))
    ap.add_argument("--only", default="")
    a = ap.parse_args(argv)
    OUT = a.outdir
    os.makedirs(OUT, exist_ok=True)
    say(f"wave-hologram {W.__version__} studies (numpy {np.__version__})")
    t0 = time.time()
    steps = {
        "1": s1_piston, "2": lambda: s2_comb(a.quick), "3": s3_carrier,
        "4": s4_bias, "5": s5_refine, "6": lambda: s6_noise(a.quick),
        "7": s7_phase_objects, "8": s8_two_plane, "9": s9_identity,
        "10": s10_fly,
    }
    only = [s for s in a.only.split(",") if s] or list(steps)
    for k in only:
        steps[k]()
    say(f"\ntotal {time.time() - t0:.0f} s")
    name = "studies_output.txt" if not a.only else f"studies_output_{a.only.replace(',', '_')}.txt"
    with open(os.path.join(OUT, name), "w") as fh:
        fh.write("\n".join(LOG) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
