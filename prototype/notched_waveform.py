"""
Spectrally notched, constant-envelope radar waveform design.

Our independent reconstruction of the methodology presented by Spencer Markowitz,
"Optimal Waveform Design in the Presence of Interference" (slide footer: "Spectral
notching for chirp waveforms"), Radio Club of America Technical Symposium,
New York, November 21, 2026. Nothing here is the speaker's code; parameters
he did not state are our assumptions and are marked ASSUMPTION.

Four waveform designs are compared:
  1. lfm     - plain linear FM chirp (no notch)
  2. jump    - linear FM that jumps over the protected band (phase-continuous)
  3. ap_flat - alternating projections toward a flat spectrum with a notch
  4. convex  - convex power-spectrum design for a chosen autocorrelation,
               then alternating projections to realize it with constant envelope

Step 1 (convex): choose g_k >= 0 on an M-point grid (M >= 2N) to minimize the
peak range sidelobe s subject to
    r = IDFT(g)                       (linear in g)
    r_0 = 1                           (normalization)
    |r_m| <= s     for K <= |m| < N   (sidelobes outside the main-lobe zone)
    r_m = 0        for N <= |m| <= M/2 (autocorrelation of a length-N sequence)
    g_k <= notch_level   in protected and out-of-band bins
    g_k <= g_cap         (keeps the spectrum flat enough to realize)
Step 2 (alternating projections): FFT -> impose |X| = sqrt(g) keep phase;
IFFT -> impose |x| = 1 on the N-sample support, zero elsewhere; repeat; stop
on the envelope step so the transmitted waveform is exactly constant envelope.
"""

import argparse
import json
import time

import numpy as np
import cvxpy as cp


# ---------------------------------------------------------------- parameters
def default_params():
    return dict(
        B=150e6,          # swept bandwidth, Hz (from talk: ~1 m resolution)
        fs=300e6,         # complex sample rate, Hz (ASSUMPTION: 2x oversampling)
        N=1024,           # samples per pulse (ASSUMPTION; T = N/fs = 3.4 us)
        notch_center=30e6,  # Hz offset from band center (ASSUMPTION)
        notch_width=25e6,   # Hz (talk: narrowband user 20-30 MHz wide)
        notch_db=-30.0,     # target notch depth (talk: ~30 dB)
        oob_db=-30.0,       # out-of-band mask (ASSUMPTION)
        K=4,                # main-lobe exclusion half-width in samples (spoiling knob)
        g_cap=3.0,          # max PSD relative to flat in-band level (ASSUMPTION)
        ap_iters=2000,      # alternating-projection iterations (ASSUMPTION)
    )


def freq_grid(M, fs):
    return np.fft.fftfreq(M, d=1 / fs)


def masks(p, M):
    f = freq_grid(M, p["fs"])
    inband = np.abs(f) <= p["B"] / 2
    notch = np.abs(f - p["notch_center"]) <= p["notch_width"] / 2
    allowed = inband & ~notch
    return f, inband, notch, allowed


# ---------------------------------------------------------------- waveforms
def lfm(p):
    n = np.arange(p["N"])
    t = n / p["fs"]
    T = p["N"] / p["fs"]
    k = p["B"] / T
    phase = 2 * np.pi * (-p["B"] / 2 * t + 0.5 * k * t ** 2)
    return np.exp(1j * phase)


def jump_lfm(p):
    """Linear FM over the allowed frequencies only: instantaneous frequency jumps
    across the protected band (phase stays continuous, frequency does not)."""
    N, fs, B = p["N"], p["fs"], p["B"]
    lo, hi = p["notch_center"] - p["notch_width"] / 2, p["notch_center"] + p["notch_width"] / 2
    span = B - (hi - lo)
    u = np.arange(N) / N * span - B / 2          # sweep over a compressed axis
    finst = np.where(u < lo, u, u + (hi - lo))   # expand around the notch
    phase = 2 * np.pi * np.cumsum(finst) / fs
    return np.exp(1j * phase)


def alternating_projections(g_target, p, x0, iters):
    """Realize a constant-envelope length-N sequence whose |DFT_M|^2 approximates g_target."""
    N = p["N"]
    M = len(g_target)
    mag = np.sqrt(np.maximum(g_target, 0) * N)   # Parseval scaling for |x| = 1
    x = np.zeros(M, complex)
    x[:N] = x0
    for _ in range(iters):
        X = np.fft.fft(x)
        X = mag * np.exp(1j * np.angle(X))
        y = np.fft.ifft(X)
        x = np.zeros(M, complex)
        x[:N] = np.exp(1j * np.angle(y[:N]))     # stop on the envelope step
    return x[:N]


def flat_notched_target(p, M):
    f, inband, notch, allowed = masks(p, M)
    g = np.full(M, 10 ** (p["oob_db"] / 10))
    g[allowed] = 1.0
    g[notch & inband] = 10 ** (p["notch_db"] / 10)
    return g * M / g.sum()


def convex_target(p, M, solver="CLARABEL", verbose=False):
    N, K = p["N"], p["K"]
    f, inband, notch, allowed = masks(p, M)
    n_in = allowed.sum()
    flat = M / n_in                              # level if energy spread evenly on allowed bins
    notch_lvl = 10 ** (p["notch_db"] / 10) * flat
    oob_lvl = 10 ** (p["oob_db"] / 10) * flat

    lags = np.arange(0, M // 2 + 1)
    E = np.exp(2j * np.pi * np.outer(lags, np.arange(M)) / M) / M   # r_m = (1/M) sum g_k e^{+j2pi km/M}
    Er, Ei = E.real, E.imag

    g = cp.Variable(M, nonneg=True)
    s = cp.Variable()
    rr, ri = Er @ g, Ei @ g
    side = np.arange(K, N)
    zero = np.arange(N, M // 2 + 1)
    cons = [
        rr[0] == 1,
        cp.norm(cp.vstack([rr[side], ri[side]]), 2, axis=0) <= s,
        rr[zero] == 0, ri[zero] == 0,
        g[notch & inband] <= notch_lvl,
        g[~inband] <= oob_lvl,
        g <= p["g_cap"] * flat,
    ]
    prob = cp.Problem(cp.Minimize(s), cons)
    t0 = time.time()
    prob.solve(solver=solver, verbose=verbose)
    return np.maximum(g.value, 0), float(s.value), prob.status, time.time() - t0


# ---------------------------------------------------------------- metrics
def metrics(x, p, M):
    N = p["N"]
    X = np.fft.fft(x, M)
    P = np.abs(X) ** 2
    f, inband, notch, allowed = masks(p, M)
    ref = P[allowed].mean()
    notch_depth = 10 * np.log10(P[notch & inband].mean() / ref)
    notch_worst = 10 * np.log10(P[notch & inband].max() / ref)
    r = np.fft.ifft(P)                           # aperiodic autocorrelation (M >= 2N)
    r = np.abs(r[:N]) / np.abs(r[0])
    # main-lobe: first local minimum, then PSL beyond it
    k = 1
    while k < N - 1 and r[k + 1] < r[k]:
        k += 1
    psl = 20 * np.log10(r[k:].max())
    isl = 10 * np.log10(2 * np.sum(r[k:] ** 2))
    above = np.where(r < 10 ** (-3 / 20))[0]
    hw3 = above[0] if len(above) else np.nan    # -3 dB half-width, samples
    res_m = 3e8 / 2 * (2 * hw3) / p["fs"] if np.isfinite(hw3) else np.nan
    env = np.abs(x)
    return dict(
        notch_mean_db=round(notch_depth, 1),
        notch_worst_db=round(notch_worst, 1),
        psl_db=round(psl, 1),
        isl_db=round(isl, 1),
        mainlobe_null_samples=int(k),
        approx_3db_width_m=round(res_m, 2),
        envelope_ripple=float(env.max() - env.min()),
    )


def interference_gain_db(x, x_ref, p, M):
    """OUR metric (not the speaker's): how much less power a narrowband interferer
    filling the protected band couples through the matched filter, relative to the
    reference waveform, at equal target peak. Positive = more interference tolerated."""
    f, inband, notch, allowed = masks(p, M)

    def coupling(w):
        P = np.abs(np.fft.fft(w, M)) ** 2
        peak = np.sum(np.abs(w) ** 2) ** 2      # matched-filter peak power
        return P[notch & inband].sum() / peak
    return round(10 * np.log10(coupling(x_ref) / coupling(x)), 1)


# ---------------------------------------------------------------- main
def run(p, plot=True, outdir="results"):
    M = 2 * p["N"]
    x_lfm = lfm(p)
    x_jump = jump_lfm(p)
    x_apf = alternating_projections(flat_notched_target(p, M), p, x_lfm, p["ap_iters"])
    g_cvx, s_opt, status, t_solve = convex_target(p, M)
    t0 = time.time()
    x_cvx = alternating_projections(g_cvx, p, x_lfm, p["ap_iters"])
    t_ap = time.time() - t0

    designs = dict(lfm=x_lfm, jump=x_jump, ap_flat=x_apf, convex=x_cvx)
    out = {name: metrics(x, p, M) for name, x in designs.items()}
    for name, x in designs.items():
        out[name]["interference_gain_vs_lfm_db"] = interference_gain_db(x, x_lfm, p, M)
    out["_convex_solver"] = dict(status=status, optimal_psl_db=round(20 * np.log10(s_opt), 1),
                                 solve_s=round(t_solve, 1), ap_s=round(t_ap, 2))
    out["_params"] = p

    if plot:
        make_plots(designs, g_cvx, p, M, outdir)
    return out


def make_plots(designs, g_cvx, p, M, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    f, inband, notch, allowed = masks(p, M)
    order = np.argsort(f)
    fig, ax = plt.subplots(2, 1, figsize=(10, 8))
    for name, x in designs.items():
        P = np.abs(np.fft.fft(x, M)) ** 2
        P = P / P[allowed].mean()
        ax[0].plot(f[order] / 1e6, 10 * np.log10(P[order] + 1e-12), lw=0.8, label=name)
    G = g_cvx / g_cvx[allowed].mean()
    ax[0].plot(f[order] / 1e6, 10 * np.log10(G[order] + 1e-12), "k--", lw=1, label="convex target")
    ax[0].axvspan((p["notch_center"] - p["notch_width"] / 2) / 1e6,
                  (p["notch_center"] + p["notch_width"] / 2) / 1e6, color="0.85")
    ax[0].set(xlabel="Frequency offset (MHz)", ylabel="PSD re in-band mean (dB)",
              ylim=(-60, 10), title="Transmit spectrum (shaded: protected band)")
    ax[0].legend(fontsize=8)

    N = p["N"]
    lag_m = np.arange(N) * 3e8 / 2 / p["fs"]
    for name, x in designs.items():
        P = np.abs(np.fft.fft(x, M)) ** 2
        r = np.abs(np.fft.ifft(P)[:N])
        ax[1].plot(lag_m, 20 * np.log10(r / r[0] + 1e-12), lw=0.8, label=name)
    ax[1].set(xlabel="Range offset (m)", ylabel="Matched-filter output (dB)", ylim=(-70, 2),
              xlim=(0, 60), title="Range response (autocorrelation)")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{outdir}/spectrum_and_range.png", dpi=130)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--K", type=int, default=None, help="main-lobe exclusion (spoiling) in samples")
    ap.add_argument("--iters", type=int, default=None)
    ap.add_argument("--N", type=int, default=None, help="samples per pulse")
    ap.add_argument("--outdir", default="results")
    a = ap.parse_args()
    p = default_params()
    if a.N is not None:
        p["N"] = a.N
    if a.K is not None:
        p["K"] = a.K
    if a.iters is not None:
        p["ap_iters"] = a.iters
    res = run(p, outdir=a.outdir)
    print(json.dumps(res, indent=2))
