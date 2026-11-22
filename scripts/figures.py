"""Draw the README figures from docs/data/*.csv and *.pbm.

    cmake -S . -B build && cmake --build build -j
    ./build/thermodynamics && ./build/binder && ./build/scaling && ./build/dynamics && ./build/snapshots
    python3 scripts/figures.py
"""

import csv
import json
import math
import pathlib

import matplotlib
import matplotlib.ticker
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA, DOCS = ROOT / "docs" / "data", ROOT / "docs"
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
BLUE, ORANGE = "#2a78d6", "#eb6834"
RAMP = {8: "#86b6ef", 16: "#5598e7", 32: "#2a78d6", 64: "#1c5cab"}
TC = 2 / math.log(1 + math.sqrt(2))
U_STAR = 0.61069  # Kamieniarz & Bloete 1993, square lattice, periodic boundaries

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlecolor": INK,
})


def rows(name):
    return list(csv.DictReader(open(DATA / name)))


def col(rs, key, cast=float):
    return np.array([cast(r[key]) for r in rs])


def style(ax, grid_axis="both"):
    ax.grid(axis=grid_axis, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def elliptic_k(k):
    a, b = 1.0, math.sqrt(1 - k * k)
    for _ in range(60):
        a, b = 0.5 * (a + b), math.sqrt(a * b)
    return math.pi / (2 * a)


def onsager_e(T):
    b = 1 / T
    t2 = math.tanh(2 * b) ** 2
    kappa = 2 * math.sinh(2 * b) / math.cosh(2 * b) ** 2
    if abs(2 * t2 - 1) < 1e-15:
        return -1 / math.tanh(2 * b)
    return -(1 / math.tanh(2 * b)) * (1 + (2 / math.pi) * (2 * t2 - 1) * elliptic_k(min(kappa, 1 - 1e-16)))


def onsager_m(T):
    x = 1 - math.sinh(2 / T) ** -4
    return x ** 0.125 if x > 0 else 0.0


def snapshots():
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.4), gridspec_kw={"wspace": 0.06})
    for ax, (T, label) in zip(axes, [(2.0, "T = 2.0, ordered"), (TC, "T = T_c ≈ 2.269, critical"),
                                     (2.6, "T = 2.6, disordered")]):
        raw = (DATA / f"snapshot_T{T:.3f}.pbm").read_bytes()
        header_end = raw.index(b"\n", raw.index(b"\n") + 1) + 1
        L = int(raw[3:header_end].split()[0])
        bits = np.unpackbits(np.frombuffer(raw[header_end:], dtype=np.uint8)).reshape(L, L)
        ax.imshow(bits, cmap="Greys_r", interpolation="nearest")
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(label, loc="left")
    fig.text(0.125, 0.04, "256 × 256 periodic lattice after 3,000 Wolff sweeps from a random start; "
             "white = spin up. At T_c there are clusters of every size.", fontsize=9, color=MUTED)
    fig.savefig(DOCS / "snapshots.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def thermodynamics():
    rs = rows("thermodynamics.csv")
    Ts = np.linspace(1.5, 3.5, 801)
    exact_e = np.array([onsager_e(t) for t in Ts])
    exact_m = np.array([onsager_m(t) for t in Ts])
    h = 1e-5
    exact_c = np.array([(onsager_e(t + h) - onsager_e(t - h)) / (2 * h) if abs(t - TC) > 2e-3 else np.nan for t in Ts])
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), gridspec_kw={"hspace": 0.32, "wspace": 0.22})
    panels = [("e", "e_err", "energy per spin", exact_e), ("abs_m", "abs_m_err", "magnetization |m|", exact_m),
              ("c", "c_err", "specific heat per spin", exact_c), ("chi", "chi_err", "susceptibility per spin", None)]
    for ax, (key, err, label, exact) in zip(axes.flat, panels):
        if exact is not None:
            ax.plot(Ts, exact, color=INK, lw=1.4, label="Onsager, infinite lattice", zorder=1)
        for L, color in [(16, RAMP[16]), (64, RAMP[64])]:
            sel = [r for r in rs if int(r["L"]) == L]
            sel.sort(key=lambda r: float(r["T"]))
            ax.errorbar(col(sel, "T"), col(sel, key), yerr=col(sel, err), fmt="o", ms=4, color=color,
                        lw=1, capsize=0, label=f"Monte Carlo, L = {L}", zorder=2)
        ax.axvline(TC, color=MUTED, lw=0.8, ls=(0, (4, 3)))
        ax.set_xlabel("temperature T")
        ax.set_title(label, loc="left")
        style(ax)
    axes[0, 0].legend(frameon=False, fontsize=9, loc="upper left")
    axes[1, 1].set_yscale("log")
    fig.text(0.01, 0.01, "Wolff cluster updates; dashed line at T_c = 2/ln(1+√2). The exact specific heat is the "
             "derivative of Onsager's energy and diverges logarithmically at T_c.", fontsize=9, color=MUTED)
    fig.savefig(DOCS / "thermodynamics.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def binder():
    curves, points = rows("binder_curves.csv"), rows("binder_points.csv")
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"wspace": 0.22})
    curve = {}
    for L in (8, 16, 32, 64):
        sel = [r for r in curves if int(r["L"]) == L]
        T, U = col(sel, "T"), col(sel, "binder")
        curve[L] = (T, U)
        for a in (ax, ax2):
            a.plot(T, U, color=RAMP[L], lw=2, label=f"L = {L}")
        p = [r for r in points if int(r["L"]) == L]
        for a in (ax, ax2):
            a.errorbar(col(p, "T"), col(p, "binder"), yerr=col(p, "binder_err"), fmt="o", ms=4,
                       mfc=SURFACE, mec=RAMP[L], ecolor=RAMP[L], capsize=0, zorder=3)
    crossings = []
    for a_, b_ in [(8, 16), (16, 32), (32, 64)]:
        (T, Ua), (_, Ub) = curve[a_], curve[b_]
        d = Ua - Ub
        i = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
        j = i[np.argmin(np.abs(T[i] - TC))]
        t = T[j] - d[j] * (T[j + 1] - T[j]) / (d[j + 1] - d[j])
        u = np.interp(t, T, Ua)
        crossings.append({"pair": [a_, b_], "T": t, "U": u})
        ax2.plot(t, u, "x", color=INK, ms=8, mew=1.5)
    for a in (ax, ax2):
        a.axvline(TC, color=MUTED, lw=0.8, ls=(0, (4, 3)))
        a.axhline(U_STAR, color=MUTED, lw=0.8, ls=(0, (1, 2)))
        a.set_xlabel("temperature T")
        style(a)
    ax.set_ylabel("Binder cumulant U₄ = 1 − ⟨m⁴⟩ / 3⟨m²⟩²")
    ax.legend(frameon=False, fontsize=9, loc="lower left")
    ax.set_title("The curves for every size cross at T_c", loc="left")
    ax2.set_xlim(2.255, 2.285)
    ax2.set_ylim(0.585, 0.635)
    ax2.set_title("Zoom: crossings (×) next to T_c and U* = 0.61069", loc="left")
    table = "crossing of L and 2L:\n" + "\n".join(
        f"{c['pair'][0]:>2} / {c['pair'][1]:<2}  T = {c['T']:.4f}  U = {c['U']:.4f}" for c in crossings)
    table += f"\nexact T_c = {TC:.4f}"
    ax2.text(0.03, 0.04, table, transform=ax2.transAxes, fontsize=8.5, color=INK2, family="monospace", va="bottom",
             bbox={"boxstyle": "round,pad=0.4", "fc": SURFACE, "ec": GRID})
    fig.text(0.01, -0.03, "Lines: one Wolff run per size at T_c (100,000–200,000 samples), reweighted to nearby T. "
             "Circles: independent runs (40,000 samples) with jackknife errors.", fontsize=9, color=MUTED)
    fig.savefig(DOCS / "binder.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    return crossings


def scaling():
    rs = rows("scaling.csv")
    L = col(rs, "L")
    fits = {}
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"wspace": 0.25})
    for ax, (key, err, label, exact, ex_label) in zip(axes, [
            ("abs_m", "abs_m_err", "⟨|m|⟩ at T_c", -1 / 8, "−β/ν = −1/8"),
            ("chi_m2", "chi_m2_err", "N⟨m²⟩ / T at T_c", 7 / 4, "γ/ν = 7/4")]):
        y, dy = col(rs, key), col(rs, err)
        w = y / dy  # weights 1/sigma on the log scale
        slope, icpt = np.polyfit(np.log(L), np.log(y), 1, w=w)
        cov = np.polyfit(np.log(L), np.log(y), 1, w=w, cov=True)[1]
        fits[key] = {"slope": slope, "slope_err": math.sqrt(cov[0, 0]), "exact": exact}
        xx = np.array([6, 340])
        ax.plot(xx, np.exp(icpt) * xx ** exact, color=MUTED, lw=1, ls=(0, (4, 3)),
                label=f"exact slope {ex_label}")
        ax.errorbar(L, y, yerr=dy, fmt="o", ms=6, color=BLUE, capsize=0, zorder=3,
                    label=f"fitted slope {slope:+.4f} ± {math.sqrt(cov[0, 0]):.4f}")
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xticks(L, [str(int(v)) for v in L])
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel("lattice size L")
        ax.set_title(label, loc="left")
        ax.legend(frameon=False, fontsize=9)
        style(ax)
    fig.text(0.01, -0.03, "Wolff runs at the exact T_c, L = 8 to 256 (10,000–50,000 samples). "
             "Weighted least squares on the log-log points.", fontsize=9, color=MUTED)
    fig.savefig(DOCS / "scaling.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    return fits


def dynamics():
    rs = rows("dynamics.csv")
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    fits = {}
    for algo, color in [("metropolis", ORANGE), ("wolff", BLUE)]:
        for obs, marker in [("abs_m", "o"), ("e", "s")]:
            sel = [r for r in rs if r["algorithm"] == algo and r["observable"] == obs]
            L, tau, err = col(sel, "L"), col(sel, "tau"), col(sel, "tau_err")
            z, icpt = np.polyfit(np.log(L), np.log(tau), 1)
            fits[f"{algo}_{obs}"] = {"z": z, "L": [int(v) for v in L], "tau": list(tau)}
            label = f"{'Metropolis' if algo == 'metropolis' else 'Wolff'}, {'|m|' if obs == 'abs_m' else 'energy'}: z ≈ {z:.2f}"
            ax.errorbar(L, tau, yerr=err, fmt=marker, ms=6, color=color, mfc=color if obs == "abs_m" else SURFACE,
                        capsize=0, label=label)
            xx = np.array([L.min(), L.max()])
            ax.plot(xx, np.exp(icpt) * xx ** z, color=color, lw=1, alpha=0.7)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks([8, 16, 32, 64, 128], ["8", "16", "32", "64", "128"])
    ax.set_yticks([1, 2, 5, 10, 20, 50], ["1", "2", "5", "10", "20", "50"])
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.text(0.98, 0.42, "published: z = 2.1665(12) for single-spin flips\n"
            "(Nightingale & Blöte 1996); z ≈ 0.25 for Wolff", transform=ax.transAxes, ha="right", fontsize=8.5,
            color=MUTED)
    ax.set_xlabel("lattice size L")
    ax.set_ylabel("integrated autocorrelation time τ (sweeps)")
    ax.set_title("Critical slowing down at T_c: τ ∝ L^z", loc="left")
    ax.legend(frameon=False, fontsize=9)
    style(ax)
    fig.text(0.01, -0.04, "Metropolis: 10⁶–2·10⁶ sweeps per size. Wolff: 50,000–100,000 samples, converted to sweeps "
             "with the measured spins flipped per sample.", fontsize=9, color=MUTED)
    fig.savefig(DOCS / "dynamics.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    return fits


if __name__ == "__main__":
    snapshots()
    thermodynamics()
    crossings = binder()
    fits = scaling()
    dyn = dynamics()
    summary = {"binder_crossings": crossings, "scaling": fits, "dynamics": dyn, "Tc": TC, "U_star": U_STAR}
    (DOCS / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(summary, indent=2, default=float))
