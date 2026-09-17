"""
converge.py -- Mesh-convergence study of the mutual-inductance matrix M.

Runs run.py for several values of --seg while holding every other mesh
parameter fixed (--seg-far, --seg-max, --fine-radius, --grid), then fits

    M_ij(h) = M_ij(inf) + a_ij h^p_ij

coefficient by coefficient on the three finest meshes and extrapolates to
h -> 0 (Richardson extrapolation with observed order). The eigenvalues of
M(inf)^-1 are recomputed from the extrapolated matrix.

Usage (all other options are forwarded verbatim to run.py):

    python converge.py device.gds --seg 1.0,0.5,0.25 --Lsq 0.05 \
        --seg-far 2.5 --seg-max 15 --fine-radius 70 --grid 8 --out conv

If --seg-far, --seg-max or --fine-radius are omitted they are derived from
the SMALLEST seg with the defaults of build_model (5 seg, 30 seg,
200 seg / 1.5) and held fixed for all runs. --grid should be given
explicitly, otherwise the pFFT grid follows the median triangle size and
changes from one run to the next (a warning is printed).

With only two seg values the order cannot be estimated and --order (default
1) is used. Outputs <out>_convergence.md and <out>_convergence.npz.
"""

from __future__ import annotations

import argparse
import sys
import time
import numpy as np

try:
    from . import run as runmod
except ImportError:
    import run as runmod

PHI0 = 2.067833848e-15
FIXED = ("--seg-far", "--seg-max", "--fine-radius", "--grid")


# ----------------------------------------------------------------------
def _has_opt(argv, name):
    return any(a == name or a.startswith(name + "=") for a in argv)


def _capture_run(argv):
    """Call run.main(argv) and return whatever run.run returned."""
    box = {}
    orig = runmod.run

    def wrapped(*a, **kw):
        out = orig(*a, **kw)
        box["ret"] = out
        return out

    runmod.run = wrapped
    try:
        runmod.main(argv)
    finally:
        runmod.run = orig
    if "ret" not in box:
        raise RuntimeError("run.run was not called by run.main")
    return box["ret"]


def _extract(ret):
    """(model, sols, res) -> dict with M (2D, henry), names, energies."""
    model, sols, res = ret
    if res is None:
        raise RuntimeError("run returned no result (mesh-only mode ?)")
    M = np.atleast_2d(np.asarray(res["M"], dtype=float))
    if M.shape[0] == 1 and M.shape[1] > 1 and not isinstance(sols, list):
        M = M.T                          # legacy single-source layout
    out = dict(M=M)
    asm = model["assembly"]
    out["n_tri"] = int(asm.n_tri)
    out["h_med"] = float(np.median(asm.h_tri))
    out["h_min"] = float(asm.h_tri.min())
    sl = sols if isinstance(sols, (list, tuple)) else [sols]
    out["E"] = np.array([float(getattr(s, "energy", np.nan)) for s in sl])
    out["iters"] = np.array([int(s.info.get("iterations", -1)) for s in sl])
    for key in ("source_names", "source_labels", "sources"):
        if key in res:
            out["src"] = [str(x) for x in res[key]]
            break
    for key in ("surface_names", "surface_labels", "surfaces", "loop_names"):
        if key in res:
            out["surf"] = [str(x) for x in res[key]]
            break
    out.setdefault("src", [f"S{j}" for j in range(M.shape[1])])
    out.setdefault("surf", [f"L{i}" for i in range(M.shape[0])])
    return out


# ----------------------------------------------------------------------
def fit_order(h, m, p_min=0.05, p_max=4.0):
    """Observed order p from three (h, m) pairs, h decreasing.

    Solves  (m1 - m2)/(m2 - m3) = (h1^p - h2^p)/(h2^p - h3^p)  for p.
    Returns (p, ok). ok is False when the sequence is not monotone or no
    root exists in [p_min, p_max].
    """
    from scipy.optimize import brentq
    h1, h2, h3 = h
    m1, m2, m3 = m
    d12, d23 = m1 - m2, m2 - m3
    if d23 == 0 or d12 * d23 <= 0:
        return np.nan, False
    target = d12 / d23

    def f(p):
        return (h1 ** p - h2 ** p) / (h2 ** p - h3 ** p) - target

    try:
        if f(p_min) * f(p_max) > 0:
            return np.nan, False
        return brentq(f, p_min, p_max, xtol=1e-6), True
    except ValueError:
        return np.nan, False


def extrapolate(h, m, p):
    """M_inf from the two finest points and an order p."""
    h2, h3 = h[-2], h[-1]
    a = (m[-2] - m[-1]) / (h2 ** p - h3 ** p)
    return m[-1] - a * h3 ** p


def eig_inverse(M):
    """Eigenvalues of M^-1 in mA/Phi0, sorted by increasing modulus."""
    if M.shape[0] != M.shape[1]:
        return None
    try:
        lam = np.linalg.eigvals(np.linalg.inv(M)) * PHI0 * 1e3
    except np.linalg.LinAlgError:
        return None
    return lam[np.argsort(np.abs(lam))]


def _fmt_c(x):
    return f"{x.real:.5g}" if abs(x.imag) < 1e-9 * max(abs(x.real), 1e-300) \
        else f"{x.real:.5g}{x.imag:+.3g}j"


# ----------------------------------------------------------------------
def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    p = argparse.ArgumentParser(
        description="Mesh-convergence study of M over several --seg values. "
                    "Unknown options are forwarded to run.py.")
    p.add_argument("gds")
    p.add_argument("--seg", required=True,
                   help="comma-separated list of edge lengths (um)")
    p.add_argument("--out", default="scdc")
    p.add_argument("--order", type=float, default=1.0,
                   help="assumed order when only two seg values are given")
    a, rest = p.parse_known_args(argv)

    segs = sorted({float(s) for s in a.seg.split(",")}, reverse=True)
    if len(segs) < 2:
        p.error("give at least two seg values")
    s_min = segs[-1]

    fixed = []
    if not _has_opt(rest, "--seg-far"):
        fixed += ["--seg-far", f"{5.0 * s_min:g}"]
    if not _has_opt(rest, "--seg-max"):
        fixed += ["--seg-max", f"{30.0 * s_min:g}"]
    if not _has_opt(rest, "--fine-radius"):
        fixed += ["--fine-radius", f"{200.0 * s_min / 1.5:g}"]
    if fixed:
        print("  mesh parameters derived from the smallest seg and held "
              "fixed: " + " ".join(fixed))
    if not _has_opt(rest, "--grid"):
        print("  WARNING: --grid not given, the pFFT grid will follow the "
              "median triangle size and change between runs")
    if _has_opt(rest, "--mesh-only"):
        p.error("--mesh-only is incompatible with a convergence study")

    runs = []
    for s in segs:
        prefix = f"{a.out}_seg{s:g}"
        argv_i = [a.gds, "--seg", f"{s:g}", "--out", prefix] + rest + fixed
        print("\n" + "=" * 72)
        print(f"  seg = {s:g} um   ->   {prefix}_*")
        print("=" * 72)
        t0 = time.time()
        r = _extract(_capture_run(argv_i))
        r["seg"] = s
        r["time"] = time.time() - t0
        r["prefix"] = prefix
        runs.append(r)
        if len(runs) > 1 and runs[-1]["M"].shape != runs[0]["M"].shape:
            raise RuntimeError("M changed shape between runs, the layer-3 or "
                               "layer-5 polygons are not resolved identically")

    h = np.array([r["seg"] for r in runs])
    Ms = np.array([r["M"] for r in runs])            # (nrun, ni, nj)
    ni, nj = Ms.shape[1:]
    P = np.full((ni, nj), np.nan)
    OK = np.zeros((ni, nj), bool)
    Minf = np.empty((ni, nj))
    for i in range(ni):
        for j in range(nj):
            m = Ms[:, i, j]
            if len(h) >= 3:
                pij, ok = fit_order(h[-3:], m[-3:])
            else:
                pij, ok = np.nan, False
            P[i, j], OK[i, j] = pij, ok
            p_use = pij if ok else a.order
            Minf[i, j] = extrapolate(h, m, p_use)
    lam_inf = eig_inverse(Minf)
    lam_fin = eig_inverse(Ms[-1])

    # ------------------------------------------------------------ report
    src, surf = runs[-1]["src"], runs[-1]["surf"]
    L = []
    L.append(f"# scdc mesh-convergence report: {a.gds}\n")
    L.append(f"Generated {time.strftime('%Y-%m-%d %H:%M')}. Edge length "
             f"`--seg` varied over {', '.join(f'{s:g}' for s in segs)} um, "
             "all other mesh parameters fixed"
             + (" (" + " ".join(fixed) + ")" if fixed else "") + ".\n")
    L.append("Model M_ij(h) = M_ij(inf) + a_ij h^p_ij, p_ij fitted on the "
             "three finest meshes by solving\n")
    L.append("    (M(h1) - M(h2)) / (M(h2) - M(h3)) = (h1^p - h2^p) / (h2^p - h3^p)\n")
    L.append(f"then M(inf) = M(h3) - a h3^p. Where the sequence is not monotone "
             f"or no order in [0.05, 4] fits, the assumed order p = {a.order:g} "
             "is used and the entry is flagged.\n")

    L.append("## Runs\n")
    L.append("| seg (um) | triangles | h_min | h_median (pFFT grid if --grid absent) | "
             + " | ".join(f"L({s}) = 2E/I^2 (pH)" for s in src)
             + " | CG iterations | time (s) |")
    L.append("|---" * (5 + len(src)) + "|---|")
    for r in runs:
        # L needs the current, unknown here, so print 2E in pJ-equivalent
        # only when --current is the default 1 mA
        Lp = " | ".join(f"{2 * e / 1e-6 * 1e12:.4f}" for e in r["E"])
        L.append(f"| {r['seg']:g} | {r['n_tri']} | {r['h_min']:.3g} | "
                 f"{r['h_med']:.3g} | {Lp} | "
                 f"{'/'.join(str(k) for k in r['iters'])} | {r['time']:.0f} |")
    L.append("\nL = 2E/I^2 assumes I = 1 mA. E decreases monotonically with "
             "refinement when the problem is well posed (variational bound).\n")

    L.append("## M(h) per coefficient (pH)\n")
    hdr = "| surface | source | " + " | ".join(f"h = {s:g}" for s in segs) \
          + " | p | M(inf) | last step (%) | extrapolation (%) |"
    L.append(hdr)
    L.append("|---" * (len(segs) + 6) + "|")
    for i in range(ni):
        for j in range(nj):
            m = Ms[:, i, j] * 1e12
            mi = Minf[i, j] * 1e12
            step = 100 * (m[-1] - m[-2]) / abs(m[-1]) if m[-1] else np.nan
            ext = 100 * (mi - m[-1]) / abs(m[-1]) if m[-1] else np.nan
            ps = f"{P[i, j]:.2f}" if OK[i, j] else f"{a.order:g} (assumed)"
            L.append(f"| {surf[i]} | {src[j]} | "
                     + " | ".join(f"{x:.6g}" for x in m)
                     + f" | {ps} | {mi:.6g} | {step:+.2f} | {ext:+.2f} |")
    L.append("\n`last step` is the relative change between the two finest "
             "meshes, `extrapolation` the relative correction from the finest "
             "mesh to M(inf). Coefficients with |M| small compared to the "
             "diagonal carry large relative errors and their order is not "
             "meaningful.\n")

    L.append("## Extrapolated matrix M(inf) (pH)\n")
    L.append("| surface \\ source | " + " | ".join(src) + " |")
    L.append("|---" * (nj + 1) + "|")
    for i in range(ni):
        L.append(f"| {surf[i]} | " + " | ".join(f"{x*1e12:.6g}" for x in Minf[i])
                 + " |")

    if lam_inf is not None:
        L.append("\n## Eigenvalues of M^-1 (mA/Phi0)\n")
        L.append("Sorted by increasing modulus, finest mesh against "
                 "extrapolated matrix.\n")
        L.append("| k | finest mesh | M(inf) | change (%) |")
        L.append("|---|---|---|---|")
        for k, (lf, li) in enumerate(zip(lam_fin, lam_inf)):
            ch = 100 * (abs(li) - abs(lf)) / abs(lf)
            L.append(f"| {k} | {_fmt_c(lf)} | {_fmt_c(li)} | {ch:+.2f} |")
        cond = np.linalg.cond(Minf)
        L.append(f"\nCondition number of M(inf): {cond:.4g}.\n")
    else:
        L.append("\nM is not square or is singular, eigenvalues of M^-1 not "
                 "computed.\n")

    L.append("## Output files\n")
    for r in runs:
        L.append(f"- `{r['prefix']}_*` (seg = {r['seg']:g} um)")
    L.append(f"- `{a.out}_convergence.npz`")
    L.append(f"- `{a.out}_convergence.md`")

    md = "\n".join(L) + "\n"
    with open(f"{a.out}_convergence.md", "w") as f:
        f.write(md)
    np.savez_compressed(f"{a.out}_convergence.npz", seg=h, M=Ms, p=P,
                        p_ok=OK, M_inf=Minf, n_tri=[r["n_tri"] for r in runs],
                        E=np.array([r["E"] for r in runs]),
                        sources=np.array(src), surfaces=np.array(surf))
    print("\n" + md)
    print(f"written: {a.out}_convergence.md, {a.out}_convergence.npz")


if __name__ == "__main__":
    main()
