import csv
import json
import os
import random
import subprocess
import sys
from _path import ROOT, RESULTS
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from latsieve.common import nfs_box, next_prime
from latsieve.reduction import lll3w
from latsieve.zfs import select_frame, zonotope_fibers, zfs_frame

FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.6, "font.family": "serif", "mathtext.fontset": "cm", "axes.spines.top": False, "axes.spines.right": False})
COLORS = {"ZFS": "#2a78d6", "LS": "#eb6834", "PS": "#1baf7a", "GGV": "#eda100", "HAKT": "#e87ba4", "NTV": "#008300", "MR": "#6250d6", "SE": "#e34948"}
MARK = {"ZFS": "o", "LS": "s", "PS": "D", "GGV": "^", "HAKT": "v", "NTV": "P", "MR": "X", "SE": "*"}
ORDER = ["ZFS", "LS", "PS", "GGV", "HAKT", "NTV", "MR", "SE"]
HEUR = {"GGV", "HAKT", "NTV", "MR"}
GRID = "#d9d8d4"


def save(fig, name, tight=True):
    eps = os.path.join(FIG, name + ".eps")
    if tight:
        fig.savefig(eps, format="eps", bbox_inches="tight", pad_inches=0.04)
    else:
        fig.savefig(eps, format="eps")
    plt.close(fig)
    subprocess.run(["epstopdf", eps], check=True)


def read(name):
    with open(os.path.join(RESULTS, name)) as fh:
        return list(csv.DictReader(fh))


def fig_benchmark():
    rows = read("bench_H64.csv")
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.35))
    for m in ORDER:
        pts = sorted((int(r["log2r"]), float(r["time_ms"]), float(r["completeness"]), int(r["peak_bytes"])) for r in rows if r["method"] == m)
        x = [p[0] for p in pts]
        ls = "--" if m in HEUR else "-"
        lw = 1.8 if m == "ZFS" else 1.0
        axs[0].plot(x, [p[1] for p in pts], ls, color=COLORS[m], marker=MARK[m], ms=3.2, lw=lw, label=m)
        axs[1].plot(x, [100 * p[2] for p in pts], ls, color=COLORS[m], marker=MARK[m], ms=3.2, lw=lw, label=m)
        axs[2].plot(x, [p[3] / 1024 for p in pts], ls, color=COLORS[m], marker=MARK[m], ms=3.2, lw=lw, label=m)
    axs[0].set_yscale("log")
    axs[0].set_ylabel("time per lattice (ms)")
    axs[1].set_ylabel("points found (%)")
    axs[2].set_yscale("log")
    axs[2].set_ylabel("peak working memory (KiB)")
    for a, t in zip(axs, ("(a)", "(b)", "(c)")):
        a.set_xlabel(r"$\log_2 r$")
        a.set_title(t, loc="left")
        a.grid(True, color=GRID, lw=0.4)
        a.set_xticks(range(4, 27, 4))
    for a in (7, 14):
        for ax in axs:
            ax.axvline(a, color="#888780", lw=0.6, ls=":")
    axs[1].set_ylim(-3, 103)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=8, frameon=False, bbox_to_anchor=(0.5, 1.07), handlelength=2.2, columnspacing=1.0)
    fig.tight_layout()
    save(fig, "benchmark")


def fig_scaling():
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.2))
    for a, H in zip(axs, (32, 64, 128)):
        rows = read(f"bench_H{H}.csv")
        for m in ORDER:
            pts = sorted((int(r["log2r"]), float(r["time_ms"])) for r in rows if r["method"] == m)
            ls = "--" if m in HEUR else "-"
            a.plot([p[0] for p in pts], [p[1] for p in pts], ls, color=COLORS[m], marker=MARK[m], ms=2.8, lw=1.6 if m == "ZFS" else 0.9, label=m)
        a.set_yscale("log")
        a.set_xlabel(r"$\log_2 r$")
        a.set_title(f"({'abc'[(32, 64, 128).index(H)]}) $H={2*H}\\times{2*H}\\times{H}$", loc="left")
        a.grid(True, color=GRID, lw=0.4)
        a.set_xticks(range(4, 27, 4))
    axs[0].set_ylabel("time per lattice (ms)")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=8, frameon=False, bbox_to_anchor=(0.5, 1.08))
    fig.tight_layout()
    save(fig, "scaling")


def fig_expected():
    rows = read("expected_nodes.csv")
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.3))
    fcol = {"ZFS": COLORS["ZFS"], "LLL": "#6250d6", "LS": COLORS["LS"], "PS": COLORS["PS"]}
    fmk = {"ZFS": "o", "LLL": "^", "LS": "s", "PS": "D"}
    for f in ("ZFS", "LLL", "LS", "PS"):
        sel = [r for r in rows if r["frame"] == f]
        axs[0].loglog([float(r["pred_planes"]) for r in sel], [float(r["meas_planes"]) for r in sel], fmk[f], color=fcol[f], ms=3.5, ls="none", label=f)
        axs[1].loglog([float(r["pred_lines"]) for r in sel], [float(r["meas_lines"]) for r in sel], fmk[f], color=fcol[f], ms=3.5, ls="none", label=f)
    for a, lab in ((axs[0], "planes"), (axs[1], "lines")):
        lo, hi = 0.02, 3000
        a.plot([lo, hi], [lo, hi], color="#52514e", lw=0.6)
        a.set_xlabel(f"predicted mean {lab}")
        a.set_ylabel(f"measured mean {lab}")
        a.grid(True, color=GRID, lw=0.4)
    axs[0].set_title("(a)", loc="left")
    axs[1].set_title("(b)", loc="left")
    axs[0].legend(frameon=False, loc="upper left")
    sel = [r for r in rows if r["frame"] == "ZFS"]
    E = [float(r["pred_points"]) for r in sel]
    tot = [float(r["meas_planes"]) + float(r["meas_lines"]) + float(r["meas_points"]) for r in sel]
    axs[2].loglog(E, tot, "o", color=COLORS["ZFS"], ms=3.5, label="ZFS nodes, measured")
    xs = [10 ** (k / 20) for k in range(-40, 81)]
    axs[2].loglog(xs, [x + 2.02 * x ** (2 / 3) + 2.02 * x ** (1 / 3) for x in xs], color="#52514e", lw=0.8, label=r"$E+2.02E^{2/3}+2.02E^{1/3}$")
    axs[2].loglog(xs, xs, color="#888780", lw=0.6, ls=":", label=r"$E$")
    axs[2].set_xlabel(r"density $E$")
    axs[2].set_ylabel("mean number of nodes")
    axs[2].set_title("(c)", loc="left")
    axs[2].grid(True, color=GRID, lw=0.4)
    axs[2].legend(frameon=False, loc="lower right")
    fig.tight_layout()
    save(fig, "expected")


def fig_transitions():
    part = json.load(open(os.path.join(RESULTS, "transitions_partition.json")))
    hist = read("transitions_hist.csv")
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.45), gridspec_kw={"width_ratios": [1, 1, 1.35]})
    ramp = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
    for a, tag, title in ((axs[0], "family", "(a) diagonal family"), (axs[1], "random", "(b) a random lattice")):
        g = part[tag]["grid"]
        n = len(part[tag]["vectors"])
        cmap = ListedColormap([ramp[round(i * (len(ramp) - 1) / max(1, n - 1))] for i in range(n)])
        I = len(g)
        J = len(g[0])
        a.imshow([[g[x][y] for x in range(I)] for y in range(J)], origin="lower", cmap=cmap, vmin=-0.5, vmax=n - 0.5, interpolation="nearest")
        for x in range(I):
            for y in range(J):
                v = g[x][y]
                a.text(x, y, str(v + 1), ha="center", va="center", fontsize=5, color="#0b0b0b" if v < n / 2 else "#ffffff")
        a.set_title(f"{title}, $|T|={n}$", loc="left")
        a.set_xlabel(r"$x_0$")
        a.set_ylabel(r"$x_1$")
        a.set_xticks(range(0, I, 3))
        a.set_yticks(range(0, J, 3))
        for s in a.spines.values():
            s.set_visible(False)
    styles = {16: ("#86b6ef", "o"), 32: ("#2a78d6", "s"), 64: ("#104281", "D")}
    for I in (16, 32, 64):
        cnt = {}
        tot = 0
        for r in hist:
            if int(r["I"]) == I:
                k = int(r["T"])
                cnt[k] = cnt.get(k, 0) + int(r["count"])
                tot += int(r["count"])
        ks = sorted(cnt)
        axs[2].plot(ks, [100 * cnt[k] / tot for k in ks], "-", marker=styles[I][1], color=styles[I][0], ms=3, lw=1.0, label=f"$I=J={I}$")
    axs[2].axvline(3, color="#888780", lw=0.6, ls=":")
    axs[2].text(3.3, 18, "three gap\nbound", fontsize=6, color="#52514e")
    axs[2].set_xlabel(r"number of transition vectors $|T|$")
    axs[2].set_ylabel("lattices (%)")
    axs[2].set_title("(c) distribution over random lattices", loc="left")
    axs[2].grid(True, color=GRID, lw=0.4)
    axs[2].legend(frameon=False, loc="upper right")
    axs[2].set_xlim(2, 26)
    fig.tight_layout()
    save(fig, "transitions")


def fig_factorbase():
    rows = read("factorbase_H64_B262144.csv")
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.2))
    names = [m for m in ORDER if any(r["method"] == m for r in rows)]
    t = [float(next(r for r in rows if r["method"] == m)["total_time_s"]) for m in names]
    rec = [100 * float(next(r for r in rows if r["method"] == m)["candidate_recall"]) for m in names]
    mass = [100 * float(next(r for r in rows if r["method"] == m)["hit_mass_fraction"]) for m in names]
    y = list(range(len(names)))
    axs[0].barh(y, t, color=[COLORS[m] for m in names], height=0.62)
    axs[0].set_xscale("log")
    axs[0].set_yticks(y)
    axs[0].set_yticklabels(names)
    axs[0].invert_yaxis()
    axs[0].set_xlabel("total sieving time (s)")
    for i, v in enumerate(t):
        axs[0].text(v * 1.08, i, f"{v:.2f}", va="center", fontsize=6.5, color="#0b0b0b")
    axs[0].set_xlim(0.3, 300)
    axs[0].set_title("(a)", loc="left")
    w = 0.38
    axs[1].barh([i - w / 2 for i in y], mass, height=w, color="#2a78d6", label="hits found")
    axs[1].barh([i + w / 2 for i in y], rec, height=w, color="#eb6834", label="candidate cells found")
    axs[1].set_yticks(y)
    axs[1].set_yticklabels(names)
    axs[1].invert_yaxis()
    axs[1].set_xlabel("fraction of the exact sieve (%)")
    axs[1].set_xlim(0, 112)
    axs[1].legend(frameon=False, loc="lower center", bbox_to_anchor=(0.55, 1.0), ncol=2)
    axs[1].set_title("(b)", loc="left")
    for a in axs:
        a.grid(True, axis="x", color=GRID, lw=0.4)
    fig.tight_layout()
    save(fig, "factorbase")


def plane_polygon(n, val, lo, hi):
    corners = [(X, Y, Z) for X in (lo[0], hi[0]) for Y in (lo[1], hi[1]) for Z in (lo[2], hi[2])]
    pts = []
    for p in corners:
        for q in corners:
            if p < q and sum(p[i] != q[i] for i in range(3)) == 1:
                fp = sum(n[i] * p[i] for i in range(3)) - val
                fq = sum(n[i] * q[i] for i in range(3)) - val
                if fp == 0:
                    pts.append(tuple(float(v) for v in p))
                if (fp < 0 < fq) or (fq < 0 < fp):
                    t = fp / (fp - fq)
                    pts.append(tuple(p[i] + t * (q[i] - p[i]) for i in range(3)))
    if len(pts) < 3:
        return []
    cx = [sum(p[i] for p in pts) / len(pts) for i in range(3)]
    import math
    nn = math.sqrt(sum(v * v for v in n))
    un = [v / nn for v in n]
    ref = [1.0, 0.0, 0.0] if abs(un[0]) < 0.9 else [0.0, 1.0, 0.0]
    e1 = [ref[1] * un[2] - ref[2] * un[1], ref[2] * un[0] - ref[0] * un[2], ref[0] * un[1] - ref[1] * un[0]]
    e2 = [un[1] * e1[2] - un[2] * e1[1], un[2] * e1[0] - un[0] * e1[2], un[0] * e1[1] - un[1] * e1[0]]
    def ang(p):
        d = [p[i] - cx[i] for i in range(3)]
        return math.atan2(sum(d[i] * e2[i] for i in range(3)), sum(d[i] * e1[i] for i in range(3)))
    pts = sorted(set(pts), key=ang)
    return pts + [pts[0]]


def fig_3d():
    box = nfs_box(16, 16, 32)
    lo, hi = box.lo, box.hi
    r, l0, l1 = 1049, 130, 908
    basis = ((r, 0, 0), (l0, 1, 0), (l1, 0, 1))
    L = [hi[i] - lo[i] for i in range(3)]
    a, c, e = zfs_frame(r, l0, l1, box)
    e, klo, khi, slabs = zonotope_fibers(a, c, e, lo, hi)
    n = (a[1] * c[2] - a[2] * c[1], a[2] * c[0] - a[0] * c[2], a[0] * c[1] - a[1] * c[0])
    D = sum(n[i] * e[i] for i in range(3))
    fig = plt.figure(figsize=(7.0, 3.3))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    corners = [(X, Y, Z) for X in (lo[0], hi[0]) for Y in (lo[1], hi[1]) for Z in (lo[2], hi[2])]
    edges = [(p, q) for p in corners for q in corners if sum(p[i] != q[i] for i in range(3)) == 1 and p < q]
    for axx in (ax, ax2):
        axx.add_collection3d(Line3DCollection(edges, colors="#888780", linewidths=0.6))
    pal = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e34948", "#eda100"]
    allpts = []
    nlines = 0
    for idx, k in enumerate(range(klo, khi + 1)):
        col = pal[idx % len(pal)]
        poly = plane_polygon(n, D * k, lo, hi)
        if poly:
            ax2.plot([p[0] for p in poly], [p[1] for p in poly], [p[2] for p in poly], color=col, lw=0.8, ls="--")
        jl = -(1 << 62)
        jh = 1 << 62
        for f, tc, g, h in slabs:
            t = tc * k
            jl = max(jl, -((t - g) // f))
            jh = min(jh, (h - t) // f)
        for j in range(jl, jh + 1):
            p = [j * c[i] + k * e[i] for i in range(3)]
            sl = -(1 << 62)
            sh = 1 << 62
            for i in range(3):
                if a[i] > 0:
                    sl = max(sl, -((p[i] - lo[i]) // a[i]))
                    sh = min(sh, (hi[i] - p[i]) // a[i])
                elif a[i] < 0:
                    sl = max(sl, -((hi[i] - p[i]) // (-a[i])))
                    sh = min(sh, (p[i] - lo[i]) // (-a[i]))
            nlines += 1
            if sl <= sh:
                pts = [tuple(p[i] + s * a[i] for i in range(3)) for s in range(sl, sh + 1)]
                allpts += [(q, col) for q in pts]
                if len(pts) > 1:
                    ax2.plot([pts[0][0], pts[-1][0]], [pts[0][1], pts[-1][1]], [pts[0][2], pts[-1][2]], color=col, lw=1.2)
    xs = [p[0][0] for p in allpts]
    ys = [p[0][1] for p in allpts]
    zs = [p[0][2] for p in allpts]
    ax.scatter(xs, ys, zs, s=9, c="#0b0b0b", depthshade=False)
    ax2.scatter(xs, ys, zs, s=9, c=[p[1] for p in allpts], depthshade=False)
    for axx, title in ((ax, "(a) the points of $\\Lambda$ in $H$"), (ax2, "(b) planes and lines visited by ZFS")):
        axx.set_xlabel("$x_0$", labelpad=-5)
        axx.set_ylabel("$x_1$", labelpad=-5)
        axx.set_zlabel("$x_2$", labelpad=-6)
        axx.tick_params(pad=-3, labelsize=6)
        axx.view_init(elev=16, azim=-62)
        axx.set_title(title, loc="left", fontsize=8)
        axx.set_box_aspect((1, 1, 1))
    fig.subplots_adjust(left=0.0, right=0.93, wspace=0.1, top=0.93, bottom=0.03)
    save(fig, "lattice3d", tight=False)
    return r, l0, l1, a, c, e, khi - klo + 1, nlines, len(allpts)


def fig_nodes():
    rows = read("nodes_H64.csv")
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.35))
    lab = {"LS": "LS frame", "PS": "PS frame", "SE": "SE tree", "ZFS": "ZFS frame"}
    for m in ("ZFS", "LS", "PS", "SE"):
        sel = [r for r in rows if r["method"] == m]
        x = [int(r["log2r"]) for r in sel]
        E = [float(r["E"]) for r in sel]
        tot = [float(r["points"]) + float(r["lines"]) + float(r["planes"]) for r in sel]
        axs[0].semilogy(x, tot, marker=MARK[m], color=COLORS[m], lw=1.2, ms=3.5, label=lab[m])
        axs[1].loglog(E, [t / e for t, e in zip(tot, E)], marker=MARK[m], color=COLORS[m], lw=1.2, ms=3.5, label=lab[m])
    sel = [r for r in rows if r["method"] == "ZFS"]
    x = [int(r["log2r"]) for r in sel]
    E = [float(r["E"]) for r in sel]
    axs[0].semilogy(x, E, color="#888780", lw=0.8, ls=":", label=r"$E$")
    xs = [10 ** (k / 20) for k in range(-40, 101)]
    axs[1].loglog(xs, [1 + 2.02 * v ** (-1 / 3) + 2.02 * v ** (-2 / 3) for v in xs], color="#52514e", lw=0.8, ls="--", label=r"$1+2.02E^{-1/3}+2.02E^{-2/3}$")
    for v in (7, 14):
        axs[0].axvline(v, color="#888780", lw=0.6, ls=":")
    axs[0].set_xlabel(r"$\log_2 r$")
    axs[0].set_ylabel("expected nodes per lattice")
    axs[0].set_xticks(range(4, 27, 4))
    axs[1].set_xlabel(r"density $E$")
    axs[1].set_ylabel("expected nodes per point")
    axs[0].set_title("(a)", loc="left")
    axs[1].set_title("(b)", loc="left")
    for a in axs:
        a.grid(True, color=GRID, lw=0.4)
    h0, l0 = axs[0].get_legend_handles_labels()
    h1, l1 = axs[1].get_legend_handles_labels()
    fig.legend(h0 + h1[-1:], l0 + l1[-1:], loc="upper center", ncol=6, frameon=False, bbox_to_anchor=(0.5, 1.07), handlelength=2.2, columnspacing=1.2)
    fig.tight_layout()
    save(fig, "nodes")


if __name__ == "__main__":
    which = sys.argv[1:] or ["benchmark", "scaling", "expected", "nodes", "transitions", "factorbase", "3d"]
    for w in which:
        if w == "benchmark":
            fig_benchmark()
        elif w == "scaling":
            fig_scaling()
        elif w == "expected":
            fig_expected()
        elif w == "nodes":
            fig_nodes()
        elif w == "transitions":
            fig_transitions()
        elif w == "factorbase":
            fig_factorbase()
        elif w == "3d":
            print(fig_3d())
