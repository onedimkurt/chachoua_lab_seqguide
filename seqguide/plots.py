"""Static figures for the HTML report (matplotlib, embedded as PNG)."""

from __future__ import annotations

import base64
import io

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrow, Patch  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"
SERIES = "#2a78d6"
FEATURE_COLORS = {"CDS": "#2a78d6", "promoter": "#eb6834", "rep_origin": "#1baf7a", "terminator": "#6250d6"}
FEATURE_OTHER = "#a3a29c"
SEVERITY = {"high": ("#d03b3b", "v", "Likely disrupts function"),
            "moderate": ("#ec835a", "D", "May affect function"),
            "low": ("#7d7c76", "o", "Probably harmless")}

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
})


def to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _depth_axis(ax, depth, min_depth, length):
    x = np.arange(1, len(depth) + 1)
    ax.fill_between(x, depth, step="mid", color=SERIES, alpha=0.25, linewidth=0)
    ax.plot(x, depth, color=SERIES, linewidth=1.2, drawstyle="steps-mid")
    if min_depth:
        ax.axhline(min_depth, color=INK2, linewidth=0.8, linestyle=(0, (4, 3)))
        ax.text(length, min_depth, f" minimum {min_depth}×", va="bottom", ha="right", color=INK2, fontsize=8)
    ax.set_ylabel("Reads covering\neach position")
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_xlim(0, length)


def plasmid_figure(ref, variants, depth=None, min_depth=30) -> str:
    L = len(ref.seq)
    rows = 2 if depth is not None else 1
    fig, axes = plt.subplots(rows, 1, figsize=(10, 1.9 + 1.6 * rows), sharex=True,
                             gridspec_kw={"height_ratios": [1.6, 1.4][:rows] if rows == 2 else [1]})
    axes = np.atleast_1d(axes)
    ax = axes[0]
    ax.set_xlim(0, L)
    ax.set_ylim(-1.6, 2.0)
    ax.axhline(0, color=INK2, linewidth=1.2)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    used_types = set()
    label_rows = {1: [], -1: []}  # per strand: right edge of the last label on each row
    for i, f in enumerate(sorted(ref.features, key=lambda f: f.start)):
        if f.type in ("gene", "source", "misc_feature") and any(g.type == "CDS" for g in ref.features):
            continue
        color = FEATURE_COLORS.get(f.type, FEATURE_OTHER)
        used_types.add(f.type if f.type in FEATURE_COLORS else "other")
        y = 0.45 if f.strand >= 0 else -0.45
        for s, e, strand in f.parts:
            w = e - s
            head = min(w * 0.35, L * 0.012)
            if strand >= 0:
                ax.add_patch(FancyArrow(s, y, w, 0, width=0.42, head_width=0.42, head_length=head,
                                        length_includes_head=True, color=color, linewidth=0))
            else:
                ax.add_patch(FancyArrow(e, y, -w, 0, width=0.42, head_width=0.42, head_length=head,
                                        length_includes_head=True, color=color, linewidth=0))
        mid = (f.start + f.end) / 2
        if (f.end - f.start) > L * 0.035:
            half = len(f.label) * 0.0042 * L
            lrows = label_rows[1 if y > 0 else -1]
            level = next((k for k, edge in enumerate(lrows) if mid - half > edge), len(lrows))
            if level < 2:
                if level == len(lrows):
                    lrows.append(0)
                lrows[level] = mid + half
                off = 0.42 + 0.32 * level
                ax.text(mid, y + (off if y > 0 else -off), f.label, ha="center", va="bottom" if y > 0 else "top",
                        fontsize=8, color=INK)
    # variants as markers above the line
    for v in variants:
        color, marker, _ = SEVERITY.get(v.severity or "low", SEVERITY["low"])
        x = v.ref_pos + max(len(v.ref), 1) / 2
        if v.structural or v.size > 30:
            ax.add_patch(plt.Rectangle((v.ref_pos, 1.47), max(len(v.ref), L * 0.004), 0.25, color=color, linewidth=0))
        ax.plot([x], [1.62], marker=marker, color=color, markersize=8, markeredgecolor=SURFACE, markeredgewidth=1.5)
        if len(variants) <= 8 and v.severity in ("high", "moderate"):
            lbl = v.describe()
            lbl = lbl if len(lbl) <= 22 else lbl[:20] + "…"
            ax.text(x, 1.76, lbl, ha="center", va="bottom", fontsize=7.5, color=INK2)
    handles = [Patch(color=FEATURE_COLORS.get(t, FEATURE_OTHER),
                     label={"CDS": "Protein-coding gene", "rep_origin": "Origin of replication",
                            "promoter": "Promoter", "terminator": "Terminator", "other": "Other feature"}.get(t, t))
               for t in ["CDS", "promoter", "rep_origin", "terminator", "other"] if t in used_types]
    sev_used = {v.severity or "low" for v in variants}
    handles += [Line2D([], [], marker=SEVERITY[s][1], color=SEVERITY[s][0], linestyle="", markersize=7, label=SEVERITY[s][2])
                for s in ("high", "moderate", "low") if s in sev_used]
    if handles:
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.32 if rows == 1 else 1.28),
                  ncol=min(4, len(handles)), frameon=False, fontsize=8)
    if rows == 2:
        _depth_axis(axes[1], np.asarray(depth), min_depth, L)
    axes[-1].set_xlabel(f"Position on reference {ref.name} (bp)")
    return to_b64(fig)


def coverage_figure(depth, name, min_depth=30) -> str:
    fig, ax = plt.subplots(figsize=(10, 2.2))
    _depth_axis(ax, np.asarray(depth), min_depth, len(depth))
    ax.set_xlabel(f"Position on {name} (bp)")
    return to_b64(fig)


def targets_figure(depth: dict, refs: list, variants: dict, min_depth: int, max_panels: int = 12) -> str:
    refs = refs[:max_panels]
    n = len(refs)
    cols = 1 if n == 1 else 2
    nrows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(nrows, cols, figsize=(10, 1.9 * nrows), squeeze=False)
    for ax, ref in zip(axes.flat, refs):
        d = depth.get(ref.name, np.zeros(len(ref.seq)))
        _depth_axis(ax, d, min_depth, len(ref.seq))
        ax.set_title(ref.name, loc="left", fontsize=9, color=INK)
        top = max(1, float(np.max(d)) if len(d) else 1)
        for v in variants.get(ref.name, []):
            color, marker, _ = SEVERITY.get(v.severity or "low", SEVERITY["low"])
            ax.plot([v.ref_pos + 1], [top * 1.05], marker=marker, color=color, markersize=7, markeredgecolor=SURFACE)
        ax.set_ylim(0, top * 1.15)
        ax.set_ylabel("Depth")
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    fig.tight_layout()
    return to_b64(fig)


def length_hist_figure(files: list) -> str:
    fig, ax = plt.subplots(figsize=(10, 2.2))
    colors = ["#2a78d6", "#eb6834", "#1baf7a"]
    for i, f in enumerate(files[:3]):
        h = f.get("length_hist") or {}
        if not h.get("counts"):
            continue
        edges = np.asarray(h["edges"])
        ax.stairs(h["counts"], edges, fill=False, color=colors[i], linewidth=1.6, label=f["name"])
    ax.set_xlabel("Read length (bp)")
    ax.set_ylabel("Number of reads")
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    if len(files) > 1:
        ax.legend(frameon=False, fontsize=8)
    return to_b64(fig)
