from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, PathPatch, Polygon, Rectangle
from matplotlib.path import Path as MplPath
import numpy as np


OUT_DIR = Path(__file__).resolve().parent


def arrow(ax, start, end, color, lw=1.5, mutation_scale=12, alpha=1.0, rad=0.0):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=mutation_scale,
            linewidth=lw,
            color=color,
            alpha=alpha,
            connectionstyle=f"arc3,rad={rad}",
        )
    )


def plume_patch(x0, y0, width, height, color, alpha):
    verts = [
        (x0, y0),
        (x0 - 0.35 * width, y0 + 0.30 * height),
        (x0 - 0.55 * width, y0 + 0.70 * height),
        (x0 - 0.15 * width, y0 + height),
        (x0 + 0.20 * width, y0 + height),
        (x0 + 0.55 * width, y0 + 0.70 * height),
        (x0 + 0.35 * width, y0 + 0.30 * height),
        (x0, y0),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE3,
        MplPath.CURVE3,
        MplPath.CURVE3,
        MplPath.LINETO,
        MplPath.CURVE3,
        MplPath.CURVE3,
        MplPath.CLOSEPOLY,
    ]
    return PathPatch(MplPath(verts, codes), facecolor=color, edgecolor="none", alpha=alpha)


def main():
    fig, ax = plt.subplots(figsize=(11.2, 6.4))
    fig.patch.set_facecolor("white")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")

    # Atmosphere and haze.
    ax.add_patch(Rectangle((0, 4.62), 10, 1.38, facecolor="#dbe8ef", edgecolor="none"))
    haze_x = np.linspace(0, 10, 300)
    for y, amp, color, alpha in [
        (5.55, 0.035, "#9eb3c0", 0.55),
        (5.20, 0.028, "#b7c6cf", 0.50),
        (4.90, 0.020, "#c7d4db", 0.70),
    ]:
        haze_y = y + amp * np.sin(1.7 * haze_x + y)
        ax.plot(haze_x, haze_y, color=color, lw=2.0, alpha=alpha)

    # Hydrocarbon pond.
    surface = 4.10
    bottom = 1.38
    water_x = np.linspace(0.25, 9.75, 400)
    wave = surface + 0.05 * np.sin(2.9 * water_x) + 0.025 * np.sin(7.1 * water_x)
    pond_poly = [(0.25, bottom), *zip(water_x, wave), (9.75, bottom)]
    ax.add_patch(Polygon(pond_poly, closed=True, facecolor="#263f67", edgecolor="none", alpha=0.96))
    ax.fill_between(water_x, wave - 0.08, wave + 0.005, color="#8fb7d6", alpha=0.45)

    # Vertical concentration bands in liquid.
    for i, alpha in enumerate(np.linspace(0.10, 0.35, 7)):
        y = bottom + 0.28 + i * 0.34
        ax.add_patch(Rectangle((0.35, y), 9.3, 0.08, facecolor="#69b3c9", edgecolor="none", alpha=alpha))

    # Sediment, evaporite, and pore bed.
    ax.add_patch(Rectangle((0, 0), 10, bottom, facecolor="#4b3a2d", edgecolor="none"))
    ax.add_patch(Rectangle((0.25, bottom - 0.16), 9.5, 0.18, facecolor="#c8b26b", edgecolor="none", alpha=0.86))
    rng = np.random.default_rng(4)
    for _ in range(85):
        x = rng.uniform(0.25, 9.75)
        y = rng.uniform(0.18, bottom - 0.18)
        r = rng.uniform(0.025, 0.08)
        ax.add_patch(Circle((x, y), r, facecolor=rng.choice(["#6b5846", "#7e6b53", "#3f3329"]), edgecolor="none", alpha=0.95))
    for x in np.linspace(0.8, 9.1, 9):
        ax.plot([x, x + 0.20 * np.sin(x)], [bottom - 0.15, bottom - 0.78], color="#2d241d", lw=1.0, alpha=0.55)

    # Cryovolcanic pulse plume.
    ax.add_patch(plume_patch(7.25, 0.75, 1.65, 3.75, "#7fd3e6", 0.38))
    ax.add_patch(plume_patch(7.25, 0.75, 0.95, 3.28, "#e8f7fb", 0.52))
    ax.add_patch(Circle((7.25, 0.72), 0.42, facecolor="#d8f5fb", edgecolor="#287d8e", lw=1.1, alpha=0.95))
    for y in np.linspace(1.30, 4.00, 7):
        ax.plot([6.75, 7.75], [y, y + 0.08 * np.sin(3 * y)], color="#e7fbff", lw=1.2, alpha=0.55)

    # Atmospheric deposition arrows.
    for x, rad in [(1.5, -0.10), (2.7, 0.06), (4.1, -0.04), (5.4, 0.08)]:
        arrow(ax, (x, 5.22), (x + 0.18, 4.28), "#8a5a1f", lw=1.1, mutation_scale=10, alpha=0.75, rad=rad)
    ax.text(2.9, 5.38, "organic haze, nitriles,\nhydrocarbon rainout", ha="center", va="center", fontsize=9.2, color="#4a3b2b")

    # Evaporation arrows.
    for x in [3.6, 4.6, 5.8]:
        arrow(ax, (x, 4.25), (x - 0.18, 4.88), "#527b99", lw=1.0, mutation_scale=9, alpha=0.65, rad=0.05)
    ax.text(5.0, 4.95, "evaporation\nconcentrates solutes", ha="center", va="bottom", fontsize=8.6, color="#31536a")

    # Surface chemistry arrows.
    for x in [1.4, 2.3, 3.2, 4.1]:
        arrow(ax, (x, 1.65), (x + 0.45, 1.22), "#d19a2e", lw=1.2, mutation_scale=9, alpha=0.8)
    ax.text(2.75, 0.98, "evaporite / sediment / tholin surface\npolymerization and retention", ha="center", va="center", fontsize=9.0, color="#f3d89b")

    # Pulse labels.
    arrow(ax, (8.55, 3.35), (7.72, 2.95), "#1b6c7c", lw=1.3, mutation_scale=11)
    ax.text(8.68, 3.45, "cryovolcanic or impact pulse\nH2O/NH3 microphase + heat\nactivates optional aq chemistry", ha="left", va="center", fontsize=9.1, color="#164e5b")

    # Phase boundary labels.
    ax.text(0.55, 4.30, "z = 0 atmosphere/pond boundary", ha="left", va="bottom", fontsize=8.8, color="#17324a")
    ax.text(0.55, 1.52, "z = L reactive bottom boundary", ha="left", va="bottom", fontsize=8.8, color="#f1d58a")
    ax.text(0.55, 3.10, "methane / ethane / N2 liquid\n1D concentration field C_i(z,t)", ha="left", va="center", fontsize=10.0, color="#dceeff")

    # Product phase annotations.
    ax.add_patch(Circle((4.55, 2.05), 0.19, facecolor="#d94675", edgecolor="white", lw=0.8, alpha=0.88))
    ax.add_patch(Circle((4.85, 2.20), 0.13, facecolor="#f0b84f", edgecolor="white", lw=0.8, alpha=0.90))
    ax.add_patch(Circle((5.05, 1.98), 0.16, facecolor="#b064c3", edgecolor="white", lw=0.8, alpha=0.88))
    ax.text(4.92, 1.66, "dense organic\nphase candidates", ha="center", va="top", fontsize=8.7, color="#f7e7ff")

    ax.text(
        5,
        5.90,
        "Titan hydrocarbon pond with transient cryovolcanic aqueous pulse",
        ha="center",
        va="top",
        fontsize=15,
        fontweight="bold",
        color="#111820",
    )
    ax.text(
        5,
        0.035,
        "Baseline solvent is hydrocarbon liquid; aqueous HCN/formamide/amide chemistry is activated only during transient geological pulses.",
        ha="center",
        va="bottom",
        fontsize=8.7,
        color="#263238",
    )

    for suffix in ("svg", "pdf", "png"):
        fig.savefig(
            OUT_DIR / f"titan_pond_cryovolcanic_pulse.{suffix}",
            bbox_inches="tight",
            dpi=300,
        )
    plt.close(fig)


if __name__ == "__main__":
    main()
