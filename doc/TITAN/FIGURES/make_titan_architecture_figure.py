from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


OUT_DIR = Path(__file__).resolve().parent


def add_box(ax, xy, width, height, title, body, facecolor, edgecolor="#263238"):
    box = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.02,rounding_size=0.035",
        linewidth=1.3,
        edgecolor=edgecolor,
        facecolor=facecolor,
    )
    ax.add_patch(box)
    x, y = xy
    ax.text(
        x + width / 2,
        y + height - 0.09,
        title,
        ha="center",
        va="top",
        fontsize=11,
        fontweight="bold",
        color="#182026",
    )
    ax.text(
        x + 0.04,
        y + height - 0.19,
        body,
        ha="left",
        va="top",
        fontsize=8.8,
        color="#263238",
        linespacing=1.24,
    )


def add_arrow(ax, start, end, label=None, rad=0.0):
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=12,
        linewidth=1.2,
        color="#37474f",
        connectionstyle=f"arc3,rad={rad}",
    )
    ax.add_patch(arrow)
    if label:
        lx = 0.5 * (start[0] + end[0])
        ly = 0.5 * (start[1] + end[1])
        ax.text(
            lx,
            ly + 0.035,
            label,
            ha="center",
            va="bottom",
            fontsize=8.2,
            color="#37474f",
        )


def main():
    fig, ax = plt.subplots(figsize=(10.8, 6.6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.text(
        0.5,
        0.965,
        "Titan atmosphere-to-protocell computational architecture",
        ha="center",
        va="top",
        fontsize=15,
        fontweight="bold",
        color="#111820",
    )
    ax.text(
        0.5,
        0.925,
        "0D atmospheric reservoir feeds a 1D methane/ethane pond with optional cryovolcanic aqueous pulses",
        ha="center",
        va="top",
        fontsize=9.5,
        color="#455a64",
    )

    add_box(
        ax,
        (0.06, 0.66),
        0.26,
        0.20,
        "0D Titan atmosphere",
        "Fixed surface T, P\nGas-phase KIDA/Hebrard network\nPhotolysis forcing J_i\nOutput: organic reservoir",
        "#d9edf7",
    )
    add_box(
        ax,
        (0.39, 0.66),
        0.25,
        0.20,
        "Feedstock boundary",
        "Rainout and deposition\nHCN, HC3N, C2Hx, C3Hx\nAerosols/tholin material\nTop boundary at z = 0",
        "#e8f5e9",
    )
    add_box(
        ax,
        (0.71, 0.66),
        0.23,
        0.20,
        "Optional pulse",
        "Cryovolcanism or impact\nH2O/NH3 microphase\nHeat and minerals\nActivates aq chemistry",
        "#fff3cd",
    )

    add_box(
        ax,
        (0.14, 0.33),
        0.32,
        0.22,
        "1D methane/ethane pond",
        "Concentrations C_i(z,t)\nEvaporation and rainfall\nDiffusion/mixing with depth\nHydrocarbon solvent is baseline",
        "#e1f0ff",
    )
    add_box(
        ax,
        (0.55, 0.33),
        0.32,
        0.22,
        "Surface and polymer chemistry",
        "Evaporites, sediments, tholins\nGeneric P_HCN_n growth\nNitrile/hydrocarbon polymers\nAq branch only during pulses",
        "#f8e8f2",
    )

    add_box(
        ax,
        (0.24, 0.08),
        0.22,
        0.15,
        "Dense phases",
        "Organic-rich droplets\nAdsorbed layers\nPolymer pools",
        "#f1f8e9",
    )
    add_box(
        ax,
        (0.55, 0.08),
        0.22,
        0.15,
        "Protocell candidates",
        "Persistence\nRetention\nGrowth/fragmentation metrics",
        "#ede7f6",
    )

    add_arrow(ax, (0.32, 0.76), (0.39, 0.76), "reservoir")
    add_arrow(ax, (0.64, 0.76), (0.71, 0.76), "geologic forcing")
    add_arrow(ax, (0.51, 0.66), (0.32, 0.55), "top input", rad=0.05)
    add_arrow(ax, (0.80, 0.66), (0.70, 0.55), "aq pulse", rad=-0.05)
    add_arrow(ax, (0.46, 0.44), (0.55, 0.44), "local reactions")
    add_arrow(ax, (0.37, 0.33), (0.35, 0.23), "phase allocation")
    add_arrow(ax, (0.67, 0.33), (0.66, 0.23), "organization")
    add_arrow(ax, (0.46, 0.155), (0.55, 0.155), "selection")

    ax.text(
        0.5,
        0.015,
        "Canonical inputs: titan_gas_phase_hebrard2013.txt + titan_photolysis_selected.txt + titan_liquid_surface_polymerization.txt",
        ha="center",
        va="bottom",
        fontsize=8.2,
        color="#455a64",
    )

    for suffix in ("svg", "pdf", "png"):
        fig.savefig(
            OUT_DIR / f"titan_atmosphere_pond_protocell_architecture.{suffix}",
            bbox_inches="tight",
            dpi=300,
        )
    plt.close(fig)


if __name__ == "__main__":
    main()
