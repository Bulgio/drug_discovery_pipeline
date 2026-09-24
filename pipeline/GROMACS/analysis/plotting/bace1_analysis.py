import matplotlib
"""
BACE1 Inhibitors - MD Simulation Analysis Suite
================================================
Matteo Bulgini - Politecnico di Torino - Tesi Magistrale 2025/2026

Usage:
    python bace1_analysis.py --input your_data.json --output ./figures

Generates:
    - Fig 1: Comparative ΔG ranking (MM-PBSA) for all 27 compounds
    - Fig 2: RMSD protein & ligand heatmap comparison
    - Fig 3: Energy decomposition (VdW, EEL, EPB, ENPOLAR) stacked bar
    - Fig 4: H-bond count and pocket contacts comparison
    - Fig 5: Pocket RMSF comparison across compounds
    - Fig 6: Stability classification summary table
    - Fig 7: Per-compound multi-panel figure (RMSF residue + SASA pocket + energy decomposition)
    - Fig 8: Correlation matrix of key metrics
    - Fig 9: Boltz score vs MM-PBSA scatter
"""

import json
import argparse
import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from scipy import stats

# ─────────────────────────────────────────────
# GLOBAL STYLE  (Nature/JMedChem look)
# ─────────────────────────────────────────────
PALETTE_MAIN   = "#178CDA"   # blue
PALETTE_ACCENT = "#CA2F1D"   # red
PALETTE_NEUTRAL= "#7F8C8D"   # grey
PALETTE_GREEN  = "#27AE60"
PALETTE_ORANGE = "#E67E22"

ENERGY_COLORS = {
    "delta_VDWAALS": "#178CDA",
    "delta_EEL":     "#8E44AD",
    "delta_EPB":     "#CA2F1D",
    "delta_ENPOLAR": "#27AE60",
}

plt.rcParams.update({
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":          9,
    "axes.titlesize":     10,
    "axes.labelsize":     9,
    "xtick.labelsize":    8,
    "ytick.labelsize":    8,
    "axes.linewidth":     0.8,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.grid":          True,
    "grid.alpha":         0.35,
    "grid.linewidth":     0.5,
    "legend.fontsize":    8,
    "legend.framealpha":  0.85,
    "figure.dpi":         150,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
    "savefig.facecolor":  "white",
})

POCKET_RESIDUES = {48, 244}   # catalytic dyad (mature protein numbering)
POCKET_RESIDUES_DISPLAY = {32: 48, 228: 244}
# ─────────────────────────────────────────────
# DATA LOADING & FLATTENING
# ─────────────────────────────────────────────

def load_data(path: str) -> list[dict]:
    with open(path) as f:
        return json.load(f)


def flatten_to_df(data: list[dict]) -> pd.DataFrame:
    """Extract scalar summary metrics into a flat DataFrame."""
    rows = []
    for c in data:
        g = c.get("gromacs", {})
        mm = g.get("mmpbsa", {}).get("summary_kJ_mol", {})
        hb = g.get("hbonds", {})
        ct = g.get("contacts", {})
        gy = g.get("gyrate", {})
        sp = g.get("sasa_protein", {})
        rp = g.get("rmsf", {}).get("rmsf_pocket_A", {})
        rl = g.get("rmsd_ligand",  {}).get("rmsd_ligand_A",  {})
        rr = g.get("rmsd_protein", {}).get("rmsd_protein_A", {})
        hm = g.get("hbmap", {}).get("hbmap", {})

        rows.append({
            "id":                  c["id"],
            "reference":           c.get("reference", ""),
            "combined_score":      c.get("combined_score", np.nan),
            "confidence_placing":  c.get("confidence_placing", np.nan),
            "affinity_placing":    c.get("affinity_placing", np.nan),
            "gromacs_rank":        g.get("rank", np.nan),
            # MM-PBSA
            "dG_total":            mm.get("delta_TOTAL",   {}).get("mean", np.nan),
            "dG_total_std":        mm.get("delta_TOTAL",   {}).get("std",  np.nan),
            "dG_vdw":              mm.get("delta_VDWAALS", {}).get("mean", np.nan),
            "dG_eel":              mm.get("delta_EEL",     {}).get("mean", np.nan),
            "dG_epb":              mm.get("delta_EPB",     {}).get("mean", np.nan),
            "dG_enpolar":          mm.get("delta_ENPOLAR", {}).get("mean", np.nan),
            "dG_ggas":             mm.get("delta_GGAS",    {}).get("mean", np.nan),
            "dG_gsolv":            mm.get("delta_GSOLV",   {}).get("mean", np.nan),
            # RMSD
            "rmsd_lig_mean":       rl.get("mean", np.nan),
            "rmsd_lig_std":        rl.get("std",  np.nan),
            "rmsd_prot_mean":      rr.get("mean", np.nan),
            "rmsd_prot_std":       rr.get("std",  np.nan),
            # H-bonds
            "hbonds_mean":         hb.get("hbonds",           {}).get("mean", np.nan),
            "hbonds_std":          hb.get("hbonds",           {}).get("std",  np.nan),
            "hbond_pairs_mean":    hb.get("hbond_pairs_035nm",{}).get("mean", np.nan),
            # Contacts
            "min_dist_mean":       ct.get("min_distance_protein_lig_nm", {}).get("mean", np.nan),
            # Gyration
            "rg_mean":             gy.get("radius_of_gyration_nm", {}).get("mean", np.nan),
            "rg_std":              gy.get("radius_of_gyration_nm", {}).get("std",  np.nan),
            # SASA
            "sasa_total_mean":     sp.get("sasa_total_nm2", {}).get("mean", np.nan),
            # RMSF pocket
            "rmsf_pocket_mean":    rp.get("mean", np.nan),
            "rmsf_pocket_std":     rp.get("std",  np.nan),
            # Hbmap
            "hbmap_max_occ":       hm.get("max_occupancy",  np.nan),
            "hbmap_mean_occ":      hm.get("mean_occupancy", np.nan),
            "hbmap_n_unique":      hm.get("n_hbonds_unique",np.nan),
        })
    df = pd.DataFrame(rows)
    df = df.sort_values("dG_total").reset_index(drop=True)
    return df


def stability_label(row) -> str:
    """Rule-based classification: stable / marginal / unstable."""
    if (row["dG_total"] < -25
            and row["rmsd_lig_mean"] < 4.0
            and row["hbonds_mean"] > 0.8):
        return "Stable"
    elif row["dG_total"] < 0:
        return "Marginal"
    else:
        return "Unstable"


# ─────────────────────────────────────────────
# FIGURE 1 — ΔG RANKING
# ─────────────────────────────────────────────

def fig_dg_ranking(df: pd.DataFrame, outdir: str):
    fig, ax = plt.subplots(figsize=(8, 5.5))

    colors = [PALETTE_ACCENT if v > -25 else PALETTE_MAIN for v in df["dG_total"]]
    bars = ax.barh(
        range(len(df)), df["dG_total"],
        xerr=df["dG_total_std"],
        color=colors, alpha=0.85,
        error_kw={"elinewidth": 0.8, "capsize": 2, "ecolor": PALETTE_NEUTRAL},
        height=0.7
    )

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels([f"{r['id']}" for _, r in df.iterrows()], fontsize=7)
    ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
    ax.axvline(-25, color=PALETTE_ORANGE, linewidth=0.8, linestyle=":", alpha=0.8,
               label="Threshold −25 kJ/mol")
    ax.set_xlabel("ΔG$_{bind}$ MM-PBSA (kJ/mol)")
    ax.set_title("Binding Free Energy — All Compounds (MM-PBSA)", fontweight="bold")
    ax.legend(loc="lower right")

    # annotate top 5
    for i, row in df.head(5).iterrows():
        ax.annotate(f"{row['dG_total']:.1f}",
                    xy=(row["dG_total"] - 0.5, i),
                    va="center", ha="right", fontsize=6.5, color="white", fontweight="bold")

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=PALETTE_MAIN,   label="ΔG < −25 kJ/mol (strong)"),
        Patch(facecolor=PALETTE_ACCENT, label="ΔG > −25 kJ/mol (weak/unfav)"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", framealpha=0.9)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig1_dg_ranking.pdf"))
    fig.savefig(os.path.join(outdir, "fig1_dg_ranking.png"))
    plt.close(fig)
    print("  ✓ Fig 1: ΔG ranking")


# ─────────────────────────────────────────────
# FIGURE 2 — RMSD HEATMAP
# ─────────────────────────────────────────────

def fig_rmsd_heatmap(df: pd.DataFrame, outdir: str):
    fig, axes = plt.subplots(1, 2, figsize=(9, 5.5))

    for ax, col, title, cmap in zip(
        axes,
        ["rmsd_lig_mean", "rmsd_prot_mean"],
        ["Ligand RMSD (Å)", "Protein RMSD (Å)"],
        ["YlOrRd", "Blues"]
    ):
        vals = df[col].values.reshape(-1, 1)
        im = ax.imshow(vals, aspect="auto", cmap=cmap,
                       vmin=df[col].min() * 0.9, vmax=df[col].max() * 1.05)
        ax.set_yticks(range(len(df)))
        ax.set_yticklabels([f"{r['id']} ({r['reference']})" for _, r in df.iterrows()], fontsize=6.5)
        ax.set_xticks([])
        ax.set_title(title, fontweight="bold")

        for i, v in enumerate(df[col].values):
            std = df[col.replace("_mean", "_std")].values[i]
            ax.text(0, i, f"{v:.1f}±{std:.1f}", ha="center", va="center",
                    fontsize=5.5, color="black" if v < df[col].mean() else "white")

        plt.colorbar(im, ax=ax, shrink=0.8, label="Å")

    fig.suptitle("RMSD Comparison — Ligand and Protein", fontweight="bold", y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig2_rmsd_heatmap.pdf"))
    fig.savefig(os.path.join(outdir, "fig2_rmsd_heatmap.png"))
    plt.close(fig)
    print("  ✓ Fig 2: RMSD heatmap")


# ─────────────────────────────────────────────
# FIGURE 3 — ENERGY DECOMPOSITION
# ─────────────────────────────────────────────

def fig_energy_decomposition(df: pd.DataFrame, outdir: str):
    top = df.head(min(15, len(df))).copy()
    ids = [f"{r['id']}\n({r['reference']})" for _, r in top.iterrows()]

    components = ["dG_vdw", "dG_eel", "dG_epb", "dG_enpolar"]
    labels     = ["ΔVdW", "ΔEEL", "ΔEPB", "ΔENpolar"]
    colors_bar = [PALETTE_MAIN, PALETTE_ACCENT, PALETTE_ORANGE, PALETTE_GREEN]

    fig, ax = plt.subplots(figsize=(10, 5))
    x    = np.arange(len(top))
    w    = 0.18
    offs = np.linspace(-(len(components)-1)/2 * w, (len(components)-1)/2 * w, len(components))

    for comp, lbl, col, off in zip(components, labels, colors_bar, offs):
        ax.bar(x + off, top[comp], width=w, label=lbl, color=col, alpha=0.85)

    ax.axhline(0, color="black", linewidth=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(ids, fontsize=6.5, rotation=45, ha="right")
    ax.set_ylabel("ΔG component (kJ/mol)")
    ax.set_title("MM-PBSA Energy Decomposition — Top Compounds", fontweight="bold")
    ax.legend(ncol=4, loc="upper right")

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig3_energy_decomposition.pdf"))
    fig.savefig(os.path.join(outdir, "fig3_energy_decomposition.png"))
    plt.close(fig)
    print("  ✓ Fig 3: Energy decomposition")


# ─────────────────────────────────────────────
# FIGURE 4 — H-BONDS & CONTACTS
# ─────────────────────────────────────────────

def fig_hbonds_contacts(df: pd.DataFrame, outdir: str):
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5))

    # Panel A: H-bond mean count
    ax = axes[0]
    colors = [PALETTE_GREEN if v > 1.0 else PALETTE_NEUTRAL for v in df["hbonds_mean"]]
    ax.barh(range(len(df)), df["hbonds_mean"], xerr=df["hbonds_std"],
            color=colors, alpha=0.85,
            error_kw={"elinewidth": 0.6, "capsize": 1.5, "ecolor": "#555"}, height=0.7)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels([f"{r['id']}" for _, r in df.iterrows()], fontsize=6.5)
    ax.axvline(1.0, color=PALETTE_ORANGE, linewidth=0.8, linestyle="--", alpha=0.8, label="≥1 H-bond")
    ax.set_xlabel("Mean H-bonds (count)")
    ax.set_title("Protein–Ligand H-bonds", fontweight="bold")
    ax.legend(fontsize=7)

    # Panel B: H-bond pairs at 3.5 Å
    ax = axes[1]
    ax.barh(range(len(df)), df["hbond_pairs_mean"],
            color=PALETTE_MAIN, alpha=0.75, height=0.7)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels([])
    ax.set_xlabel("H-bond pairs ≤ 3.5 Å (count)")
    ax.set_title("H-bond Pairs", fontweight="bold")

    # Panel C: min distance protein-ligand
    ax = axes[2]
    cmap_dist = plt.colormaps.get_cmap("RdYlGn_r")
    norm_d = plt.Normalize(df["min_dist_mean"].min(), df["min_dist_mean"].max())
    cols_d = [cmap_dist(norm_d(v)) for v in df["min_dist_mean"]]
    ax.barh(range(len(df)), df["min_dist_mean"] * 10,  # nm -> Å
            color=cols_d, alpha=0.85, height=0.7)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels([])
    ax.axvline(2.0, color=PALETTE_ACCENT, linewidth=0.8, linestyle="--",
               alpha=0.8, label="2.0 Å contact threshold")
    ax.set_xlabel("Min protein–ligand distance (Å)")
    ax.set_title("Closest Contact", fontweight="bold")
    ax.legend(fontsize=7)

    fig.suptitle("Protein–Ligand Interaction Metrics", fontweight="bold", y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig4_hbonds_contacts.pdf"))
    fig.savefig(os.path.join(outdir, "fig4_hbonds_contacts.png"))
    plt.close(fig)
    print("  ✓ Fig 4: H-bonds & contacts")


# ─────────────────────────────────────────────
# FIGURE 5 — POCKET RMSF COMPARISON
# ─────────────────────────────────────────────

def fig_pocket_rmsf(data: list[dict], df: pd.DataFrame, outdir: str):
    """Heatmap of per-residue RMSF for the active site pocket across all compounds."""
    # collect pocket residues present in at least one compound
    all_pocket_res = set()
    for c in data:
        pr = c.get("gromacs", {}).get("rmsf", {}).get("rmsf_pocket_A", {}).get("residues", {})
        all_pocket_res.update(int(k) for k in pr.keys())
    pocket_res = sorted(all_pocket_res)

    # build matrix: rows = compounds (sorted by dG), cols = pocket residues
    compound_ids = df["id"].tolist()
    matrix = np.full((len(compound_ids), len(pocket_res)), np.nan)

    id_to_compound = {c["id"]: c for c in data}
    for i, cid in enumerate(compound_ids):
        c = id_to_compound.get(cid, {})
        pr = c.get("gromacs", {}).get("rmsf", {}).get("rmsf_pocket_A", {}).get("residues", {})
        for j, res in enumerate(pocket_res):
            if str(res) in pr:
                matrix[i, j] = pr[str(res)]

    fig, ax = plt.subplots(figsize=(max(6, len(pocket_res) * 0.55), max(4, len(compound_ids) * 0.28)))

    cmap = LinearSegmentedColormap.from_list("rmsf", ["#1A5276", "#2ECC71", "#F39C12", "#C0392B"])
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=np.nanmax(matrix))

    ax.set_xticks(range(len(pocket_res)))
    xlabels = []
    actual_to_display = {v: k for k, v in POCKET_RESIDUES_DISPLAY.items()}
    for r in pocket_res:
        if r in actual_to_display:
            lbl = rf"* = Asp{str(r)}"   
        else:
            lbl = str(r)
        xlabels.append(lbl)
    ax.set_xticklabels(xlabels, rotation=90, fontsize=7)
    ax.set_yticks(range(len(compound_ids)))
    ax.set_yticklabels([f"{cid} ({df.loc[df['id']==cid,'reference'].values[0]})"
                        for cid in compound_ids], fontsize=6.5)
    ax.set_xlabel("Pocket Residue (* = catalytic dyad Asp32/Asp228, mature numbering)")
    ax.set_title("Pocket Residue RMSF Across All Compounds (Å)\n(compounds sorted by ΔG, most negative = top)",
                 fontweight="bold")

    cbar = plt.colorbar(im, ax=ax, shrink=0.6, label="RMSF (Å)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig5_pocket_rmsf_heatmap.pdf"))
    fig.savefig(os.path.join(outdir, "fig5_pocket_rmsf_heatmap.png"))
    plt.close(fig)
    print("  ✓ Fig 5: Pocket RMSF heatmap")


# ─────────────────────────────────────────────
# FIGURE 6 — STABILITY TABLE
# ─────────────────────────────────────────────

def fig_stability_table(df: pd.DataFrame, outdir: str):
    df2 = df.copy()
    df2["stability"] = df2.apply(stability_label, axis=1)

    color_map = {"Stable": PALETTE_GREEN, "Marginal": PALETTE_ORANGE, "Unstable": PALETTE_ACCENT}
    cols_show  = ["id", "reference", "dG_total", "rmsd_lig_mean", "hbonds_mean", "rmsf_pocket_mean", "stability"]
    col_labels = ["Compound ID", "DrugBank Ref", "ΔG (kJ/mol)", "RMSD lig (Å)", "H-bonds", "RMSF pocket (Å)", "Stability"]

    fig, ax = plt.subplots(figsize=(13, max(4, len(df2) * 0.32 + 1.5)))
    ax.axis("off")

    table_data = []
    for _, row in df2.iterrows():
        table_data.append([
            row["id"], row["reference"],
            f"{row['dG_total']:.2f}",
            f"{row['rmsd_lig_mean']:.2f}",
            f"{row['hbonds_mean']:.2f}",
            f"{row['rmsf_pocket_mean']:.3f}",
            row["stability"]
        ])

    tbl = ax.table(
        cellText=table_data,
        colLabels=col_labels,
        cellLoc="center", loc="center"
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)
    tbl.scale(1, 1.35)

    # header style
    for j in range(len(col_labels)):
        tbl[(0, j)].set_facecolor("#2C3E50")
        tbl[(0, j)].set_text_props(color="white", fontweight="bold")

    # row coloring by stability
    for i, row in enumerate(table_data):
        stability = row[-1]
        fc = {"Stable": "#D5F5E3", "Marginal": "#FDEBD0", "Unstable": "#FADBD8"}[stability]
        for j in range(len(col_labels)):
            tbl[(i+1, j)].set_facecolor(fc)
            if j == len(col_labels) - 1:
                tbl[(i+1, j)].set_text_props(color=color_map[stability], fontweight="bold")

    ax.set_title("Compound Stability Classification Summary",
                 fontsize=11, fontweight="bold", pad=10, y=0.98)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig6_stability_table.pdf"))
    fig.savefig(os.path.join(outdir, "fig6_stability_table.png"))
    plt.close(fig)
    print("  ✓ Fig 6: Stability table")

    # also save CSV
    df2[cols_show].to_csv(os.path.join(outdir, "stability_summary.csv"), index=False)
    print("  ✓     stability_summary.csv saved")


# ─────────────────────────────────────────────
# FIGURE 7 — PER-COMPOUND MULTI-PANEL
# ─────────────────────────────────────────────

def fig_per_compound(data: list[dict], df: pd.DataFrame, outdir: str):
    per_dir = os.path.join(outdir, "per_compound")
    os.makedirs(per_dir, exist_ok=True)

    id_to_compound = {c["id"]: c for c in data}

    for _, row in df.iterrows():
        cid = row["id"]
        c   = id_to_compound.get(cid, {})
        g   = c.get("gromacs", {})

        # --- collect data ---
        # RMSF per residue
        rmsf_all = g.get("rmsf", {}).get("rmsf_residues_A", {}).get("per_residue", {})
        pocket_res_ids = set(g.get("rmsf", {}).get("rmsf_pocket_A", {}).get("residues", {}).keys())
        res_nums  = sorted(int(k) for k in rmsf_all.keys())
        rmsf_vals = [rmsf_all[str(r)] for r in res_nums]
        pocket_mask = [str(r) in pocket_res_ids for r in res_nums]

        # SASA per residue for pocket
        sasa_pocket = g.get("sasa_residue", {}).get("sasa_pocket_nm2", {}).get("residues", {})
        sasa_res  = sorted(int(k) for k in sasa_pocket.keys())
        sasa_vals = [sasa_pocket[str(r)] for r in sasa_res]

        # Energy decomposition
        mm  = g.get("mmpbsa", {}).get("summary_kJ_mol", {})
        eng_labels = ["ΔVdW", "ΔEEL", "ΔEPB", "ΔENpolar", "ΔG total"]
        eng_keys   = ["delta_VDWAALS", "delta_EEL", "delta_EPB", "delta_ENPOLAR", "delta_TOTAL"]
        eng_vals   = [mm.get(k, {}).get("mean", 0) for k in eng_keys]
        eng_errs   = [mm.get(k, {}).get("std",  0) for k in eng_keys]
        eng_colors = [PALETTE_MAIN, PALETTE_ACCENT, PALETTE_ORANGE, PALETTE_GREEN, PALETTE_NEUTRAL]

        # --- figure ---
        fig = plt.figure(figsize=(12, 7))
        gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.38)

        # Title
        ref     = row["reference"]
        dg      = row["dG_total"]
        stab    = stability_label(row)
        stab_c  = {"Stable": PALETTE_GREEN, "Marginal": PALETTE_ORANGE, "Unstable": PALETTE_ACCENT}[stab]
        fig.suptitle(
            f"Compound {cid}  |  Reference: {ref}  |  ΔG = {dg:.2f} kJ/mol  |  "
            f"Stability: {stab}",
            fontsize=10, fontweight="bold", color="black"
        )

        # Panel A — RMSF per residue
        ax_a = fig.add_subplot(gs[0, :2])
        ax_a.plot(res_nums, rmsf_vals, color=PALETTE_NEUTRAL, linewidth=0.7, zorder=1)
        # highlight pocket residues
        pocket_x = [r for r, m in zip(res_nums, pocket_mask) if m]
        pocket_y = [rmsf_all[str(r)] for r in pocket_x]
        ax_a.scatter(pocket_x, pocket_y, color=PALETTE_ACCENT, s=22, zorder=3,
                     label="Pocket residues")
        # mark catalytic dyad
        for display_num, actual_idx in POCKET_RESIDUES_DISPLAY.items():
            if str(actual_idx) in rmsf_all:
                ax_a.axvline(actual_idx, color=PALETTE_GREEN, linewidth=1.0,
                            linestyle="--", alpha=0.7)
                ax_a.annotate(f"Asp{display_num}",
                            xy=(actual_idx, max(rmsf_vals) * 0.85),
                            fontsize=7, color=PALETTE_GREEN, ha="center")
        ax_a.axhline(3.0, color=PALETTE_ORANGE, linewidth=0.8, linestyle=":",
                     alpha=0.8, label="3.0 Å threshold")
        ax_a.set_xlabel("Residue number")
        ax_a.set_ylabel("RMSF (Å)")
        ax_a.set_title("Per-Residue RMSF", fontweight="bold")
        ax_a.legend(fontsize=7)

        # Panel B — Energy decomposition
        ax_b = fig.add_subplot(gs[0, 2])
        bars = ax_b.bar(eng_labels, eng_vals, yerr=eng_errs,
                        color=eng_colors, alpha=0.85, width=0.6,
                        error_kw={"elinewidth": 0.8, "capsize": 3, "ecolor": "#555"})
        ax_b.axhline(0, color="black", linewidth=0.7)
        ax_b.set_ylabel("kJ/mol")
        ax_b.set_title("Energy Decomposition", fontweight="bold")
        ax_b.tick_params(axis="x", rotation=30)
        # highlight ΔG total bar
        bars[-1].set_edgecolor("black")
        bars[-1].set_linewidth(1.2)

        # Panel C — SASA pocket residues
        ax_c = fig.add_subplot(gs[1, 0])
        if sasa_res:
            cmap_s = plt.colormaps.get_cmap("Blues")
            norm_s = plt.Normalize(0, max(sasa_vals) + 0.01)
            ax_c.bar([str(r) for r in sasa_res], sasa_vals,
                     color=[cmap_s(norm_s(v)) for v in sasa_vals], alpha=0.9)
            ax_c.axhline(0.20, color=PALETTE_ORANGE, linewidth=0.8, linestyle="--",
                         alpha=0.8, label="Buried threshold 0.20 nm²")
            ax_c.set_xlabel("Pocket residue")
            ax_c.set_ylabel("SASA (nm²)")
            ax_c.set_title("Pocket Residue SASA", fontweight="bold")
            ax_c.legend(fontsize=6.5)
        else:
            ax_c.text(0.5, 0.5, "No SASA data", ha="center", va="center", transform=ax_c.transAxes)
            ax_c.set_title("Pocket Residue SASA", fontweight="bold")

        # Panel D — H-bond summary text box
        ax_d = fig.add_subplot(gs[1, 1])
        ax_d.axis("off")
        hb       = g.get("hbonds", {})
        hm_data  = g.get("hbmap", {}).get("hbmap", {})
        ct       = g.get("contacts", {})
        gy_data  = g.get("gyrate", {}).get("radius_of_gyration_nm", {})
        text_lines = [
            f"H-bonds (mean):        {row['hbonds_mean']:.2f} ± {row['hbonds_std']:.2f}",
            f"H-bond pairs ≤3.5Å:    {row['hbond_pairs_mean']:.2f}",
            f"H-bonds unique:        {hm_data.get('n_hbonds_unique', '–')}",
            f"Max H-bond occupancy:  {hm_data.get('max_occupancy', 0):.3f}",
            f"Mean H-bond occupancy: {hm_data.get('mean_occupancy', 0):.3f}",
            "",
            f"Min P–L distance:      {row['min_dist_mean']*10:.2f} Å",
            f"Rg (mean):             {row['rg_mean']:.4f} ± {row['rg_std']:.4f} nm",
            f"SASA protein (mean):   {row['sasa_total_mean']:.2f} nm²",
            "",
            f"RMSD ligand (mean):    {row['rmsd_lig_mean']:.2f} ± {row['rmsd_lig_std']:.2f} Å",
            f"RMSD protein (mean):   {row['rmsd_prot_mean']:.2f} ± {row['rmsd_prot_std']:.2f} Å",
            f"RMSF pocket (mean):    {row['rmsf_pocket_mean']:.3f} ± {row['rmsf_pocket_std']:.3f} Å",
        ]
        ax_d.text(0.05, 0.95, "\n".join(text_lines),
                  transform=ax_d.transAxes,
                  va="top", ha="left", fontsize=7.5, family="monospace",
                  bbox=dict(boxstyle="round,pad=0.5", facecolor="#EBF5FB", edgecolor="#2980B9", alpha=0.8))
        ax_d.set_title("Simulation Summary", fontweight="bold")

        # Panel E — Rg text + stability badge
        ax_e = fig.add_subplot(gs[1, 2])
        ax_e.axis("off")
        badge_text = f"  {stab}  "
        ax_e.text(0.5, 0.65, badge_text, transform=ax_e.transAxes,
                  ha="center", va="center", fontsize=16, fontweight="bold",
                  color="white",
                  bbox=dict(boxstyle="round,pad=0.6", facecolor=stab_c, edgecolor="white", alpha=0.95))
        ax_e.text(0.5, 0.35,
                  f"ΔG = {dg:.2f} kJ/mol\n"
                  f"RMSD lig = {row['rmsd_lig_mean']:.2f} Å\n"
                  f"H-bonds = {row['hbonds_mean']:.2f}",
                  transform=ax_e.transAxes, ha="center", va="center",
                  fontsize=8.5, family="monospace",
                  bbox=dict(boxstyle="round,pad=0.4", facecolor="#F9F9F9", edgecolor="#CCC", alpha=0.9))
        ax_e.set_title("Classification", fontweight="bold")

        fname = os.path.join(per_dir, f"compound_{cid}")
        fig.savefig(fname + ".pdf")
        fig.savefig(fname + ".png")
        plt.close(fig)

    print(f"  ✓ Fig 7: Per-compound panels saved to {per_dir}/")


# ─────────────────────────────────────────────
# FIGURE 8 — CORRELATION MATRIX
# ─────────────────────────────────────────────

def fig_correlation_matrix(df: pd.DataFrame, outdir: str):
    cols = ["dG_total", "rmsd_lig_mean", "rmsd_prot_mean",
            "hbonds_mean", "hbond_pairs_mean", "min_dist_mean",
            "rg_mean", "sasa_total_mean", "rmsf_pocket_mean",
            "hbmap_max_occ", "combined_score"]
    col_labels = ["ΔG total", "RMSD lig", "RMSD prot",
                  "H-bonds", "H-bond pairs", "Min dist",
                  "Rg", "SASA prot", "RMSF pocket",
                  "Max Hb occ", "Boltz score"]

    sub = df[cols].dropna()
    corr = sub.corr()
    corr.columns = col_labels
    corr.index   = col_labels

    mask = np.triu(np.ones_like(corr, dtype=bool), k=1)

    fig, ax = plt.subplots(figsize=(8, 7))
    cmap = LinearSegmentedColormap.from_list("corr", ["#C0392B", "white", "#2980B9"])
    sns.heatmap(corr, mask=mask, cmap=cmap, vmin=-1, vmax=1,
                annot=True, fmt=".2f", annot_kws={"size": 6.5},
                linewidths=0.4, ax=ax, square=True,
                cbar_kws={"shrink": 0.7, "label": "Pearson r"})
    ax.set_title("Pearson Correlation Matrix — Key MD Metrics", fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig8_correlation_matrix.pdf"))
    fig.savefig(os.path.join(outdir, "fig8_correlation_matrix.png"))
    plt.close(fig)
    print("  ✓ Fig 8: Correlation matrix")


# ─────────────────────────────────────────────
# FIGURE 9 — BOLTZ vs MM-PBSA SCATTER
# ─────────────────────────────────────────────

def fig_boltz_vs_mmpbsa(df: pd.DataFrame, outdir: str):
    sub = df.dropna(subset=["dG_total", "combined_score"])
    if len(sub) < 3:
        print("  ⚠ Fig 9: not enough data for scatter, skipping")
        return

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    for ax, xcol, xlabel in zip(
        axes,
        ["affinity_placing", "combined_score"],
        ["Boltz-2 Affinity Rank ", "Boltz-2 Combined Score"]
    ):
        valid = sub.dropna(subset=[xcol])
        if len(valid) < 3:
            continue

        ax.scatter(valid[xcol], valid["dG_total"],
                   c=valid["rmsd_lig_mean"], cmap="RdYlGn_r",
                   s=55, alpha=0.85, edgecolors="grey", linewidths=0.4,
                   zorder=3)

        # regression
        slope, intercept, r, p, se = stats.linregress(valid[xcol], valid["dG_total"])
        xfit = np.linspace(valid[xcol].min(), valid[xcol].max(), 100)
        ax.plot(xfit, slope * xfit + intercept, color=PALETTE_MAIN,
                linewidth=1.2, linestyle="--", alpha=0.8,
                label=f"r = {r:.2f}  p = {p:.3f}")

        # annotate top compounds
        top5 = valid.nsmallest(5, "dG_total")
        for _, row in top5.iterrows():
            ax.annotate(row["id"], xy=(row[xcol], row["dG_total"]),
                        fontsize=5.5, ha="left", va="bottom",
                        color="#333", xytext=(3, 2), textcoords="offset points")

        sm = plt.cm.ScalarMappable(cmap="RdYlGn_r",
                                    norm=plt.Normalize(valid["rmsd_lig_mean"].min(),
                                                       valid["rmsd_lig_mean"].max()))
        sm.set_array([])
        plt.colorbar(sm, ax=ax, shrink=0.8, label="RMSD ligand (Å)")

        ax.set_xlabel(xlabel)
        ax.set_ylabel("ΔG$_{bind}$ MM-PBSA (kJ/mol)")
        ax.legend(fontsize=7)
        ax.axhline(0, color="black", linewidth=0.6, linestyle=":")

    axes[0].set_title("Boltz-2 Affinity Rank vs MM-PBSA ΔG", fontweight="bold")
    axes[1].set_title("Boltz-2 Combined Score vs MM-PBSA ΔG", fontweight="bold")

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig9_boltz_vs_mmpbsa.pdf"))
    fig.savefig(os.path.join(outdir, "fig9_boltz_vs_mmpbsa.png"))
    plt.close(fig)
    print("  ✓ Fig 9: Boltz vs MM-PBSA scatter")


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="BACE1 MD Analysis — Publication-Ready Figures"
    )
    parser.add_argument("--input",  "-i", required=True,
                        help="Path to JSON file with all 27 compound results")
    parser.add_argument("--output", "-o", default="./figures",
                        help="Output directory for figures (default: ./figures)")
    parser.add_argument("--no-per-compound", action="store_true",
                        help="Skip per-compound panels (faster)")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    print(f"\n{'─'*55}")
    print(f"  BACE1 Analysis Suite — Bulgini Thesis 2025/2026")
    print(f"{'─'*55}")
    print(f"  Input : {args.input}")
    print(f"  Output: {args.output}")
    print(f"{'─'*55}\n")

    print("Loading data...")
    data = load_data(args.input)
    df   = flatten_to_df(data)
    print(f"  Loaded {len(df)} compounds.\n")
    print(f"  ΔG range: {df['dG_total'].min():.2f} to {df['dG_total'].max():.2f} kJ/mol\n")

    print("Generating figures:")
    fig_dg_ranking(df, args.output)
    fig_rmsd_heatmap(df, args.output)
    fig_energy_decomposition(df, args.output)
    fig_hbonds_contacts(df, args.output)
    fig_pocket_rmsf(data, df, args.output)
    fig_stability_table(df, args.output)
    if not args.no_per_compound:
        fig_per_compound(data, df, args.output)
    fig_correlation_matrix(df, args.output)
    fig_boltz_vs_mmpbsa(df, args.output)

    print(f"\n{'─'*55}")
    print(f"  All figures saved to: {args.output}/")
    print(f"{'─'*55}\n")


if __name__ == "__main__":
    main()
