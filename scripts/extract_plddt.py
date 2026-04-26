#!/usr/bin/env python3
"""
extract_plddt.py
Extracts pLDDT per residue from AF2 ranked_0.pdb files.
Compares domain pLDDT alone vs in chimera context.

Usage:
    python scripts/extract_plddt.py

Output:
    results/plddt_summary.csv
    results/plddt_GFP_alone.png
    results/plddt_IHFbeta_alone.png
    results/plddt_F002_chimera.png
    results/plddt_domain_comparison.png
"""

import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── paths ─────────────────────────────────────────────────────────────────────
WORKDIR = "/rds/user/jm2564/hpc-work/chimeric_biosecurity_benchmark"

# AF2 appends the FASTA sequence name as a subdirectory
JOBS = {
    "GFP_alone": {
        "pdb": f"{WORKDIR}/af2_results/GFP_alone/GFP_alone/ranked_0.pdb",
        "json": f"{WORKDIR}/af2_results/GFP_alone/GFP_alone/ranking_debug.json",
        "length": 238,
        "type": "baseline",
    },
    "IHFbeta_alone": {
        "pdb": f"{WORKDIR}/af2_results/IHFbeta_alone/IHFbeta_alone/ranked_0.pdb",
        "json": f"{WORKDIR}/af2_results/IHFbeta_alone/IHFbeta_alone/ranking_debug.json",
        "length": 99,
        "type": "baseline",
    },
    "F002_chimera": {
        # AF2 uses the FASTA header name >F002_IHFbeta_GFP_chimera
        "pdb": f"{WORKDIR}/af2_results/F002_chimera/F002_IHFbeta_GFP_chimera/ranked_0.pdb",
        "json": f"{WORKDIR}/af2_results/F002_chimera/F002_IHFbeta_GFP_chimera/ranking_debug.json",
        "length": 343,
        "type": "chimera",
    },
}

# F002 domain boundaries (1-indexed, inclusive)
# IHFbeta: 1-99 | linker: 100-105 | GFP: 106-343
CHIMERA_DOMAINS = {
    "F002_chimera": {
        "IHFbeta": (1, 99),
        "linker":  (100, 105),
        "GFP":     (106, 343),
    }
}

OUTPUT_DIR = f"{WORKDIR}/results/plddt"
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ── helpers ───────────────────────────────────────────────────────────────────

def parse_plddt(pdb_path):
    """Extract per-residue pLDDT from CA atoms (B-factor column)."""
    plddt = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                plddt.append(float(line[60:66].strip()))
    return np.array(plddt)


def stats(arr, label=""):
    return {
        "label": label,
        "n": len(arr),
        "mean": round(float(np.mean(arr)), 2),
        "median": round(float(np.median(arr)), 2),
        "min": round(float(np.min(arr)), 2),
        "max": round(float(np.max(arr)), 2),
        "frac_gt70": round(float(np.mean(arr >= 70)), 3),
        "frac_gt90": round(float(np.mean(arr >= 90)), 3),
    }


def load_model_scores(json_path):
    """Load per-model mean pLDDT from ranking_debug.json."""
    with open(json_path) as f:
        data = json.load(f)
    return data.get("plddts", {})


# ── per-job plot ──────────────────────────────────────────────────────────────

def plot_plddt(name, plddt, domains=None, save_path=None):
    fig, ax = plt.subplots(figsize=(12, 4))
    residues = np.arange(1, len(plddt) + 1)

    # domain shading for chimeras
    colors_used = []
    if domains:
        domain_colors = {"IHFbeta": "#dbeafe", "linker": "#f3f4f6", "GFP": "#dcfce7"}
        for domain, (start, end) in domains.items():
            color = domain_colors.get(domain, "#f9fafb")
            ax.axvspan(start, end, alpha=0.3, color=color, label=f"{domain} ({start}–{end})")
            colors_used.append(domain)

    ax.plot(residues, plddt, lw=0.9, color="#1d4ed8", zorder=3)
    ax.axhline(90, color="#16a34a", ls="--", lw=0.8, label="90 — very high")
    ax.axhline(70, color="#d97706", ls="--", lw=0.8, label="70 — confident")
    ax.axhline(50, color="#dc2626", ls="--", lw=0.8, label="50 — low")

    ax.set_xlabel("Residue position", fontsize=11)
    ax.set_ylabel("pLDDT", fontsize=11)
    ax.set_title(f"{name}  (mean pLDDT = {np.mean(plddt):.2f})", fontsize=12, fontweight="bold")
    ax.set_ylim(0, 100)
    ax.set_xlim(1, len(plddt))
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"  Saved: {save_path}")
    plt.close()


# ── domain comparison plot ────────────────────────────────────────────────────

def plot_domain_comparison(results, save_path=None):
    """
    Bar chart: IHFbeta alone vs IHFbeta-in-chimera
               GFP alone     vs GFP-in-chimera
    """
    if "GFP_alone" not in results or "IHFbeta_alone" not in results or "F002_chimera" not in results:
        print("  Skipping domain comparison — not all results available yet")
        return

    plddt_f002 = results["F002_chimera"]["plddt"]
    domains = CHIMERA_DOMAINS["F002_chimera"]

    ihfbeta_in_chimera = plddt_f002[domains["IHFbeta"][0]-1 : domains["IHFbeta"][1]]
    gfp_in_chimera     = plddt_f002[domains["GFP"][0]-1    : domains["GFP"][1]]

    categories = ["IHFbeta\nalone", "IHFbeta\nin F002", "GFP\nalone", "GFP\nin F002"]
    values = [
        np.mean(results["IHFbeta_alone"]["plddt"]),
        np.mean(ihfbeta_in_chimera),
        np.mean(results["GFP_alone"]["plddt"]),
        np.mean(gfp_in_chimera),
    ]
    bar_colors = ["#3b82f6", "#1d4ed8", "#22c55e", "#15803d"]

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(categories, values, color=bar_colors, width=0.5, edgecolor="white")

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f"{val:.1f}", ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax.axhline(90, color="#16a34a", ls="--", lw=0.8, label="90 — very high confidence")
    ax.axhline(70, color="#d97706", ls="--", lw=0.8, label="70 — confident")
    ax.set_ylabel("Mean pLDDT", fontsize=11)
    ax.set_title("Domain pLDDT: alone vs in chimera F002\n(lower = more disrupted by fusion)",
                 fontsize=11, fontweight="bold")
    ax.set_ylim(0, 105)
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"  Saved: {save_path}")
    plt.close()


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    results = {}
    summary_rows = []

    print("\n" + "="*60)
    print("pLDDT Extraction — Chimeric Biosecurity Benchmark")
    print("="*60)

    for name, info in JOBS.items():
        pdb_path = info["pdb"]
        print(f"\n── {name} ──────────────────────────")

        if not os.path.exists(pdb_path):
            print(f"  MISSING: {pdb_path}")
            print(f"  Skipping — job may still be running")
            continue

        plddt = parse_plddt(pdb_path)
        s = stats(plddt, label=name)
        results[name] = {"plddt": plddt, "stats": s}

        print(f"  Residues:     {s['n']}")
        print(f"  Mean pLDDT:   {s['mean']}")
        print(f"  Median:       {s['median']}")
        print(f"  Min/Max:      {s['min']} / {s['max']}")
        print(f"  ≥70 pLDDT:    {s['frac_gt70']*100:.1f}%")
        print(f"  ≥90 pLDDT:    {s['frac_gt90']*100:.1f}%")

        # model scores
        if os.path.exists(info["json"]):
            scores = load_model_scores(info["json"])
            print(f"  Model scores: {', '.join(f'{v:.1f}' for v in scores.values())}")

        # per-job plot
        domains = CHIMERA_DOMAINS.get(name)
        plot_plddt(name, plddt, domains=domains,
                   save_path=f"{OUTPUT_DIR}/plddt_{name}.png")

        summary_rows.append(s)

    # domain comparison plot
    print(f"\n── Domain comparison ───────────────────────────────")
    plot_domain_comparison(results, save_path=f"{OUTPUT_DIR}/plddt_domain_comparison.png")

    # CSV summary
    csv_path = f"{OUTPUT_DIR}/plddt_summary.csv"
    with open(csv_path, "w") as f:
        f.write("job,n_residues,mean_pLDDT,median_pLDDT,min_pLDDT,max_pLDDT,frac_gt70,frac_gt90\n")
        for s in summary_rows:
            f.write(f"{s['label']},{s['n']},{s['mean']},{s['median']},{s['min']},{s['max']},{s['frac_gt70']},{s['frac_gt90']}\n")
    print(f"\n✓ Summary CSV: {csv_path}")

    print("\n" + "="*60)
    print("INTERPRETATION GUIDE")
    print("="*60)
    print("pLDDT ≥ 90  : very high confidence, well-structured")
    print("pLDDT 70–90 : confident, likely correct")
    print("pLDDT 50–70 : low confidence, may be disordered")
    print("pLDDT < 50  : very low, likely unstructured")
    print()
    print("For chimera benchmark:")
    print("  FAILS chimera → expect domain pLDDT DROP vs alone")
    print("  WORKS chimera → domain pLDDT should be PRESERVED")
    print("="*60)


if __name__ == "__main__":
    main()
