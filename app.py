#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ChloroCodon – Codon Usage Analyzer
-----------------------------------
A Streamlit web application for comprehensive codon usage analysis from GenBank files.

Features:
- Genome statistics (A%, T%, G%, C%, GC%, GC1%, GC2%, GC3%, A3s/T3s/C3s/G3s/GC3s, ENC, CAI, CBI, FOP)
- RSCU analysis (64 codons with their RSCU values and amino acid)
- Optimal codon analysis (high/low ENC groups, ΔRSCU), Sharp & Li CAI, Wright 1990 CBI
- ENC‑GC3 plot (with expected curve)
- PR2 bias plot (A3/(A3+T3) vs G3/(G3+C3))
- Neutrality plot (GC12 vs GC3 with regression)
- Correspondence Analysis (COA) of codon usage (PCA approximation)
- RSCU heatmap across genes with hierarchical clustering
- Correlation matrix heatmap for ENC, GC1, GC2, and GC3

All plots are publication-ready and can be saved in selected formats (PNG, PDF, SVG, TIFF).
All raw data are exported as CSV files and one Excel workbook.

Dependencies:
    pip install streamlit biopython pandas numpy matplotlib scipy scikit-learn openpyxl pillow
"""

import threading
import sys
import os
import subprocess
from datetime import datetime
from pathlib import Path
import re
import random
import io
import zipfile
import tempfile
import shutil
import base64
import html
from collections import defaultdict

import streamlit as st

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import linregress
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from Bio import SeqIO
from Bio.Data import CodonTable

try:
    from PIL import Image, ImageTk
    PIL_AVAILABLE = True
except Exception:
    Image = None
    ImageTk = None
    PIL_AVAILABLE = False


# ----------------------------------------------------------------------
#  Application theme
# ----------------------------------------------------------------------
APP_THEME = {
    "bg": "#f3f8f5",
    "bg_deep": "#071a12",
    "card": "#ffffff",
    "card_alt": "#f7fbf8",
    "glass": "rgba(255, 255, 255, 0.82)",
    "primary": "#1f7a4d",
    "primary_dark": "#0e3d28",
    "primary_light": "#dff3e8",
    "accent": "#d99a2b",
    "accent_dark": "#9a6416",
    "text": "#13251c",
    "muted": "#627267",
    "border": "#d8e8de",
    "table_even": "#ffffff",
    "table_odd": "#f4faf6",
    "table_selected": "#bfe8d1",
    "danger": "#b91c1c",
}

DEVELOPER_FOOTER = """
<div class="cc-footer-inner">
    <div class="cc-footer-brand">ChloroCodon</div>
    <div class="cc-footer-meta">GenBank parsing · Codon-bias analytics · Streamlit UI · Pandas/NumPy tables · Matplotlib figures · SciPy regression/clustering · scikit-learn PCA · Excel/ZIP exports</div>
    <div class="cc-footer-copy">© 2026 · Developed by Sheikh Sunzid Ahmed and M. Oliur Rahman</div>
    <div class="cc-footer-lab">Plant Taxonomy and Ethnobotany Laboratory, Department of Botany, University of Dhaka</div>
</div>
"""

# Fixed export DPI. The DPI control is intentionally hidden from the GUI.
EXPORT_DPI = 300

# Put your downloaded banner file beside this script using this exact name.
# The app will fall back to the solid green header if the image is missing
# or Pillow is not available.
BANNER_IMAGE_NAME = "banner_image.png"
SUMMARY_BACKGROUND_IMAGE_NAME = "background.png"
HEADER_HEIGHT = 96


# ----------------------------------------------------------------------
#  Safe plot themes
# ----------------------------------------------------------------------
# Important: The GUI does not rely on global plt.style for embedded figures.
# Reusing global styles inside Tkinter can make labels disappear, especially
# after switching to dark styles. These palettes are applied explicitly to
# every new Figure/Axes object.
PLOT_PALETTES = {
    "default": {
        "figure_bg": "#ffffff", "axes_bg": "#ffffff", "text": "#222222",
        "grid": "#d9d9d9", "marker": "#2b8c6b", "edge": "#ffffff",
        "line": "#111111", "accent": "#c95f29", "cmap": "viridis",
    },
    "seaborn-v0_8": {
        "figure_bg": "#ffffff", "axes_bg": "#f6f7fb", "text": "#202020",
        "grid": "#d0d5dd", "marker": "#4c72b0", "edge": "#ffffff",
        "line": "#1f2937", "accent": "#dd8452", "cmap": "viridis",
    },
    "ggplot": {
        "figure_bg": "#ffffff", "axes_bg": "#f0f0f0", "text": "#222222",
        "grid": "#ffffff", "marker": "#d55e00", "edge": "#ffffff",
        "line": "#111111", "accent": "#0072b2", "cmap": "plasma",
    },
    "tableau-colorblind10": {
        "figure_bg": "#ffffff", "axes_bg": "#ffffff", "text": "#222222",
        "grid": "#d8d8d8", "marker": "#117733", "edge": "#ffffff",
        "line": "#332288", "accent": "#cc6677", "cmap": "cividis",
    },
    "dark_background": {
        "figure_bg": "#151515", "axes_bg": "#1f1f1f", "text": "#f2f2f2",
        "grid": "#555555", "marker": "#70d6ff", "edge": "#111111",
        "line": "#ffffff", "accent": "#ffb703", "cmap": "viridis",
    },
    "fivethirtyeight": {
        "figure_bg": "#ffffff", "axes_bg": "#f0f0f0", "text": "#222222",
        "grid": "#c8c8c8", "marker": "#30a2da", "edge": "#ffffff",
        "line": "#fc4f30", "accent": "#e5ae38", "cmap": "viridis",
    },
    "bmh": {
        "figure_bg": "#ffffff", "axes_bg": "#eeeeee", "text": "#222222",
        "grid": "#cccccc", "marker": "#348abd", "edge": "#ffffff",
        "line": "#111111", "accent": "#a60628", "cmap": "magma",
    },
    "classic": {
        "figure_bg": "#ffffff", "axes_bg": "#ffffff", "text": "#111111",
        "grid": "#dddddd", "marker": "#1f77b4", "edge": "#ffffff",
        "line": "#000000", "accent": "#d62728", "cmap": "viridis",
    },
    "Solarize_Light2": {
        "figure_bg": "#fdf6e3", "axes_bg": "#fdf6e3", "text": "#073642",
        "grid": "#eee8d5", "marker": "#268bd2", "edge": "#fdf6e3",
        "line": "#073642", "accent": "#dc322f", "cmap": "viridis",
    },
    "seaborn-v0_8-paper": {
        "figure_bg": "#ffffff", "axes_bg": "#ffffff", "text": "#222222",
        "grid": "#d9d9d9", "marker": "#4c72b0", "edge": "#ffffff",
        "line": "#111111", "accent": "#dd8452", "cmap": "viridis",
    },
    "seaborn-v0_8-notebook": {
        "figure_bg": "#ffffff", "axes_bg": "#ffffff", "text": "#222222",
        "grid": "#d9d9d9", "marker": "#55a868", "edge": "#ffffff",
        "line": "#111111", "accent": "#c44e52", "cmap": "viridis",
    },
    "seaborn-v0_8-talk": {
        "figure_bg": "#ffffff", "axes_bg": "#ffffff", "text": "#222222",
        "grid": "#d9d9d9", "marker": "#8172b3", "edge": "#ffffff",
        "line": "#111111", "accent": "#ccb974", "cmap": "plasma",
    },
    "seaborn-v0_8-poster": {
        "figure_bg": "#ffffff", "axes_bg": "#ffffff", "text": "#222222",
        "grid": "#d9d9d9", "marker": "#64b5cd", "edge": "#ffffff",
        "line": "#111111", "accent": "#dd8452", "cmap": "viridis",
    },
}


# ----------------------------------------------------------------------
#  Core analysis functions
# ----------------------------------------------------------------------

# Standard DNA codon order (64 codons)
CODON_ORDER = [
    # 64 DNA codons, including stop codons.
    # Stop codons are displayed in the RSCU tab with amino acid symbol "*".
    "TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
    "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG",
    "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG",
    "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG",
    "TAT", "TAC", "TAA", "TAG", "TGA", "TGT", "TGC", "TGG",
    "CAT", "CAC", "CAA", "CAG",
    "AAT", "AAC", "AAA", "AAG",
    "GAT", "GAC", "GAA", "GAG",
    "CGT", "CGC", "CGA", "CGG", "AGT", "AGC", "AGA", "AGG",
    "GGT", "GGC", "GGA", "GGG"
]
DNA_BASES = set("ACGT")
AA_TO_CODONS = None
STOP_CODONS = set()


def get_codon_table(table_id=11):
    """Return an NCBI unambiguous DNA codon table."""
    try:
        return CodonTable.unambiguous_dna_by_id[int(table_id)]
    except Exception as exc:
        raise ValueError(f"Unsupported or unknown NCBI genetic code table: {table_id}") from exc


def init_codon_maps(table_id=11):
    """Initialize global AA_TO_CODONS and STOP_CODONS."""
    global AA_TO_CODONS, STOP_CODONS
    table = get_codon_table(table_id)
    aa_map = defaultdict(list)
    for codon, aa in table.forward_table.items():
        aa_map[aa].append(codon)
    AA_TO_CODONS = {aa: sorted(codons, key=lambda c: CODON_ORDER.index(c)) for aa, codons in aa_map.items()}
    STOP_CODONS = set(table.stop_codons)


def parse_genbank(file_path, min_codons=30, table_id=11, trim_partial_tail=True):
    """Parse GenBank file, extract and filter CDS features."""
    records = []
    qc_rows = []

    def qualifier_first(q, key, default=""):
        values = q.get(key)
        if not values:
            return default
        return str(values[0])

    def is_pseudo(q):
        return "pseudo" in q or "pseudogene" in q

    with open(file_path, "r", encoding="utf-8", errors="replace") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            organism = record.annotations.get("organism", "")
            for idx, feature in enumerate(record.features, start=1):
                if feature.type != "CDS":
                    continue
                q = feature.qualifiers
                gene = qualifier_first(q, "gene")
                locus_tag = qualifier_first(q, "locus_tag")
                product = qualifier_first(q, "product")
                protein_id = qualifier_first(q, "protein_id")
                label = gene or locus_tag or protein_id or product or f"CDS_{idx}"

                if is_pseudo(q):
                    qc_rows.append({"label": label, "status": "skipped", "reason": "pseudo/pseudogene"})
                    continue

                try:
                    transl_table = int(qualifier_first(q, "transl_table", str(table_id)))
                except ValueError:
                    transl_table = table_id

                try:
                    seq = str(feature.extract(record.seq)).upper().replace("U", "T")
                except Exception as exc:
                    qc_rows.append({"label": label, "status": "skipped", "reason": f"extract failed: {exc}"})
                    continue

                try:
                    codon_start = int(qualifier_first(q, "codon_start", "1"))
                except ValueError:
                    codon_start = 1
                codon_start = min(max(codon_start, 1), 3)
                if codon_start > 1:
                    seq = seq[codon_start - 1:]

                if not set(seq).issubset(DNA_BASES):
                    qc_rows.append({"label": label, "status": "skipped", "reason": "ambiguous/non-ACGT bases"})
                    continue

                if len(seq) < 3:
                    qc_rows.append({"label": label, "status": "skipped", "reason": "too short"})
                    continue

                remainder = len(seq) % 3
                if remainder:
                    if trim_partial_tail:
                        seq = seq[:len(seq) - remainder]
                    else:
                        qc_rows.append({"label": label, "status": "skipped", "reason": "length not divisible by 3"})
                        continue

                try:
                    table = get_codon_table(transl_table)
                except ValueError as exc:
                    qc_rows.append({"label": label, "status": "skipped", "reason": str(exc)})
                    continue

                codons_all = [seq[i:i+3] for i in range(0, len(seq), 3)]

                # Keep a copy with the terminal stop codon for the RSCU table.
                # For ENC, PR2, neutrality, COA, and genome statistics, terminal
                # stop codons are excluded because they do not encode amino acids.
                terminal_stop_codon = ""
                sequence_with_stop = seq

                codons = codons_all[:]
                if codons and codons[-1] in table.stop_codons:
                    terminal_stop_codon = codons[-1]
                    codons = codons[:-1]
                    seq = "".join(codons)
                else:
                    sequence_with_stop = seq

                # After removing a possible terminal stop, any remaining stop codon
                # is internal and the CDS is skipped.
                if any(c in table.stop_codons for c in codons):
                    qc_rows.append({"label": label, "status": "skipped", "reason": "internal stop codon"})
                    continue

                if len(codons) < min_codons:
                    qc_rows.append({"label": label, "status": "skipped", "reason": f"short (< {min_codons} codons)"})
                    continue

                records.append({
                    "label": label,
                    "gene": gene,
                    "locus_tag": locus_tag,
                    "product": product,
                    "protein_id": protein_id,
                    "organism": organism,
                    "record_id": record.id,
                    "feature_index": idx,
                    "location": str(feature.location),
                    "strand": "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else "",
                    "transl_table": transl_table,
                    "sequence": seq,
                    "sequence_with_stop": sequence_with_stop,
                    "terminal_stop_codon": terminal_stop_codon,
                    "length_codons": len(codons)
                })
                qc_rows.append({"label": label, "status": "accepted", "reason": ""})

    # Remove duplicate genes (keep longest) and report duplicated IR copies.
    grouped = defaultdict(list)
    for rec in records:
        key = (rec["gene"] or rec["locus_tag"] or rec["protein_id"] or rec["label"])
        grouped[key].append(rec)

    unique_records = []
    duplicate_rows = []
    for key, copies in grouped.items():
        kept = max(copies, key=lambda r: r["length_codons"])
        unique_records.append(kept)
        removed = [r for r in copies if r is not kept]
        if len(copies) > 1:
            duplicate_rows.append({
                "gene_key": key,
                "copies_found": len(copies),
                "kept_label": kept.get("label", ""),
                "kept_gene": kept.get("gene", ""),
                "kept_length_codons": kept.get("length_codons", ""),
                "kept_location": kept.get("location", ""),
                "removed_labels": "; ".join(r.get("label", "") for r in removed),
                "removed_lengths_codons": "; ".join(str(r.get("length_codons", "")) for r in removed),
                "removed_locations": "; ".join(r.get("location", "") for r in removed),
                "reason": "duplicate gene/CDS copy removed; longest CDS retained"
            })

    qc_df = pd.DataFrame(qc_rows)
    duplicate_df = pd.DataFrame(duplicate_rows)
    return unique_records, qc_df, duplicate_df


def compute_codon_counts(seq, table_id=11):
    """Return a dict of codon counts for all 64 codons (sense + stop)."""
    table = get_codon_table(table_id)
    counts = {codon: 0 for codon in CODON_ORDER}
    valid_codons = set(table.forward_table.keys()) | set(table.stop_codons)
    for i in range(0, len(seq)-2, 3):
        c = seq[i:i+3]
        if c in valid_codons:
            counts[c] += 1
    return counts


def gc_content(seq):
    """Overall GC fraction."""
    if not seq:
        return np.nan
    return (seq.count("G") + seq.count("C")) / len(seq)


def positional_gc(seq):
    """GC1, GC2, GC3 fractions."""
    if len(seq) < 3:
        return np.nan, np.nan, np.nan
    pos1 = seq[0::3]
    pos2 = seq[1::3]
    pos3 = seq[2::3]
    def gc_frac(s):
        return (s.count("G") + s.count("C")) / len(s) if len(s) else np.nan
    return gc_frac(pos1), gc_frac(pos2), gc_frac(pos3)


def effective_number_of_codons(codon_counts, table_id=11):
    """ENC using Wright's homozygosity approach."""
    global AA_TO_CODONS
    if AA_TO_CODONS is None:
        init_codon_maps(table_id)
    family_f = {2: [], 3: [], 4: [], 6: []}
    for aa, codons in AA_TO_CODONS.items():
        k = len(codons)
        if k not in family_f:
            continue
        values = np.array([float(codon_counts.get(c, 0)) for c in codons], dtype=float)
        n = values.sum()
        if n <= 1:
            continue
        p = values / n
        s = float(np.sum(p**2))
        f = (n * s - 1.0) / (n - 1.0)
        if np.isfinite(f) and f > 0:
            family_f[k].append(f)
    avg_f = {k: (float(np.mean(v)) if v else 1.0) for k, v in family_f.items()}
    f2 = avg_f.get(2, 1.0)
    f3 = avg_f.get(3, 1.0)
    f4 = avg_f.get(4, 1.0)
    f6 = avg_f.get(6, 1.0)
    enc = 2.0 + 9.0/f2 + 1.0/f3 + 5.0/f4 + 3.0/f6
    if not np.isfinite(enc):
        return np.nan
    return min(max(enc, 20.0), 61.0)


def expected_enc(gc3):
    """Expected ENC under the GC3-only model."""
    if pd.isna(gc3):
        return np.nan
    s = float(gc3)
    denom = s**2 + (1.0 - s)**2
    if denom == 0:
        return np.nan
    return 2.0 + s + 29.0 / denom


def rscu_from_counts(codon_counts, table_id=11):
    """
    Compute RSCU for all 64 codons.
    Sense codons grouped by amino acid; stop codons grouped together.
    Returns (rscu_dict, codon_to_aa_dict).
    """
    global AA_TO_CODONS, STOP_CODONS
    if AA_TO_CODONS is None or not STOP_CODONS:
        init_codon_maps(table_id)

    rscu = {c: 0.0 for c in CODON_ORDER}
    codon_to_aa = {}

    # Sense codons
    for aa, codons in AA_TO_CODONS.items():
        total = sum(codon_counts.get(c, 0) for c in codons)
        n = len(codons)
        if total == 0:
            for c in codons:
                rscu[c] = 0.0
                codon_to_aa[c] = aa
            continue
        for c in codons:
            rscu[c] = (codon_counts.get(c, 0) * n) / total
            codon_to_aa[c] = aa

    # Stop codons
    stop_codons = [c for c in CODON_ORDER if c in STOP_CODONS]
    total_stop = sum(codon_counts.get(c, 0) for c in stop_codons)
    n_stop = len(stop_codons)
    if total_stop > 0:
        for c in stop_codons:
            rscu[c] = (codon_counts.get(c, 0) * n_stop) / total_stop
            codon_to_aa[c] = "*"
    else:
        for c in stop_codons:
            rscu[c] = 0.0
            codon_to_aa[c] = "*"

    return rscu, codon_to_aa


def aggregate_rscu(records, table_id=11):
    """
    Compute overall RSCU from all genes.

    For the RSCU publication table, terminal stop codons are counted when they
    are present in the GenBank CDS feature. Stop codons are grouped together
    and displayed with amino acid symbol "*".
    """
    total_counts = defaultdict(int)
    for rec in records:
        seq_for_rscu = rec.get("sequence_with_stop", rec["sequence"])
        counts = compute_codon_counts(seq_for_rscu, table_id)
        for c, v in counts.items():
            total_counts[c] += v
    return rscu_from_counts(total_counts, table_id)


def third_position_stats(seq):
    """A3, T3, G3, C3 counts and biases."""
    pos3 = seq[2::3]
    if not pos3:
        return {"A3":0, "T3":0, "G3":0, "C3":0, "AT_bias":np.nan, "GC_bias":np.nan}
    A = pos3.count("A")
    T = pos3.count("T")
    G = pos3.count("G")
    C = pos3.count("C")
    at = A + T
    gc = G + C
    return {
        "A3": A, "T3": T, "G3": G, "C3": C,
        "AT_bias": A/at if at else np.nan,
        "GC_bias": G/gc if gc else np.nan
    }


def compute_per_gene_metrics(rec, table_id=11):
    """Compute all metrics for one CDS."""
    # seq excludes terminal stop codon and is used for biological metrics.
    seq = rec["sequence"]

    # seq_for_rscu keeps a terminal stop codon if the GenBank CDS location has it.
    # This allows TAA/TAG/TGA to appear in the RSCU tab as "*" codons.
    seq_for_rscu = rec.get("sequence_with_stop", seq)

    codon_counts = compute_codon_counts(seq, table_id)
    codon_counts_for_rscu = compute_codon_counts(seq_for_rscu, table_id)

    # total_codons excludes terminal stop codons.
    total_codons = sum(codon_counts.values())

    gc = gc_content(seq)
    gc1, gc2, gc3 = positional_gc(seq)
    gc12 = (gc1 + gc2) / 2 if not np.isnan(gc1) and not np.isnan(gc2) else np.nan
    enc = effective_number_of_codons(codon_counts, table_id)
    enc_exp = expected_enc(gc3)
    third = third_position_stats(seq)
    rscu, _ = rscu_from_counts(codon_counts_for_rscu, table_id)
    return {
        "label": rec["label"],
        "gene": rec["gene"],
        "locus_tag": rec["locus_tag"],
        "protein_id": rec["protein_id"],
        "terminal_stop_codon": rec.get("terminal_stop_codon", ""),
        "organism": rec["organism"],
        "length_codons": rec["length_codons"],
        "total_codons": total_codons,
        "GC": gc,
        "GC1": gc1,
        "GC2": gc2,
        "GC3": gc3,
        "GC12": gc12,
        "ENC": enc,
        "ENC_exp": enc_exp,
        "A3": third["A3"],
        "T3": third["T3"],
        "G3": third["G3"],
        "C3": third["C3"],
        "AT_bias": third["AT_bias"],
        "GC_bias": third["GC_bias"],
        "codon_counts": codon_counts,
        "rscu": rscu
    }



def aggregate_genome_stats(metrics_list):
    """
    Compute enriched overall genome/codon-composition statistics.

    Notes
    -----
    - A%, T%, C%, G%, GC%, GC1%, GC2%, and GC3% are calculated from all
      accepted non-stop CDS codons after duplicate-gene removal.
    - A3s, T3s, C3s, G3s, and GC3s are calculated from synonymous codon
      families only, excluding single-codon amino acids and stop codons.
    - CAI and CBI are added later after optimal-codon and reference-set
      calculations are completed.
    """
    if not metrics_list:
        return {}

    if AA_TO_CODONS is None or not STOP_CODONS:
        init_codon_maps(11)

    total_codons = sum(m["total_codons"] for m in metrics_list)
    total_counts = defaultdict(int)
    for m in metrics_list:
        for c, v in m["codon_counts"].items():
            total_counts[c] += v

    base_totals = {"A": 0, "T": 0, "G": 0, "C": 0}
    pos1 = {"A": 0, "T": 0, "G": 0, "C": 0}
    pos2 = {"A": 0, "T": 0, "G": 0, "C": 0}
    pos3 = {"A": 0, "T": 0, "G": 0, "C": 0}

    for codon, v in total_counts.items():
        if codon in STOP_CODONS:
            continue
        if len(codon) != 3:
            continue
        base_totals[codon[0]] += v
        base_totals[codon[1]] += v
        base_totals[codon[2]] += v
        pos1[codon[0]] += v
        pos2[codon[1]] += v
        pos3[codon[2]] += v

    total_bases = sum(base_totals.values())
    if total_bases == 0 or total_codons == 0:
        return {}

    def pct(num, den):
        return (num / den * 100.0) if den else np.nan

    syn3 = {"A": 0, "T": 0, "G": 0, "C": 0}
    for aa, codons in AA_TO_CODONS.items():
        if len(codons) <= 1:
            continue
        for codon in codons:
            if codon in STOP_CODONS:
                continue
            syn3[codon[2]] += total_counts.get(codon, 0)
    syn3_total = sum(syn3.values())

    enc_values = [m["ENC"] for m in metrics_list if not np.isnan(m["ENC"])]

    return {
        "Total_bp": total_bases,
        "Total_codons": total_codons,
        "Num_genes": len(metrics_list),

        "A%": pct(base_totals["A"], total_bases),
        "T%": pct(base_totals["T"], total_bases),
        "C%": pct(base_totals["C"], total_bases),
        "G%": pct(base_totals["G"], total_bases),
        "GC%": pct(base_totals["G"] + base_totals["C"], total_bases),

        "GC1%": pct(pos1["G"] + pos1["C"], sum(pos1.values())),
        "GC2%": pct(pos2["G"] + pos2["C"], sum(pos2.values())),
        "GC3%": pct(pos3["G"] + pos3["C"], sum(pos3.values())),

        "A3s": pct(syn3["A"], syn3_total),
        "T3s": pct(syn3["T"], syn3_total),
        "C3s": pct(syn3["C"], syn3_total),
        "G3s": pct(syn3["G"], syn3_total),
        "GC3s": pct(syn3["G"] + syn3["C"], syn3_total),

        "Avg_ENC": float(np.mean(enc_values)) if enc_values else np.nan,
    }

def optimal_codon_analysis(metrics_list, top_frac=0.1, delta_threshold=0.08):
    """Identify optimal and high-frequency codons."""
    enc_sorted = sorted(metrics_list, key=lambda m: m["ENC"] if not np.isnan(m["ENC"]) else 999)
    n = len(enc_sorted)
    if n < 2:
        return [], [], pd.Series()
    n_top = max(1, int(n * top_frac))
    n_low = max(1, int(n * top_frac))
    high_group = enc_sorted[:n_top]
    low_group = enc_sorted[-n_low:]

    def avg_rscu(group):
        avg = {c: 0.0 for c in CODON_ORDER}
        n = len(group)
        if n == 0:
            return avg
        for m in group:
            for c, v in m["rscu"].items():
                avg[c] += v / n
        return avg

    high_avg = avg_rscu(high_group)
    low_avg = avg_rscu(low_group)

    delta = {c: high_avg[c] - low_avg[c] for c in CODON_ORDER}

    # Optimal codon and HF codon lists should use sense codons only.
    # Stop codons are displayed in the RSCU tab, but are not treated as
    # optimal amino-acid codons.
    sense_codons = [c for c in CODON_ORDER if c not in STOP_CODONS]
    opt = [c for c in sense_codons if delta.get(c, 0.0) >= delta_threshold]
    hf = [c for c in sense_codons if high_avg.get(c, 0.0) > 1.0]
    return opt, hf, pd.Series(delta)


def _high_expression_reference_group(metrics_list, top_frac=0.10):
    """Return the low-ENC gene set used as a high-expression proxy."""
    valid = [m for m in metrics_list if not np.isnan(m.get("ENC", np.nan))]
    if not valid:
        return []
    enc_sorted = sorted(valid, key=lambda m: m["ENC"])
    n_ref = max(1, int(len(enc_sorted) * top_frac))
    return enc_sorted[:n_ref]


def build_cai_weights_sharp_li(metrics_list, top_frac=0.10, table_id=11, zero_weight=1e-6):
    """
    Build Sharp & Li CAI relative-adaptiveness weights from the low-ENC
    reference set.

    w_ij = f_ij / max(f_i)

    Stop codons and one-codon amino acids are excluded.
    """
    if AA_TO_CODONS is None or not STOP_CODONS:
        init_codon_maps(table_id)

    ref_group = _high_expression_reference_group(metrics_list, top_frac=top_frac)
    ref_counts = defaultdict(int)
    for m in ref_group:
        for codon, value in m["codon_counts"].items():
            ref_counts[codon] += int(value)

    weights = {}
    for aa, codons in AA_TO_CODONS.items():
        codons = [c for c in codons if c not in STOP_CODONS]
        if len(codons) <= 1:
            continue

        values = np.array([ref_counts.get(c, 0) for c in codons], dtype=float)
        max_value = float(values.max()) if len(values) else 0.0
        if max_value <= 0:
            for c in codons:
                weights[c] = 1.0
            continue

        for c, v in zip(codons, values):
            w = float(v / max_value)
            weights[c] = max(w, zero_weight)

    return weights


def calculate_cai_sharp_li_from_counts(codon_counts, cai_weights):
    """
    Calculate Sharp & Li CAI as the geometric mean of relative adaptiveness.

    CAI = exp(mean(log(w_i)))
    """
    log_weights = []
    for codon, count in codon_counts.items():
        if codon not in cai_weights:
            continue
        n = int(count)
        if n <= 0:
            continue
        w = max(float(cai_weights.get(codon, 1.0)), 1e-12)
        log_weights.extend([np.log(w)] * n)

    if not log_weights:
        return np.nan

    return float(np.exp(np.mean(log_weights)))


def calculate_cai_sharp_li(metrics_list, cai_weights):
    """Genome-level Sharp & Li CAI from aggregated non-stop CDS codon counts."""
    total_counts = defaultdict(int)
    for m in metrics_list:
        for c, v in m["codon_counts"].items():
            total_counts[c] += int(v)
    return calculate_cai_sharp_li_from_counts(total_counts, cai_weights)


def calculate_cbi_wright_from_counts(codon_counts, optimal_codons, table_id=11):
    """
    Calculate Wright's Codon Bias Index (CBI).

    CBI = (N_opt - N_exp) / (N - N_exp)
    """
    if AA_TO_CODONS is None or not STOP_CODONS:
        init_codon_maps(table_id)

    opt_set = set(optimal_codons or [])
    if not opt_set:
        return 0.0

    n_total = 0.0
    n_opt = 0.0
    n_exp = 0.0

    for aa, codons in AA_TO_CODONS.items():
        codons = [c for c in codons if c not in STOP_CODONS]
        if len(codons) <= 1:
            continue

        opt_in_family = [c for c in codons if c in opt_set]
        if not opt_in_family:
            continue

        family_total = float(sum(codon_counts.get(c, 0) for c in codons))
        if family_total <= 0:
            continue

        n_total += family_total
        n_opt += float(sum(codon_counts.get(c, 0) for c in opt_in_family))
        n_exp += family_total * (len(opt_in_family) / len(codons))

    denom = n_total - n_exp
    if denom == 0:
        return np.nan

    return float((n_opt - n_exp) / denom)


def calculate_cbi_wright(metrics_list, optimal_codons, table_id=11):
    """Genome-level Wright 1990 CBI from aggregated codon counts."""
    total_counts = defaultdict(int)
    for m in metrics_list:
        for c, v in m["codon_counts"].items():
            total_counts[c] += int(v)
    return calculate_cbi_wright_from_counts(total_counts, optimal_codons, table_id=table_id)


def calculate_fop_from_counts(codon_counts, optimal_codons, table_id=11):
    """Frequency of optimal codons over families where optimal codons are defined."""
    if AA_TO_CODONS is None or not STOP_CODONS:
        init_codon_maps(table_id)

    opt_set = set(optimal_codons or [])
    if not opt_set:
        return 0.0

    n_total = 0
    n_opt = 0
    for aa, codons in AA_TO_CODONS.items():
        codons = [c for c in codons if c not in STOP_CODONS]
        if len(codons) <= 1:
            continue

        opt_in_family = [c for c in codons if c in opt_set]
        if not opt_in_family:
            continue

        n_total += sum(codon_counts.get(c, 0) for c in codons)
        n_opt += sum(codon_counts.get(c, 0) for c in opt_in_family)

    return float(n_opt / n_total) if n_total else 0.0


def calculate_fop(metrics_list, optimal_codons, table_id=11):
    """Genome-level FOP."""
    total_counts = defaultdict(int)
    for m in metrics_list:
        for c, v in m["codon_counts"].items():
            total_counts[c] += int(v)
    return calculate_fop_from_counts(total_counts, optimal_codons, table_id=table_id)



def aggregate_codon_counts(records, table_id=11, include_terminal_stop=True):
    """Aggregate codon counts across accepted, duplicate-filtered CDS records."""
    total_counts = {codon: 0 for codon in CODON_ORDER}
    for rec in records:
        seq = rec.get("sequence_with_stop", rec["sequence"]) if include_terminal_stop else rec["sequence"]
        counts = compute_codon_counts(seq, table_id=table_id)
        for codon in CODON_ORDER:
            total_counts[codon] += int(counts.get(codon, 0))
    return total_counts


def build_rscu_publication_table(records, table_id=11):
    """Build the publication-style RSCU table: Codon, amino acid, count, RFSC, RSCU, preferred."""
    total_counts = aggregate_codon_counts(records, table_id=table_id, include_terminal_stop=True)
    rscu_dict, codon_to_aa = rscu_from_counts(total_counts, table_id=table_id)
    total_codons = sum(total_counts.values())
    rows = []
    for codon in CODON_ORDER:
        aa = codon_to_aa.get(codon, "")
        count = int(total_counts.get(codon, 0))
        rfsc_value = (count / total_codons) if total_codons else 0.0
        rscu_value = float(rscu_dict.get(codon, 0.0))
        preferred = "Yes" if aa != "*" and rscu_value > 1.0 else "No"
        rows.append({
            "Codon": codon,
            "Amino_acid": aa,
            "Count": count,
            "RFSC": rfsc_value,
            "RSCU": rscu_value,
            "Preferred_RSCU_gt_1": preferred,
        })
    return pd.DataFrame(rows)


def build_amino_acid_usage_table(records, table_id=11):
    """Build amino-acid usage table from accepted CDS, including a stop summary row."""
    if AA_TO_CODONS is None or not STOP_CODONS:
        init_codon_maps(table_id)
    counts = aggregate_codon_counts(records, table_id=table_id, include_terminal_stop=True)
    rows = []
    total_sense = 0
    for aa, codons in sorted(AA_TO_CODONS.items()):
        n = int(sum(counts.get(c, 0) for c in codons))
        total_sense += n
        rows.append({
            "Amino_acid": aa,
            "Codons": ", ".join(codons),
            "Count": n,
            "Frequency_percent": np.nan,
        })
    stop_codons = [c for c in CODON_ORDER if c in STOP_CODONS]
    stop_count = int(sum(counts.get(c, 0) for c in stop_codons))
    rows.append({"Amino_acid": "*", "Codons": ", ".join(stop_codons), "Count": stop_count, "Frequency_percent": np.nan})
    denom = total_sense if total_sense else 0
    for row in rows:
        if row["Amino_acid"] != "*":
            row["Frequency_percent"] = (row["Count"] / denom * 100.0) if denom else 0.0
        else:
            row["Frequency_percent"] = np.nan
    return pd.DataFrame(rows)


def build_stop_codon_usage_table(records, table_id=11):
    """Build stop-codon usage summary from terminal stop codons retained for RSCU."""
    if AA_TO_CODONS is None or not STOP_CODONS:
        init_codon_maps(table_id)
    stop_codons = [c for c in CODON_ORDER if c in STOP_CODONS]
    counts = {c: 0 for c in stop_codons}
    for rec in records:
        stop = rec.get("terminal_stop_codon", "")
        if stop in counts:
            counts[stop] += 1
    total = sum(counts.values())
    most_used = max(counts, key=counts.get) if total else ""
    rows = []
    for codon in stop_codons:
        rows.append({
            "Stop_codon": codon,
            "Amino_acid_symbol": "*",
            "Count": counts[codon],
            "Percent": (counts[codon] / total * 100.0) if total else 0.0,
            "Most_used": "Yes" if codon == most_used and total else "No",
        })
    return pd.DataFrame(rows)


def build_per_gene_rscu_matrix(metrics_list):
    """Return a gene-by-codon RSCU matrix for downstream COA/heatmap work."""
    rows = []
    for m in metrics_list:
        row = {"label": m.get("label", ""), "gene": m.get("gene", "")}
        for codon in CODON_ORDER:
            row[codon] = m.get("rscu", {}).get(codon, 0.0)
        rows.append(row)
    return pd.DataFrame(rows)




def build_rscu_heatmap_matrix(metrics_list, table_id=11):
    """
    Build a genes × all-codons RSCU matrix for heatmap analysis.

    All 64 codons are included. Stop codons are retained from the
    per-gene RSCU dictionaries and are displayed in the heatmap labels
    with an asterisk (*) symbol.
    """
    if AA_TO_CODONS is None or not STOP_CODONS:
        init_codon_maps(table_id)

    heatmap_codons = list(CODON_ORDER)
    rows = []
    used_labels = defaultdict(int)

    for m in metrics_list:
        base_label = str(m.get("gene") or m.get("label") or "gene")
        used_labels[base_label] += 1
        label = base_label if used_labels[base_label] == 1 else f"{base_label}_{used_labels[base_label]}"
        row = {"label": label, "gene": m.get("gene", "")}
        for codon in heatmap_codons:
            row[codon] = float(m.get("rscu", {}).get(codon, 0.0))
        rows.append(row)

    if not rows:
        return pd.DataFrame(columns=["label", "gene"] + heatmap_codons)
    return pd.DataFrame(rows)


def build_correlation_metrics_table(metrics_list):
    """Build the per-gene table used for codon-bias correlation heatmap."""
    rows = []
    for m in metrics_list:
        rows.append({
            "label": m.get("label", ""),
            "gene": m.get("gene", ""),
            "ENC": m.get("ENC", np.nan),
            "GC1": m.get("GC1", np.nan),
            "GC2": m.get("GC2", np.nan),
            "GC3": m.get("GC3", np.nan),
        })
    return pd.DataFrame(rows)


def build_bias_correlation_matrix(correlation_metrics_df, method="pearson"):
    """Return Pearson/Spearman correlation matrix for ENC, GC1, GC2, and GC3."""
    metric_cols = ["ENC", "GC1", "GC2", "GC3"]
    if correlation_metrics_df is None or correlation_metrics_df.empty:
        return pd.DataFrame(index=metric_cols, columns=metric_cols, dtype=float)
    df = correlation_metrics_df.copy()
    for col in metric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df[metric_cols].corr(method=method)

def generate_methods_text(data, min_codons=30, file_path=""):
    """Generate reusable methods text for paper/thesis reporting."""
    stats = data.get("genome_stats", {}) if data else {}
    n_genes = stats.get("Num_genes", "")
    qc_df = data.get("qc_df") if data else None
    accepted = int((qc_df["status"] == "accepted").sum()) if qc_df is not None and not qc_df.empty else ""
    skipped = int(len(qc_df) - accepted) if isinstance(accepted, int) and qc_df is not None else ""
    return f"""ChloroCodon analysis methods
Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Input file: {Path(file_path).name if file_path else ''}

Coding sequences were extracted from GenBank CDS features. Pseudogene/pseudo features, ambiguous CDS sequences, CDS with internal stop codons, and CDS shorter than {min_codons} codons were excluded. Terminal stop codons were retained only for stop-codon/RSCU reporting and excluded from ENC, PR2, neutrality, COA, and coding-composition statistics. Duplicate gene copies were collapsed by retaining the longest CDS copy.

Accepted CDS before duplicate removal: {accepted}
Skipped CDS: {skipped}
CDS used after duplicate removal: {n_genes}

RSCU was calculated by normalizing each codon count by the expected count under equal synonymous codon usage within the same amino-acid family. Stop codons were grouped together and displayed with amino-acid symbol '*'. ENC was estimated using Wright's homozygosity approach and plotted against GC3 with the expected ENC curve. PR2 was plotted as A3/(A3+T3) versus G3/(G3+C3). Neutrality analysis was performed by regressing GC12 against GC3. COA was approximated by PCA of standardized codon-frequency profiles. Optimal codons were identified using low-ENC CDS as the high-expression proxy group and high-ENC CDS as the low-expression proxy group; codons with ΔRSCU ≥ 0.08 were reported as optimal. CAI was calculated using the Sharp and Li relative-adaptiveness approach with the low-ENC reference set. CBI was calculated using Wright's codon bias index formula. FOP was calculated as the frequency of detected optimal codons within synonymous families where optimal codons were defined. A gene-by-codon RSCU matrix including all 64 codons was visualized as a heatmap ordered by average-linkage hierarchical clustering; stop codons were labelled with '*', and dendrograms were not displayed to improve readability. Pearson correlations among ENC, GC1, GC2, and GC3 were summarized using a correlation heatmap.
"""



# ----------------------------------------------------------------------
#  Batch-processing and reusable export helpers
# ----------------------------------------------------------------------

BATCH_PLOT_SCHEMES = [
    "default",
    "seaborn-v0_8",
    "seaborn-v0_8-paper",
    "seaborn-v0_8-notebook",
    "tableau-colorblind10",
    "bmh",
    "classic",
    "ggplot",
    "fivethirtyeight",
]

BATCH_HEATMAP_CMAPS = [
    "viridis",
    "cividis",
    "plasma",
    "magma",
    "inferno",
    "turbo",
    "YlGnBu",
    "BuGn",
    "GnBu",
    "PuBuGn",
    "YlOrRd",
    "OrRd",
    "PuRd",
    "BuPu",
    "Greens",
    "Blues",
    "Purples",
    "Oranges",
    "Reds",
    "cubehelix",
]

CORRELATION_CMAPS = [
    "coolwarm",
    "RdBu_r",
    "bwr",
    "seismic",
    "Spectral",
    "Spectral_r",
    "PiYG",
    "PRGn",
    "BrBG",
    "PuOr",
    "RdYlBu",
    "RdYlGn",
    "twilight",
    "twilight_shifted",
]

GENBANK_EXTENSIONS = {".gb", ".gbk", ".genbank"}
BATCH_FILE_LIMIT = 100

# Scatter-plot visibility constants.
# Strong marker outlines and high opacity prevent pale academic colormaps from
# making gene/CDS points invisible on light or grey plot backgrounds.
VISIBLE_SCATTER_EDGE = "#1f2937"
VISIBLE_SCATTER_SIZE = 74
VISIBLE_SCATTER_ALPHA = 0.96
VISIBLE_SCATTER_LINEWIDTH = 0.55
VISIBLE_SCATTER_ZORDER = 3

BATCH_SCATTER_CMAPS = [
    "viridis",
    "cividis",
    "plasma",
    "magma",
    "inferno",
    "turbo",
    "YlGnBu",
    "YlOrRd",
    "BuGn",
    "GnBu",
    "PuBuGn",
    "BuPu",
    "PuRd",
    "Greens",
    "Blues",
    "Purples",
    "Oranges",
    "Reds",
    "cubehelix",
    "rainbow",
    "nipy_spectral",
    "gist_earth",
    "terrain",
    "tab10",
    "tab20",
]


def sanitize_filename(name, fallback="ChloroCodon_Result"):
    """Return a filesystem-safe file/folder name while preserving readability."""
    clean = re.sub(r"[^A-Za-z0-9._ -]+", "_", str(name or "")).strip()
    clean = re.sub(r"\s+", "_", clean)
    clean = clean.strip("._-")
    return clean or fallback


def unique_folder(parent_dir, base_name):
    """Create and return a unique folder path under parent_dir."""
    parent = Path(parent_dir)
    base = sanitize_filename(base_name)
    candidate = parent / base
    counter = 2
    while candidate.exists():
        candidate = parent / f"{base}_{counter}"
        counter += 1
    candidate.mkdir(parents=True, exist_ok=False)
    return candidate


def analyze_genbank_file(file_path, min_codons=30, table_id=11):
    """Analyze one GenBank file and return the complete ChloroCodon data dictionary."""
    # Reinitialize codon maps for the selected NCBI table on every analysis run.
    # Without this, changing the genetic-code table in the sidebar could leave
    # the global synonymous-codon families from a previous run.
    table_id = int(table_id)
    init_codon_maps(table_id)
    records, qc_df, duplicate_df = parse_genbank(file_path, min_codons=int(min_codons), table_id=table_id)
    if not records:
        raise ValueError("No valid CDS found after filtering.")

    metrics_list = [compute_per_gene_metrics(rec, table_id=table_id) for rec in records]
    genome_stats = aggregate_genome_stats(metrics_list)
    overall_rscu = aggregate_rscu(records, table_id=table_id)
    rscu_publication_df = build_rscu_publication_table(records, table_id=table_id)
    amino_acid_usage_df = build_amino_acid_usage_table(records, table_id=table_id)
    stop_usage_df = build_stop_codon_usage_table(records, table_id=table_id)
    per_gene_rscu_df = build_per_gene_rscu_matrix(metrics_list)
    rscu_heatmap_df = build_rscu_heatmap_matrix(metrics_list, table_id=table_id)
    correlation_metrics_df = build_correlation_metrics_table(metrics_list)
    correlation_matrix_df = build_bias_correlation_matrix(correlation_metrics_df, method="pearson")
    opt_codons, hf_codons, delta_rscu = optimal_codon_analysis(metrics_list)

    cai_weights = build_cai_weights_sharp_li(metrics_list, top_frac=0.10, table_id=table_id)
    genome_cai = calculate_cai_sharp_li(metrics_list, cai_weights)
    genome_cbi = calculate_cbi_wright(metrics_list, opt_codons, table_id=table_id)
    genome_fop = calculate_fop(metrics_list, opt_codons, table_id=table_id)

    genome_stats["CAI"] = genome_cai
    genome_stats["CBI"] = genome_cbi
    genome_stats["FOP"] = genome_fop

    for m in metrics_list:
        m["CAI"] = calculate_cai_sharp_li_from_counts(m["codon_counts"], cai_weights)
        m["CBI"] = calculate_cbi_wright_from_counts(m["codon_counts"], opt_codons, table_id=table_id)
        m["FOP"] = calculate_fop_from_counts(m["codon_counts"], opt_codons, table_id=table_id)

    df_enc = pd.DataFrame([{
        "GC3": m["GC3"] * 100 if not np.isnan(m["GC3"]) else np.nan,
        "ENC": m["ENC"] if not np.isnan(m["ENC"]) else np.nan,
        "total_codons": m["total_codons"],
        "label": m["label"],
    } for m in metrics_list]).dropna()

    df_pr2 = pd.DataFrame([{
        "AT_bias": m["AT_bias"],
        "GC_bias": m["GC_bias"],
        "label": m["label"],
    } for m in metrics_list]).dropna()

    df_neutral = pd.DataFrame([{
        "GC3": m["GC3"] if not np.isnan(m["GC3"]) else np.nan,
        "GC12": m["GC12"] if not np.isnan(m["GC12"]) else np.nan,
        "label": m["label"],
    } for m in metrics_list]).dropna()

    codon_list = CODON_ORDER
    X = []
    labels = []
    for m in metrics_list:
        total = m["total_codons"]
        if total == 0:
            continue
        counts = [m["codon_counts"].get(c, 0) / total for c in codon_list]
        X.append(counts)
        labels.append(m["label"])

    if len(X) < 2:
        coa_df = pd.DataFrame({"Axis1": [], "Axis2": [], "label": []})
        coa_explained = [0, 0]
    else:
        X_arr = np.array(X)
        X_arr = np.nan_to_num(X_arr)
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_arr)
        pca = PCA(n_components=2)
        coords = pca.fit_transform(X_scaled)
        coa_df = pd.DataFrame({
            "Axis1": coords[:, 0],
            "Axis2": coords[:, 1],
            "label": labels,
        })
        coa_explained = pca.explained_variance_ratio_ * 100

    data = {
        "records": records,
        "metrics_list": metrics_list,
        "genome_stats": genome_stats,
        "overall_rscu": overall_rscu,
        "rscu_publication_df": rscu_publication_df,
        "amino_acid_usage_df": amino_acid_usage_df,
        "stop_usage_df": stop_usage_df,
        "per_gene_rscu_df": per_gene_rscu_df,
        "rscu_heatmap_df": rscu_heatmap_df,
        "correlation_metrics_df": correlation_metrics_df,
        "correlation_matrix_df": correlation_matrix_df,
        "qc_df": qc_df,
        "duplicate_df": duplicate_df,
        "methods_text": generate_methods_text({"genome_stats": genome_stats, "qc_df": qc_df}, min_codons=int(min_codons), file_path=file_path),
        "opt_codons": opt_codons,
        "hf_codons": hf_codons,
        "delta_rscu": delta_rscu,
        "cai_weights": cai_weights,
        "df_enc": df_enc,
        "df_pr2": df_pr2,
        "df_neutral": df_neutral,
        "coa_df": coa_df,
        "coa_explained": coa_explained,
        "input_file": str(file_path),
    }
    return data


def build_summary_text_from_data(data):
    """Build the same summary text used in the single-file GUI summary tab."""
    stats = data.get("genome_stats", {}) if data else {}
    if not stats:
        return "No data available."

    preferred_order = [
        "Total_bp", "Total_codons", "Num_genes",
        "A%", "T%", "C%", "G%", "GC%",
        "GC1%", "GC2%", "GC3%",
        "A3s", "T3s", "C3s", "G3s", "GC3s",
        "Avg_ENC", "CAI", "CBI", "FOP",
    ]

    lines = ["GENOME AND CODON-BIAS SUMMARY", "=" * 72]
    for key in preferred_order:
        if key not in stats:
            continue
        val = stats[key]
        if isinstance(val, float) or isinstance(val, np.floating):
            lines.append(f"{key:20s}: {val:.6f}")
        else:
            lines.append(f"{key:20s}: {val}")

    for key, val in stats.items():
        if key in preferred_order:
            continue
        if isinstance(val, float) or isinstance(val, np.floating):
            lines.append(f"{key:20s}: {val:.6f}")
        else:
            lines.append(f"{key:20s}: {val}")

    lines.append("")
    lines.append("INDEX DEFINITIONS")
    lines.append("=" * 72)
    lines.append("CAI : Sharp & Li codon adaptation index using low-ENC CDS as reference.")
    lines.append("CBI : Wright 1990 codon bias index using detected optimal codons.")
    lines.append("FOP : Frequency of optimal codons.")
    lines.append("A3s/T3s/C3s/G3s/GC3s: synonymous third-position composition.")

    qc = data.get("qc_df")
    if qc is not None:
        accepted = int((qc["status"] == "accepted").sum()) if not qc.empty and "status" in qc.columns else 0
        skipped = int(len(qc) - accepted)
        lines.append("")
        lines.append("QC SUMMARY")
        lines.append("=" * 72)
        lines.append(f"Accepted CDS before duplicate removal : {accepted}")
        lines.append(f"Skipped CDS                         : {skipped}")
        lines.append(f"CDS used after duplicate removal    : {stats.get('Num_genes', '')}")

    return "\n".join(lines)


def build_output_tables_from_data(data):
    """Create all tabular outputs in one place for CSV and Excel export."""
    tables = {}
    stats = data.get("genome_stats", {}) if data else {}
    tables["summary"] = pd.DataFrame([stats])

    df_rscu = data.get("rscu_publication_df")
    if df_rscu is not None:
        tables["rscu"] = df_rscu

    aa_df = data.get("amino_acid_usage_df")
    if aa_df is not None:
        tables["amino_acid_usage"] = aa_df

    stop_df = data.get("stop_usage_df")
    if stop_df is not None:
        tables["stop_codon_usage"] = stop_df

    duplicate_df = data.get("duplicate_df")
    if duplicate_df is not None:
        tables["duplicate_removed"] = duplicate_df

    metrics_list = data.get("metrics_list", [])
    if metrics_list:
        tables["per_gene_metrics"] = pd.DataFrame([{
            "label": m["label"], "gene": m["gene"],
            "terminal_stop_codon": m.get("terminal_stop_codon", ""),
            "GC": m["GC"], "GC1": m["GC1"], "GC2": m["GC2"], "GC3": m["GC3"], "GC12": m["GC12"],
            "ENC": m["ENC"], "ENC_expected": m.get("ENC_exp", np.nan),
            "total_codons": m["total_codons"],
            "A3": m["A3"], "T3": m["T3"], "G3": m["G3"], "C3": m["C3"],
            "AT_bias": m["AT_bias"], "GC_bias": m["GC_bias"],
            "CAI": m.get("CAI", np.nan), "CBI": m.get("CBI", np.nan), "FOP": m.get("FOP", np.nan),
        } for m in metrics_list])

    qc = data.get("qc_df")
    if qc is not None:
        tables["qc"] = qc

    opt = data.get("opt_codons", [])
    tables["optimal_codons"] = pd.DataFrame({"Optimal": opt})

    hf = data.get("hf_codons", [])
    tables["high_freq_codons"] = pd.DataFrame({"High_Frequency": hf})

    delta = data.get("delta_rscu", pd.Series(dtype=float))
    if delta is not None and not delta.empty:
        tables["delta_rscu"] = delta.rename("Delta_RSCU").reset_index().rename(columns={"index": "Codon"})

    cai_weights = data.get("cai_weights", {})
    if cai_weights:
        tables["cai_weights"] = pd.DataFrame({"Codon": list(cai_weights.keys()), "CAI_weight": list(cai_weights.values())})

    per_gene_rscu = data.get("per_gene_rscu_df")
    if per_gene_rscu is not None:
        tables["per_gene_rscu"] = per_gene_rscu

    rscu_heatmap = data.get("rscu_heatmap_df")
    if rscu_heatmap is not None:
        tables["rscu_heatmap_matrix"] = rscu_heatmap

    corr_metrics = data.get("correlation_metrics_df")
    if corr_metrics is not None:
        tables["correlation_metrics"] = corr_metrics

    corr_matrix = data.get("correlation_matrix_df")
    if corr_matrix is not None:
        tables["correlation_matrix"] = corr_matrix.reset_index().rename(columns={"index": "Metric"})

    df_enc = data.get("df_enc")
    if df_enc is not None:
        tables["enc_plot_data"] = df_enc
    df_pr2 = data.get("df_pr2")
    if df_pr2 is not None:
        tables["pr2_plot_data"] = df_pr2
    df_neutral = data.get("df_neutral")
    if df_neutral is not None:
        tables["neutrality_plot_data"] = df_neutral
    coa_df = data.get("coa_df")
    if coa_df is not None:
        tables["coa_scores"] = coa_df

    return tables


def save_excel_workbook_generic(tables, out_path):
    """Save all tables to one Excel workbook."""
    try:
        with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
            for sheet_name, df in tables.items():
                safe_name = re.sub(r"[^A-Za-z0-9_]+", "_", str(sheet_name))[:31]
                df.to_excel(writer, sheet_name=safe_name, index=False)
    except Exception:
        with pd.ExcelWriter(out_path) as writer:
            for sheet_name, df in tables.items():
                safe_name = re.sub(r"[^A-Za-z0-9_]+", "_", str(sheet_name))[:31]
                df.to_excel(writer, sheet_name=safe_name, index=False)


def get_plot_palette_by_name(name="default", cmap_override=None):
    palette = dict(PLOT_PALETTES.get(name, PLOT_PALETTES["default"]))
    if palette.get("figure_bg") == "#151515":
        palette = dict(PLOT_PALETTES["default"])
    if cmap_override:
        palette["cmap"] = cmap_override
    return palette


def style_axes_static(fig, ax, palette):
    fig.patch.set_facecolor(palette["figure_bg"])
    ax.set_facecolor(palette["axes_bg"])
    ax.title.set_color(palette["text"])
    ax.xaxis.label.set_color(palette["text"])
    ax.yaxis.label.set_color(palette["text"])
    ax.tick_params(axis="both", colors=palette["text"], labelsize=9)
    for spine in ax.spines.values():
        spine.set_color(palette["text"])
        spine.set_linewidth(0.8)
    ax.grid(True, linestyle="--", alpha=0.35, color=palette["grid"])


def style_colorbar_static(cbar, palette):
    cbar.ax.yaxis.label.set_color(palette["text"])
    cbar.ax.tick_params(colors=palette["text"], labelsize=8)
    cbar.outline.set_edgecolor(palette["text"])


def label_top_points_static(ax, df, x_col, y_col, label_top=0, palette=None, score_col=None):
    try:
        n = int(label_top)
    except Exception:
        n = 0
    if n <= 0 or df is None or df.empty or "label" not in df.columns:
        return
    palette = palette or get_plot_palette_by_name("default")
    plot_df = df.copy()
    if score_col is None or score_col not in plot_df.columns:
        plot_df["_score"] = np.arange(len(plot_df), 0, -1)
        score_col = "_score"
    plot_df = plot_df.dropna(subset=[x_col, y_col, score_col])
    if plot_df.empty:
        return
    plot_df = plot_df.sort_values(score_col, ascending=False).head(n)
    for _, row in plot_df.iterrows():
        ax.annotate(str(row.get("label", "")), (row[x_col], row[y_col]),
                    fontsize=8, color=palette["text"], alpha=0.88,
                    xytext=(4, 4), textcoords="offset points",
                    bbox=dict(boxstyle="round,pad=0.15", facecolor=palette["axes_bg"],
                              edgecolor="none", alpha=0.65))


def make_batch_figure(width, height, palette):
    fig = Figure(figsize=(max(float(width), 3.0), max(float(height), 2.5)), dpi=100)
    ax = fig.add_subplot(111)
    style_axes_static(fig, ax, palette)
    return fig, ax


def create_enc_figure(data, width, height, label_top, palette):
    df = data.get("df_enc")
    fig, ax = make_batch_figure(width, height, palette)
    if df is None or df.empty:
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", transform=ax.transAxes, color=palette["text"])
        fig.tight_layout()
        return fig
    sc = ax.scatter(df["GC3"], df["ENC"], c=df["total_codons"], cmap=palette["cmap"],
                    s=VISIBLE_SCATTER_SIZE, edgecolors=VISIBLE_SCATTER_EDGE,
                    linewidth=VISIBLE_SCATTER_LINEWIDTH, alpha=VISIBLE_SCATTER_ALPHA,
                    zorder=VISIBLE_SCATTER_ZORDER)
    x_curve = np.linspace(0.001, 0.999, 500) * 100
    s = x_curve / 100
    y_curve = 2.0 + s + 29.0 / (s**2 + (1 - s)**2)
    ax.plot(x_curve, y_curve, linestyle="--", color=palette["line"], linewidth=2.0, label="Expected ENC")
    ax.set_xlabel("GC3 (%)", fontsize=11)
    ax.set_ylabel("Effective Number of Codons (ENC)", fontsize=11)
    ax.set_title("ENC-GC3 Plot", fontsize=15, fontweight="bold", pad=12)
    ax.set_xlim(0, 100)
    ax.set_ylim(18, 63)
    enc_label_df = df.copy()
    s_obs = enc_label_df["GC3"] / 100.0
    enc_label_df["deviation"] = (enc_label_df["ENC"] - (2.0 + s_obs + 29.0 / (s_obs**2 + (1 - s_obs)**2))).abs()
    label_top_points_static(ax, enc_label_df, "GC3", "ENC", label_top=label_top, palette=palette, score_col="deviation")
    style_axes_static(fig, ax, palette)
    ax.legend(loc="upper right", frameon=True)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("CDS length (codons)")
    style_colorbar_static(cbar, palette)
    fig.tight_layout()
    return fig


def create_pr2_figure(data, width, height, label_top, palette):
    df = data.get("df_pr2")
    fig, ax = make_batch_figure(width, height, palette)
    if df is None or df.empty:
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", transform=ax.transAxes, color=palette["text"])
        fig.tight_layout()
        return fig
    lengths = np.ones(len(df))
    if data.get("df_enc") is not None:
        try:
            lengths = data["df_enc"]["total_codons"].iloc[:len(df)].to_numpy()
        except Exception:
            lengths = np.ones(len(df))
    if len(lengths) != len(df):
        lengths = np.ones(len(df))
    sc = ax.scatter(df["AT_bias"], df["GC_bias"], c=lengths, cmap=palette["cmap"],
                    s=VISIBLE_SCATTER_SIZE, alpha=VISIBLE_SCATTER_ALPHA,
                    edgecolors=VISIBLE_SCATTER_EDGE, linewidth=VISIBLE_SCATTER_LINEWIDTH,
                    zorder=VISIBLE_SCATTER_ZORDER)
    ax.axhline(0.5, linestyle="--", color=palette["accent"], linewidth=1.25)
    ax.axvline(0.5, linestyle="--", color=palette["accent"], linewidth=1.25)
    ax.set_xlabel("A3 / (A3 + T3)", fontsize=11)
    ax.set_ylabel("G3 / (G3 + C3)", fontsize=11)
    ax.set_title("PR2 Bias Plot", fontsize=15, fontweight="bold", pad=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    pr2_label_df = df.copy()
    pr2_label_df["distance"] = np.sqrt((pr2_label_df["AT_bias"] - 0.5)**2 + (pr2_label_df["GC_bias"] - 0.5)**2)
    label_top_points_static(ax, pr2_label_df, "AT_bias", "GC_bias", label_top=label_top, palette=palette, score_col="distance")
    style_axes_static(fig, ax, palette)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("CDS length (codons)")
    style_colorbar_static(cbar, palette)
    fig.tight_layout()
    return fig


def create_neutrality_figure(data, width, height, label_top, palette):
    """Create Neutrality plot with user-selectable theme and point colormap."""
    df = data.get("df_neutral")
    fig, ax = make_batch_figure(width, height, palette)
    if df is None or df.empty:
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", transform=ax.transAxes, color=palette["text"])
        fig.tight_layout()
        return fig

    clean = df.dropna()
    res = None
    color_values = pd.to_numeric(df.get("GC3", pd.Series(np.zeros(len(df)))), errors="coerce").fillna(0.0)
    color_label = "GC3"

    if len(clean) > 1 and clean["GC3"].nunique() > 1:
        res = linregress(clean["GC3"], clean["GC12"])
        fitted = res.intercept + res.slope * df["GC3"]
        color_values = (df["GC12"] - fitted).abs().fillna(0.0)
        color_label = "Absolute residual"

    sc = ax.scatter(
        df["GC3"], df["GC12"],
        c=color_values, cmap=palette.get("cmap", "viridis"),
        s=VISIBLE_SCATTER_SIZE, alpha=VISIBLE_SCATTER_ALPHA,
        edgecolors=VISIBLE_SCATTER_EDGE, linewidth=VISIBLE_SCATTER_LINEWIDTH,
        zorder=VISIBLE_SCATTER_ZORDER,
    )

    if res is not None:
        x_line = np.linspace(0, 1, 200)
        y_line = res.intercept + res.slope * x_line
        ax.plot(x_line, y_line, color=palette["line"], linewidth=2.0, label=f"R² = {res.rvalue**2:.3f}")
        ax.text(0.05, 0.95,
                f"Slope = {res.slope:.3f}\nIntercept = {res.intercept:.3f}\nR² = {res.rvalue**2:.3f}\np = {res.pvalue:.2e}",
                transform=ax.transAxes, fontsize=9.5, verticalalignment="top",
                color=palette["text"],
                bbox=dict(boxstyle="round,pad=0.35", facecolor=palette["axes_bg"],
                          edgecolor=palette["text"], alpha=0.78))
        ax.legend(loc="lower right", frameon=True)

    ax.axhline(0.5, linestyle="--", color=palette["accent"], linewidth=1.0, alpha=0.75)
    ax.axvline(0.5, linestyle="--", color=palette["accent"], linewidth=1.0, alpha=0.75)
    ax.set_xlabel("GC3", fontsize=11)
    ax.set_ylabel("GC12 (mean GC1 & GC2)", fontsize=11)
    ax.set_title("Neutrality Plot", fontsize=15, fontweight="bold", pad=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    if res is not None:
        neutral_label_df = df.copy()
        neutral_label_df["residual"] = (neutral_label_df["GC12"] - (res.intercept + res.slope * neutral_label_df["GC3"])).abs()
        label_top_points_static(ax, neutral_label_df, "GC3", "GC12", label_top=label_top, palette=palette, score_col="residual")

    style_axes_static(fig, ax, palette)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(color_label)
    style_colorbar_static(cbar, palette)
    fig.tight_layout()
    return fig


def create_coa_figure(data, width, height, label_top, palette):
    """Create COA/PCA figure with user-selectable theme and point colormap."""
    df = data.get("coa_df")
    explained = data.get("coa_explained", [0, 0])
    fig, ax = make_batch_figure(width, height, palette)
    if df is None or df.empty:
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", transform=ax.transAxes, color=palette["text"])
        fig.tight_layout()
        return fig

    color_values = np.sqrt(pd.to_numeric(df["Axis1"], errors="coerce")**2 + pd.to_numeric(df["Axis2"], errors="coerce")**2).fillna(0.0)
    sc = ax.scatter(
        df["Axis1"], df["Axis2"],
        c=color_values, cmap=palette.get("cmap", "viridis"),
        s=VISIBLE_SCATTER_SIZE, alpha=VISIBLE_SCATTER_ALPHA,
        edgecolors=VISIBLE_SCATTER_EDGE, linewidth=VISIBLE_SCATTER_LINEWIDTH,
        zorder=VISIBLE_SCATTER_ZORDER,
    )
    if len(df) < 50:
        for _, row in df.iterrows():
            ax.annotate(row["label"], (row["Axis1"], row["Axis2"]),
                        fontsize=8, alpha=0.75, color=palette["text"],
                        xytext=(3, 3), textcoords="offset points")
    ax.axhline(0, linestyle="--", color=palette["accent"], linewidth=1.0, alpha=0.75)
    ax.axvline(0, linestyle="--", color=palette["accent"], linewidth=1.0, alpha=0.75)
    coa_label_df = df.copy()
    coa_label_df["distance"] = color_values
    label_top_points_static(ax, coa_label_df, "Axis1", "Axis2", label_top=label_top, palette=palette, score_col="distance")
    ax.set_xlabel(f"Axis 1 ({explained[0]:.2f}% variance)", fontsize=11)
    ax.set_ylabel(f"Axis 2 ({explained[1]:.2f}% variance)", fontsize=11)
    ax.set_title("Correspondence Analysis (PCA approximation)", fontsize=15, fontweight="bold", pad=12)
    style_axes_static(fig, ax, palette)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Distance from origin")
    style_colorbar_static(cbar, palette)
    fig.tight_layout()
    return fig


def create_rscu_heatmap_figure(data, width, height, palette):
    df = data.get("rscu_heatmap_df")
    fig, ax = make_batch_figure(width, height, palette)
    if df is None or df.empty:
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", transform=ax.transAxes, color=palette["text"])
        fig.tight_layout()
        return fig

    heatmap_codons = [c for c in CODON_ORDER if c in df.columns]
    heatmap_data = df[heatmap_codons].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    labels = df["label"].astype(str).tolist() if "label" in df.columns else [f"Gene_{i+1}" for i in range(len(heatmap_data))]
    if heatmap_data.shape[0] < 1 or heatmap_data.shape[1] < 1:
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", transform=ax.transAxes, color=palette["text"])
        fig.tight_layout()
        return fig

    values = heatmap_data.to_numpy(dtype=float)
    if values.shape[0] >= 2:
        try:
            row_link = linkage(pdist(values, metric="euclidean"), method="average")
            row_order = dendrogram(row_link, no_plot=True)["leaves"]
        except Exception:
            row_order = list(range(values.shape[0]))
    else:
        row_order = list(range(values.shape[0]))

    if values.shape[1] >= 2:
        try:
            col_link = linkage(pdist(values.T, metric="euclidean"), method="average")
            col_order = dendrogram(col_link, no_plot=True)["leaves"]
        except Exception:
            col_order = list(range(values.shape[1]))
    else:
        col_order = list(range(values.shape[1]))

    ordered_values = values[np.ix_(row_order, col_order)]
    ordered_rows = [labels[i] for i in row_order]
    ordered_cols_raw = [heatmap_data.columns[i] for i in col_order]
    ordered_cols = [f"{c}*" if c in STOP_CODONS else c for c in ordered_cols_raw]

    fig.clf()
    ax = fig.add_subplot(111)
    fig.patch.set_facecolor(palette["figure_bg"])
    ax.set_facecolor(palette["axes_bg"])
    ax.grid(False)
    im = ax.imshow(ordered_values, aspect="auto", interpolation="nearest", cmap=palette["cmap"])
    ax.set_title("RSCU Heatmap Across Genes", fontsize=12.5, fontweight="bold", color=palette["text"], pad=7)

    n_rows = len(ordered_rows)
    n_cols = len(ordered_cols)
    x_font = 4.4 if n_cols >= 64 else 5.0
    if n_rows > 100:
        y_font = 2.6
    elif n_rows > 80:
        y_font = 3.0
    elif n_rows > 60:
        y_font = 3.4
    elif n_rows > 40:
        y_font = 4.2
    else:
        y_font = 5.2

    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(ordered_cols, rotation=90, fontsize=x_font, color=palette["text"])
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(ordered_rows, fontsize=y_font, color=palette["text"])
    ax.tick_params(axis="both", colors=palette["text"], length=1.0, pad=1.0)
    ax.set_ylabel("Genes/CDS (all labels shown)", fontsize=8.8, color=palette["text"])
    ax.set_xlabel("Codons (* = stop codon)", fontsize=8.8, color=palette["text"])

    for spine in ax.spines.values():
        spine.set_color(palette["text"])
        spine.set_linewidth(0.45)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.012)
    cbar.set_label("RSCU", color=palette["text"])
    style_colorbar_static(cbar, palette)

    max_label_len = max((len(x) for x in ordered_rows), default=4)
    left_margin = min(0.18, max(0.075, 0.060 + max_label_len * 0.0032))
    fig.subplots_adjust(left=left_margin, right=0.965, bottom=0.165, top=0.905)
    return fig


def create_correlation_figure(data, width, height, palette):
    """Create a compact correlation heatmap with the colour bar close to the matrix.

    The colour bar is appended with Matplotlib's axes divider instead of a far
    fixed-position axis. This keeps the legend visually attached to the matrix
    in both batch exports and single-mode preview figures.
    """
    corr = data.get("correlation_matrix_df")
    fig, ax = make_batch_figure(width, height, palette)
    if corr is None or corr.empty:
        ax.text(0.5, 0.5, "No data available", ha="center", va="center",
                transform=ax.transAxes, color=palette["text"])
        fig.tight_layout()
        return fig

    corr = corr.copy().astype(float)
    labels = list(corr.columns)
    values = corr.to_numpy(dtype=float)
    values_masked = np.ma.masked_invalid(values)

    fig.clf()
    fig.patch.set_facecolor(palette["figure_bg"])
    ax = fig.add_subplot(111)
    ax.set_facecolor(palette["axes_bg"])

    # Give the square matrix enough room for labels and title, while keeping
    # the colorbar close to the actual heatmap body.
    fig.subplots_adjust(left=0.135, right=0.835, bottom=0.175, top=0.850)

    im = ax.imshow(values_masked, vmin=-1, vmax=1, cmap=palette.get("cmap", "coolwarm"),
                   interpolation="nearest", aspect="equal")
    ax.set_title("Correlation Matrix of Codon-Bias Indices",
                 fontsize=13, fontweight="bold", color=palette["text"], pad=10)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", color=palette["text"], fontsize=10)
    ax.set_yticklabels(labels, color=palette["text"], fontsize=10)
    ax.tick_params(axis="both", colors=palette["text"])
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_color(palette["text"])
        spine.set_linewidth(0.8)

    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            val = values[i, j]
            txt = "" if np.isnan(val) else f"{val:.2f}"
            ax.text(j, i, txt, ha="center", va="center", color="#111111", fontsize=10,
                    bbox=dict(boxstyle="round,pad=0.12", facecolor="#ffffff",
                              edgecolor="none", alpha=0.70))

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4.5%", pad=0.08)
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Pearson r", color=palette["text"], labelpad=7)
    style_colorbar_static(cbar, palette)
    return fig

def create_batch_figures(data, width, height, label_top=0, palette_name="default", cmap_name=None,
                         scatter_cmap_name=None, correlation_palette_name=None, correlation_cmap_name=None,
                         neutrality_palette_name=None, neutrality_cmap_name=None,
                         coa_palette_name=None, coa_cmap_name=None):
    """Create all figures with separate user-selectable themes for key plots."""
    heatmap_palette = get_plot_palette_by_name(palette_name, cmap_override=cmap_name)
    scatter_palette = get_plot_palette_by_name(
        palette_name,
        cmap_override=scatter_cmap_name or cmap_name or "viridis"
    )
    correlation_palette = get_plot_palette_by_name(
        correlation_palette_name or palette_name,
        cmap_override=correlation_cmap_name or "coolwarm"
    )
    neutrality_palette = get_plot_palette_by_name(
        neutrality_palette_name or palette_name,
        cmap_override=neutrality_cmap_name or scatter_cmap_name or "viridis"
    )
    coa_palette = get_plot_palette_by_name(
        coa_palette_name or palette_name,
        cmap_override=coa_cmap_name or scatter_cmap_name or "plasma"
    )
    return {
        "RSCU_Heatmap": create_rscu_heatmap_figure(data, width, height, heatmap_palette),
        "Correlations": create_correlation_figure(data, width, height, correlation_palette),
        "ENC": create_enc_figure(data, width, height, label_top, scatter_palette),
        "PR2": create_pr2_figure(data, width, height, label_top, scatter_palette),
        "Neutrality": create_neutrality_figure(data, width, height, label_top, neutrality_palette),
        "COA": create_coa_figure(data, width, height, label_top, coa_palette),
    }


def save_analysis_package(data, output_dir, prefix_name, formats, fig_width=8.0, fig_height=6.0,
                          label_top=0, palette_name="default", cmap_name=None,
                          scatter_cmap_name=None, correlation_palette_name=None, correlation_cmap_name=None,
                          neutrality_palette_name=None, neutrality_cmap_name=None,
                          coa_palette_name=None, coa_cmap_name=None, save_figures=True):
    """Save the complete ChloroCodon output package for one analyzed file."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    safe_prefix = sanitize_filename(prefix_name, fallback="ChloroCodon_results")
    prefix = out_dir / safe_prefix

    tables = build_output_tables_from_data(data)
    for name, df in tables.items():
        csv_path = prefix.parent / f"{prefix.name}.{name}.csv"
        df.to_csv(csv_path, index=False)

    save_excel_workbook_generic(tables, prefix.parent / f"{prefix.name}.all_results.xlsx")

    summary_text = build_summary_text_from_data(data)
    with open(prefix.with_suffix(".summary.txt"), "w", encoding="utf-8") as f:
        f.write(summary_text)

    methods_text = data.get("methods_text", "")
    if methods_text:
        with open(prefix.with_suffix(".methods.txt"), "w", encoding="utf-8") as f:
            f.write(str(methods_text))

    formats = list(formats or [])
    if save_figures and formats:
        figures = create_batch_figures(data, fig_width, fig_height, label_top=label_top,
                                       palette_name=palette_name, cmap_name=cmap_name,
                                       scatter_cmap_name=scatter_cmap_name,
                                       correlation_palette_name=correlation_palette_name,
                                       correlation_cmap_name=correlation_cmap_name,
                                       neutrality_palette_name=neutrality_palette_name,
                                       neutrality_cmap_name=neutrality_cmap_name,
                                       coa_palette_name=coa_palette_name,
                                       coa_cmap_name=coa_cmap_name)
        for name, fig in figures.items():
            if fig is None:
                continue
            try:
                for fmt in formats:
                    ext = "tiff" if fmt == "tiff" else fmt
                    fname = prefix.parent / f"{prefix.name}_{name}.{ext}"
                    fig.savefig(fname, dpi=EXPORT_DPI, bbox_inches="tight")
            finally:
                try:
                    fig.clf()
                    plt.close(fig)
                except Exception:
                    pass

    return tables


# ----------------------------------------------------------------------
#  Streamlit Application
# ----------------------------------------------------------------------

st.set_page_config(
    page_title="ChloroCodon – Codon Usage Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


def _local_image_data_uri(image_name):
    """Return a base64 data URI for an image placed beside this script or in the run folder."""
    try:
        script_dir = Path(__file__).resolve().parent
    except Exception:
        script_dir = Path.cwd()

    search_paths = [script_dir / image_name, Path.cwd() / image_name]
    image_path = next((path for path in search_paths if path.exists() and path.is_file()), None)
    if image_path is None:
        return None

    suffix = image_path.suffix.lower().replace(".", "") or "png"
    mime = "jpeg" if suffix in {"jpg", "jpeg"} else suffix
    try:
        encoded = base64.b64encode(image_path.read_bytes()).decode("utf-8")
        return f"data:image/{mime};base64,{encoded}"
    except Exception:
        return None


def _render_global_css():
    """Inject the premium Streamlit theme without changing app logic."""
    # Page texture can still use background.png, but the main hero/banner uses
    # banner_image.png from the same working folder as this script.
    bg_uri = _local_image_data_uri(SUMMARY_BACKGROUND_IMAGE_NAME)
    banner_uri = _local_image_data_uri(BANNER_IMAGE_NAME)

    if bg_uri:
        app_background = f"""
            background-image:
                radial-gradient(circle at top left, rgba(223, 243, 232, 0.88), rgba(243, 248, 245, 0.92) 34%, rgba(247, 251, 248, 0.97) 72%),
                url("{bg_uri}");
            background-size: cover;
            background-position: center;
            background-attachment: fixed;
        """
    else:
        app_background = """
            background:
                radial-gradient(circle at top left, rgba(223, 243, 232, 0.90), transparent 30%),
                linear-gradient(135deg, #f3f8f5 0%, #fbfdfb 48%, #eef6f0 100%);
        """

    if banner_uri:
        hero_background = f"""
            background-image:
                linear-gradient(135deg, rgba(7, 26, 18, 0.72), rgba(31, 122, 77, 0.50)),
                url("{banner_uri}");
            background-size: cover;
            background-position: center;
        """
    else:
        hero_background = """
            background:
                radial-gradient(circle at top left, rgba(255,255,255,0.24), transparent 26%),
                linear-gradient(135deg, #071a12 0%, #0e3d28 44%, #1f7a4d 100%);
        """

    st.markdown(f"""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap');

    html, body, [class*="css"] {{
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    }}

    .stApp {{
        {app_background}
        color: #13251c;
    }}

    .main .block-container {{
        padding-top: .38rem;
        padding-bottom: 2.2rem;
        max-width: 1440px;
    }}

    section[data-testid="stSidebar"] {{
        background: linear-gradient(180deg, rgba(255,255,255,0.96), rgba(246,251,248,0.96));
        border-right: 1px solid rgba(216, 232, 222, 0.95);
        box-shadow: 12px 0 35px rgba(14, 61, 40, 0.055);
    }}

    /* Sidebar remains native/collapsible.
       Do not force width/transform here; otherwise Streamlit cannot close it. */
    section[data-testid="stSidebar"] {{
        z-index: 999998 !important;
    }}

    [data-testid="stSidebarContent"],
    section[data-testid="stSidebar"] > div:first-child {{
        overflow-y: auto !important;
    }}

    section[data-testid="stSidebar"] h1,
    section[data-testid="stSidebar"] h2,
    section[data-testid="stSidebar"] h3 {{
        color: #0e3d28;
        font-weight: 800;
        letter-spacing: -0.02em;
    }}

    section[data-testid="stSidebar"] > div:first-child {{
        padding-top: 1.35rem;
    }}

    .cc-sidebar-title-card {{
        position: relative;
        overflow: hidden;
        padding: 1rem 1rem .95rem 1rem;
        margin: .15rem 0 .95rem 0;
        border-radius: 22px;
        border: 1px solid rgba(216,232,222,0.96);
        background:
            radial-gradient(circle at top left, rgba(223,243,232,0.90), transparent 42%),
            linear-gradient(135deg, rgba(255,255,255,0.96), rgba(247,251,248,0.90));
        box-shadow: 0 16px 36px rgba(14,61,40,0.08);
    }}

    .cc-sidebar-title-card:before {{
        content: "";
        position: absolute;
        right: -26px;
        top: -34px;
        width: 92px;
        height: 92px;
        border-radius: 50%;
        background: rgba(31,122,77,0.10);
    }}

    .cc-sidebar-kicker {{
        display: inline-flex;
        align-items: center;
        gap: .35rem;
        padding: .28rem .62rem;
        margin-bottom: .55rem;
        border-radius: 999px;
        background: rgba(31,122,77,0.10);
        border: 1px solid rgba(31,122,77,0.16);
        color: #1f7a4d;
        font-size: .68rem;
        font-weight: 900;
        letter-spacing: .09em;
        text-transform: uppercase;
    }}

    .cc-sidebar-title {{
        position: relative;
        color: #0e3d28;
        font-size: 1.08rem;
        font-weight: 950;
        line-height: 1.18;
        letter-spacing: -0.035em;
    }}

    .cc-sidebar-subtitle {{
        position: relative;
        color: #627267;
        font-size: .78rem;
        font-weight: 650;
        line-height: 1.45;
        margin-top: .35rem;
    }}

    .cc-sidebar-recovery-note {{
        padding: .72rem .86rem;
        margin: -.35rem 0 1rem 0;
        border-radius: 16px;
        background: rgba(223,243,232,0.52);
        border: 1px solid rgba(31,122,77,0.16);
        color: #385244;
        font-size: .78rem;
        font-weight: 650;
        line-height: 1.45;
    }}

    /* Keep Streamlit's sidebar restore button usable.
       The previous ultra-thin header could make the sidebar arrow difficult
       to click after the browser window was minimized or the sidebar was
       manually collapsed. */
    [data-testid="stHeader"] {{
        height: 3.1rem !important;
        min-height: 3.1rem !important;
        background: rgba(255, 255, 255, 0.62) !important;
        backdrop-filter: blur(10px) !important;
        z-index: 999997 !important;
    }}

    /* Keep Streamlit's native top-right controls visible.
       Do not hide stToolbar/stStatusWidget, because that removes the
       Run/Rerun menu, Settings menu, and Deploy button from localhost. */
    [data-testid="stToolbar"],
    [data-testid="stStatusWidget"],
    [data-testid="stMainMenu"],
    #MainMenu {{
        display: flex !important;
        visibility: visible !important;
        opacity: 1 !important;
        pointer-events: auto !important;
        z-index: 999999 !important;
    }}

    [data-testid="stToolbar"] {{
        right: .85rem !important;
    }}

    [data-testid="stDecoration"] {{
        display: none !important;
    }}

    [data-testid="collapsedControl"],
    [data-testid="stSidebarCollapsedControl"],
    button[title="Open sidebar"],
    button[aria-label="Open sidebar"],
    button[title="Close sidebar"],
    button[aria-label="Close sidebar"] {{
        display: flex !important;
        visibility: visible !important;
        opacity: 1 !important;
        pointer-events: auto !important;
        z-index: 999999 !important;
    }}

    [data-testid="collapsedControl"] {{
        position: fixed !important;
        top: .55rem !important;
        left: .55rem !important;
        width: 2.15rem !important;
        height: 2.15rem !important;
        border-radius: 999px !important;
        background: rgba(255,255,255,0.96) !important;
        border: 1px solid rgba(216,232,222,0.95) !important;
        box-shadow: 0 10px 26px rgba(14,61,40,0.14) !important;
    }}

    .cc-hero {{
        {hero_background}
        position: relative;
        overflow: hidden;
        padding: 1.25rem 1.65rem 1.15rem 1.65rem;
        min-height: 158px;
        border-radius: 24px;
        color: #ffffff;
        margin-bottom: .42rem;
        box-shadow: 0 18px 46px rgba(7, 26, 18, 0.16);
        border: 1px solid rgba(255,255,255,0.16);
    }}

    .cc-hero-compact {{
        min-height: 132px;
        padding: .95rem 1.45rem .86rem 1.45rem;
        margin-bottom: .06rem;
    }}

    .cc-hero-compact p {{
        margin-top: .55rem;
        max-width: 760px;
        line-height: 1.48;
    }}

    .cc-hero-compact + div {{
        margin-top: 0 !important;
    }}


    .cc-processing-title {{
        display: inline-flex;
        align-items: center;
        margin-top: .92rem;
        margin-bottom: .42rem;
        padding: .34rem .78rem;
        border-radius: 999px;
        background: rgba(255,255,255,0.76);
        border: 1px solid rgba(216,232,222,0.95);
        color: #0e3d28;
        font-size: 1.04rem;
        font-weight: 850;
        letter-spacing: -0.025em;
        box-shadow: 0 10px 24px rgba(14,61,40,0.055);
    }}

    .cc-processing-title::before {{
        content: "";
        width: .54rem;
        height: .54rem;
        margin-right: .46rem;
        border-radius: 50%;
        background: linear-gradient(135deg, #1f7a4d, #d99a2b);
        box-shadow: 0 0 0 4px rgba(31,122,77,0.11);
    }}

    .cc-hero:after {{
        content: "";
        position: absolute;
        width: 260px;
        height: 260px;
        right: -90px;
        top: -110px;
        background: radial-gradient(circle, rgba(217,154,43,0.24), transparent 68%);
        pointer-events: none;
    }}

    .cc-hero-kicker {{
        display: inline-flex;
        align-items: center;
        gap: .45rem;
        padding: .35rem .72rem;
        border-radius: 999px;
        background: rgba(255,255,255,0.14);
        border: 1px solid rgba(255,255,255,0.22);
        font-size: .78rem;
        font-weight: 700;
        text-transform: uppercase;
        letter-spacing: .08em;
        color: #e9fff1;
        margin-bottom: .7rem;
    }}

    .cc-hero h1 {{
        margin: 0;
        max-width: 920px;
        font-size: clamp(2rem, 3.5vw, 3.45rem);
        line-height: 0.96;
        letter-spacing: -0.065em;
        font-weight: 850;
    }}

    .cc-hero p {{
        max-width: 800px;
        margin: .72rem 0 0 0;
        color: rgba(242,255,248,0.90);
        font-size: .98rem;
        line-height: 1.55;
    }}

    .cc-hero-strip {{
        display: flex;
        flex-wrap: wrap;
        gap: .7rem;
        margin-top: 1.25rem;
    }}

    .cc-pill {{
        padding: .48rem .78rem;
        border-radius: 999px;
        background: rgba(255, 255, 255, 0.14);
        border: 1px solid rgba(255, 255, 255, 0.22);
        color: #ffffff;
        font-size: .86rem;
        font-weight: 650;
        backdrop-filter: blur(8px);
    }}

    .cc-panel {{
        border: 1px solid rgba(216, 232, 222, 0.95);
        background: rgba(255, 255, 255, 0.84);
        border-radius: 24px;
        padding: 1.15rem 1.2rem;
        margin: .48rem 0 .86rem 0;
        box-shadow: 0 18px 48px rgba(14, 61, 40, 0.07);
        backdrop-filter: blur(18px);
    }}

    .cc-panel-tight {{
        padding: .72rem .88rem;
        border-radius: 18px;
    }}

    .cc-section-title {{
        margin: .2rem 0 .75rem 0;
        color: #0e3d28;
        font-weight: 800;
        letter-spacing: -0.035em;
        font-size: 1.28rem;
    }}

    .cc-muted {{
        color: #627267;
        font-size: .93rem;
        line-height: 1.55;
    }}

    .cc-feature-grid {{
        display: grid;
        grid-template-columns: repeat(4, minmax(0, 1fr));
        gap: .85rem;
        margin: 1rem 0 .2rem 0;
    }}

    .cc-feature-card {{
        background: linear-gradient(180deg, rgba(255,255,255,0.94), rgba(247,251,248,0.92));
        border: 1px solid rgba(216,232,222,0.92);
        border-radius: 20px;
        padding: 1rem;
        box-shadow: 0 14px 38px rgba(14,61,40,0.06);
    }}

    .cc-feature-icon {{
        width: 38px;
        height: 38px;
        display: inline-flex;
        align-items: center;
        justify-content: center;
        border-radius: 13px;
        background: #dff3e8;
        color: #0e3d28;
        font-size: 1.08rem;
        margin-bottom: .7rem;
        box-shadow: inset 0 0 0 1px rgba(31,122,77,0.10);
    }}

    .cc-feature-card b {{
        display: block;
        color: #13251c;
        font-weight: 800;
        margin-bottom: .2rem;
    }}

    .cc-feature-card span {{
        color: #627267;
        font-size: .86rem;
        line-height: 1.48;
    }}

    .cc-metric-card {{
        position: relative;
        overflow: hidden;
        min-height: 112px;
        background: linear-gradient(180deg, rgba(255,255,255,0.98), rgba(247,251,248,0.95));
        border: 1px solid rgba(216,232,222,0.98);
        border-radius: 22px;
        padding: .95rem 1rem;
        box-shadow: 0 18px 42px rgba(14, 61, 40, 0.07);
    }}

    .cc-metric-card:after {{
        content: "";
        position: absolute;
        width: 92px;
        height: 92px;
        right: -38px;
        bottom: -44px;
        border-radius: 50%;
        background: rgba(217,154,43,0.13);
    }}

    .cc-metric-label {{
        color: #627267;
        font-size: .76rem;
        font-weight: 800;
        text-transform: uppercase;
        letter-spacing: .075em;
        margin-bottom: .42rem;
    }}

    .cc-metric-value {{
        color: #0e3d28;
        font-size: 1.72rem;
        font-weight: 850;
        letter-spacing: -0.045em;
        line-height: 1.05;
        word-break: break-word;
    }}

    .cc-metric-note {{
        color: #7c8b81;
        font-size: .78rem;
        margin-top: .35rem;
    }}

    .cc-metric-row-gap {{
        height: 1.15rem;
        min-height: 1.15rem;
    }}

    .cc-status-good {{
        border: 1px solid rgba(31,122,77,0.18);
        background: linear-gradient(135deg, rgba(223,243,232,0.92), rgba(255,255,255,0.92));
        color: #0e3d28;
        border-radius: 18px;
        padding: .85rem 1rem;
        font-weight: 650;
    }}

    .cc-summary-box {{
        background: linear-gradient(180deg, #071a12, #0e3d28);
        color: #e9fff1;
        border-radius: 22px;
        padding: 1.15rem;
        border: 1px solid rgba(255,255,255,0.12);
        box-shadow: 0 22px 48px rgba(7, 26, 18, 0.20);
    }}

    .cc-summary-box pre {{
        margin: 0;
        white-space: pre-wrap;
        word-break: break-word;
        color: #e9fff1;
        font-size: .9rem;
        line-height: 1.55;
        font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', monospace;
    }}

    .cc-footer {{
        position: relative;
        overflow: hidden;
        margin-top: 2rem;
        padding: 1.2rem 1.35rem;
        text-align: center;
        color: #627267;
        border: 1px solid rgba(216,232,222,0.98);
        border-radius: 24px;
        background:
            radial-gradient(circle at 12% 10%, rgba(217,154,43,0.16), transparent 28%),
            linear-gradient(135deg, rgba(255,255,255,0.94), rgba(247,251,248,0.90));
        box-shadow: 0 18px 44px rgba(14,61,40,0.075);
        backdrop-filter: blur(16px);
    }}

    .cc-footer:before {{
        content: "";
        position: absolute;
        inset: 0;
        border-radius: 24px;
        pointer-events: none;
        background: linear-gradient(90deg, rgba(31,122,77,0.10), transparent 34%, rgba(217,154,43,0.12));
    }}

    .cc-footer-inner {{
        position: relative;
        z-index: 1;
        max-width: 980px;
        margin: 0 auto;
    }}

    .cc-footer-brand {{
        color: #0e3d28;
        font-size: 1.05rem;
        font-weight: 900;
        letter-spacing: -0.035em;
        margin-bottom: .35rem;
    }}

    .cc-footer-meta {{
        display: inline-block;
        color: #486056;
        font-size: .82rem;
        font-weight: 700;
        line-height: 1.5;
        padding: .5rem .85rem;
        border-radius: 999px;
        background: rgba(255,255,255,0.70);
        border: 1px solid rgba(216,232,222,0.88);
        box-shadow: inset 0 0 0 1px rgba(255,255,255,0.45);
    }}

    .cc-footer-copy {{
        color: #13251c;
        font-size: .88rem;
        font-weight: 800;
        margin-top: .75rem;
    }}

    .cc-footer-lab {{
        color: #6b7c72;
        font-size: .80rem;
        line-height: 1.45;
        margin-top: .2rem;
    }}

    div[data-testid="stFileUploader"] {{
        background: rgba(255,255,255,0.72);
        border: 1px dashed rgba(31,122,77,0.34);
        border-radius: 22px;
        padding: .55rem .75rem .75rem .75rem;
    }}

    div[data-testid="stFileUploader"] section {{
        border: 0 !important;
        background: transparent !important;
    }}

    div[data-testid="stButton"] > button,
    div[data-testid="stDownloadButton"] > button {{
        border-radius: 999px;
        min-height: 2.75rem;
        border: 1px solid rgba(31,122,77,0.18);
        font-weight: 800;
        letter-spacing: -0.01em;
        box-shadow: 0 12px 24px rgba(31, 122, 77, 0.12);
        transition: all .18s ease;
    }}

    div[data-testid="stButton"] > button:hover,
    div[data-testid="stDownloadButton"] > button:hover {{
        transform: translateY(-1px);
        box-shadow: 0 18px 36px rgba(31, 122, 77, 0.18);
    }}

    div[data-baseweb="tab-list"], div[role="radiogroup"] {{
        gap: .55rem;
    }}

    div[role="radiogroup"] label {{
        background: rgba(255,255,255,0.72);
        border: 1px solid rgba(216,232,222,0.82);
        border-radius: 999px;
        padding: .28rem .78rem;
        box-shadow: 0 8px 18px rgba(14,61,40,0.035);
    }}

    div[data-testid="stRadio"] {{
        margin-top: 0 !important;
        margin-bottom: .48rem !important;
    }}

    div[data-testid="stRadio"] > label {{
        padding-bottom: .15rem !important;
        margin-bottom: .1rem !important;
        color: #0e3d28 !important;
        font-weight: 800 !important;
    }}

    [data-testid="stDataFrame"] {{
        border-radius: 18px;
        overflow: hidden;
        border: 1px solid rgba(216,232,222,0.95);
        box-shadow: 0 18px 42px rgba(14,61,40,0.055);
    }}

    .stAlert {{
        border-radius: 18px;
    }}

    h1, h2, h3 {{
        letter-spacing: -0.035em;
    }}


    .cc-workspace-heading {{
        position: sticky;
        top: .75rem;
        z-index: 5;
        padding: 1rem 1rem .95rem 1rem;
        margin-bottom: .85rem;
        border-radius: 22px;
        border: 1px solid rgba(216,232,222,0.96);
        background:
            radial-gradient(circle at top left, rgba(223,243,232,0.82), transparent 38%),
            linear-gradient(135deg, rgba(255,255,255,0.94), rgba(247,251,248,0.88));
        box-shadow: 0 18px 40px rgba(14,61,40,0.075);
        backdrop-filter: blur(16px);
    }}

    .cc-workspace-kicker {{
        display: inline-flex;
        align-items: center;
        gap: .35rem;
        padding: .28rem .62rem;
        margin-bottom: .55rem;
        border-radius: 999px;
        background: rgba(31,122,77,0.10);
        border: 1px solid rgba(31,122,77,0.16);
        color: #1f7a4d;
        font-size: .68rem;
        font-weight: 900;
        letter-spacing: .09em;
        text-transform: uppercase;
    }}

    .cc-workspace-title {{
        color: #0e3d28;
        font-size: 1.08rem;
        font-weight: 950;
        line-height: 1.18;
        letter-spacing: -0.035em;
    }}

    .cc-workspace-subtitle {{
        color: #627267;
        font-size: .78rem;
        font-weight: 650;
        line-height: 1.45;
        margin-top: .35rem;
    }}

    div[data-testid="stExpander"] {{
        border: 1px solid rgba(216,232,222,0.96) !important;
        border-radius: 18px !important;
        background: rgba(255,255,255,0.78) !important;
        box-shadow: 0 12px 28px rgba(14,61,40,0.045);
        overflow: hidden;
    }}

    @media (max-width: 980px) {{
        .cc-feature-grid {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
        .cc-hero {{ padding: 1.15rem 1.25rem; border-radius: 22px; }}
    }}

    @media (max-width: 640px) {{
        .cc-feature-grid {{ grid-template-columns: 1fr; }}
        .cc-hero h1 {{ font-size: 2rem; }}
    }}
</style>
""", unsafe_allow_html=True)


def _render_premium_hero():
    st.markdown(
        """
        <div class="cc-hero cc-hero-compact">
            <h1>ChloroCodon</h1>
            <p>
                A premium Streamlit workspace for GenBank codon-usage analysis, publication-ready figures,
                interactive result review, and complete CSV / Excel / ZIP export packages.
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )


def _render_feature_cards():
    st.markdown(
        """
        <div class="cc-feature-grid">
            <div class="cc-feature-card">
                <div class="cc-feature-icon">📊</div>
                <b>Complete genome summary</b>
                <span>A%, T%, G%, C%, GC positions, ENC, CAI, CBI, FOP and third-position statistics.</span>
            </div>
            <div class="cc-feature-card">
                <div class="cc-feature-icon">🧪</div>
                <b>Codon-bias analysis</b>
                <span>RSCU, optimal codons, high-frequency codons, amino-acid and stop-codon usage tables.</span>
            </div>
            <div class="cc-feature-card">
                <div class="cc-feature-icon">🎨</div>
                <b>Custom visual themes</b>
                <span>Separate colormap and theme controls for scatter plots, COA, neutrality and heatmaps.</span>
            </div>
            <div class="cc-feature-card">
                <div class="cc-feature-icon">📦</div>
                <b>Clean export package</b>
                <span>Download CSV tables, Excel workbook, method notes, summary text and figure files in one ZIP.</span>
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )



def _uploaded_file_to_temp_path(uploaded_file, work_dir):
    """Persist a Streamlit UploadedFile to a temporary GenBank path and return it."""
    suffix = Path(uploaded_file.name).suffix or ".gbk"
    safe_stem = sanitize_filename(Path(uploaded_file.name).stem, fallback="uploaded_genbank")
    out_path = Path(work_dir) / f"{safe_stem}{suffix}"
    out_path.write_bytes(uploaded_file.getvalue())
    return out_path


def _zip_directory_to_bytes(folder_path):
    """Zip a folder into memory for Streamlit download buttons."""
    folder_path = Path(folder_path)
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        for file_path in sorted(folder_path.rglob("*")):
            if file_path.is_file():
                zf.write(file_path, arcname=file_path.relative_to(folder_path))
    buffer.seek(0)
    return buffer.getvalue()


def _make_one_based_index_df(df, index_name="SL"):
    """Return a display-only DataFrame with a 1-based index for Streamlit tables."""
    if df is None:
        return df
    display_df = df.copy()
    display_df.index = range(1, len(display_df) + 1)
    display_df.index.name = index_name
    return display_df


def _display_dataframe_one_based(df, *, height=None, use_container_width=True, index_name="SL"):
    """Display every Streamlit table with 1-based indexing instead of Python's 0-based index."""
    display_df = _make_one_based_index_df(df, index_name=index_name)
    kwargs = {"use_container_width": use_container_width}
    if height is not None:
        kwargs["height"] = height
    st.dataframe(display_df, **kwargs)


def _dataframe_download_button(df, filename, label):
    if df is None or df.empty:
        return
    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button(label, data=csv, file_name=filename, mime="text/csv")


def _format_metric_value(value, decimals=3, integer=False):
    """Format values for premium metric cards."""
    if value is None:
        return "—"
    try:
        if pd.isna(value):
            return "—"
    except Exception:
        pass
    try:
        if integer:
            return f"{int(value):,}"
        if isinstance(value, (float, np.floating)):
            return f"{float(value):.{decimals}f}"
        if isinstance(value, (int, np.integer)):
            return f"{int(value):,}"
    except Exception:
        return str(value)
    return str(value)


def _metric_card(label, value, note=""):
    """Render a single premium metric card."""
    st.markdown(
        f"""
        <div class="cc-metric-card">
            <div class="cc-metric-label">{label}</div>
            <div class="cc-metric-value">{value}</div>
            <div class="cc-metric-note">{note}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def _render_metric_cards(stats):
    """Show the most important summary metrics in premium cards."""
    if not stats:
        st.info("No summary statistics available.")
        return

    metrics = [
        ("Genes", _format_metric_value(stats.get("Num_genes"), integer=True), "accepted CDS after duplicate removal"),
        ("Total codons", _format_metric_value(stats.get("Total_codons"), integer=True), "non-stop codons analysed"),
        ("GC%", _format_metric_value(stats.get("GC%")), "overall coding GC content"),
        ("Avg ENC", _format_metric_value(stats.get("Avg_ENC")), "effective number of codons"),
        ("CAI", _format_metric_value(stats.get("CAI")), "Sharp & Li adaptation index"),
        ("GC1%", _format_metric_value(stats.get("GC1%")), "first codon position"),
        ("GC2%", _format_metric_value(stats.get("GC2%")), "second codon position"),
        ("GC3%", _format_metric_value(stats.get("GC3%")), "third codon position"),
        ("CBI", _format_metric_value(stats.get("CBI")), "Wright codon bias index"),
        ("FOP", _format_metric_value(stats.get("FOP")), "frequency of optimal codons"),
    ]

    first_row = metrics[:5]
    second_row = metrics[5:]

    cols = st.columns(5)
    for col, (label, value, note) in zip(cols, first_row):
        with col:
            _metric_card(label, value, note)

    # Keep the two summary rows visually separated without injecting raw grid HTML.
    st.markdown('<div class="cc-metric-row-gap"></div>', unsafe_allow_html=True)

    cols = st.columns(5)
    for col, (label, value, note) in zip(cols, second_row):
        with col:
            _metric_card(label, value, note)



def _render_table_picker(tables, key_prefix="single"):
    """Interactive table browser for all output tables."""
    if not tables:
        st.info("No tables available.")
        return
    names = list(tables.keys())
    selected = st.selectbox("Select output table", names, key=f"{key_prefix}_table_select")
    df = tables.get(selected)
    if df is None:
        st.info("Selected table is unavailable.")
        return
    st.caption(f"Rows: {len(df):,} | Columns: {len(df.columns):,}")
    _display_dataframe_one_based(df, use_container_width=True, height=460)
    _dataframe_download_button(df, f"{selected}.csv", f"Download {selected}.csv")


def _render_plot_gallery(data, settings, key_prefix="single"):
    """Generate and display the six publication-ready plots."""
    figures = create_batch_figures(
        data,
        settings["fig_width"],
        settings["fig_height"],
        label_top=settings["label_top"],
        palette_name=settings["palette_name"],
        cmap_name=settings["cmap_name"],
        scatter_cmap_name=settings["scatter_cmap_name"],
        correlation_palette_name=settings["correlation_palette_name"],
        correlation_cmap_name=settings["correlation_cmap_name"],
        neutrality_palette_name=settings["neutrality_palette_name"],
        neutrality_cmap_name=settings["neutrality_cmap_name"],
        coa_palette_name=settings["coa_palette_name"],
        coa_cmap_name=settings["coa_cmap_name"],
    )
    plot_names = list(figures.keys())
    st.markdown("<div class='cc-panel cc-panel-tight'>", unsafe_allow_html=True)
    selected_plots = st.multiselect(
        "Plots to preview",
        plot_names,
        default=plot_names,
        key=f"{key_prefix}_plot_select",
        help="Choose one or more plots to display. Exported files still follow the format settings in the sidebar.",
    )
    st.markdown("</div>", unsafe_allow_html=True)
    for name in selected_plots:
        fig = figures.get(name)
        if fig is None:
            continue
        st.markdown(f"<div class='cc-panel'><div class='cc-section-title'>{name.replace('_', ' ')}</div>", unsafe_allow_html=True)
        st.pyplot(fig, clear_figure=False, use_container_width=True)
        st.markdown("</div>", unsafe_allow_html=True)
        try:
            plt.close(fig)
        except Exception:
            pass



def _make_single_output_zip(data, uploaded_name, settings):
    """Create the exact full output package for one file and return zip bytes."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        result_dir = Path(tmp_dir) / "ChloroCodon_Single_Result"
        result_dir.mkdir(parents=True, exist_ok=True)
        prefix = sanitize_filename(Path(uploaded_name).stem, fallback="ChloroCodon_result")
        save_analysis_package(
            data,
            output_dir=result_dir,
            prefix_name=prefix,
            formats=settings["formats"],
            fig_width=settings["fig_width"],
            fig_height=settings["fig_height"],
            label_top=settings["label_top"],
            palette_name=settings["palette_name"],
            cmap_name=settings["cmap_name"],
            scatter_cmap_name=settings["scatter_cmap_name"],
            correlation_palette_name=settings["correlation_palette_name"],
            correlation_cmap_name=settings["correlation_cmap_name"],
            neutrality_palette_name=settings["neutrality_palette_name"],
            neutrality_cmap_name=settings["neutrality_cmap_name"],
            coa_palette_name=settings["coa_palette_name"],
            coa_cmap_name=settings["coa_cmap_name"],
            save_figures=settings["save_figures"],
        )
        return _zip_directory_to_bytes(result_dir)


def _run_batch_analysis_streamlit(uploaded_files, settings):
    """Run the original batch workflow inside Streamlit and return logs plus zip bytes."""
    if len(uploaded_files or []) > BATCH_FILE_LIMIT:
        raise ValueError(f"Batch upload limit exceeded. Please upload a maximum of {BATCH_FILE_LIMIT} GenBank files.")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        input_dir = tmp_dir / "inputs"
        input_dir.mkdir(parents=True, exist_ok=True)
        batch_root = tmp_dir / f"ChloroCodon_Batch_Results_{timestamp}"
        batch_root.mkdir(parents=True, exist_ok=True)

        theme_sequence = list(BATCH_PLOT_SCHEMES)
        cmap_sequence = list(BATCH_HEATMAP_CMAPS)
        scatter_cmap_sequence = list(BATCH_SCATTER_CMAPS)
        correlation_cmap_sequence = list(CORRELATION_CMAPS)
        if settings.get("auto_cycle_batch_styles", True):
            random.shuffle(theme_sequence)
            random.shuffle(cmap_sequence)
            random.shuffle(scatter_cmap_sequence)
            random.shuffle(correlation_cmap_sequence)
        else:
            theme_sequence = [settings["palette_name"]]
            cmap_sequence = [settings["cmap_name"]]
            scatter_cmap_sequence = [settings["scatter_cmap_name"]]
            correlation_cmap_sequence = [settings["correlation_cmap_name"]]

        batch_log_rows = []
        combined = {
            "Batch_Genome_Stats": [],
            "Batch_RSCU": [],
            "Batch_Amino_Acid_Usage": [],
            "Batch_Stop_Codon_Usage": [],
            "Batch_Optimal_Codons": [],
            "Batch_QC_Log": [],
            "Batch_Duplicate_Log": [],
            "Batch_Error_Log": [],
        }

        progress = st.progress(0)
        status_box = st.empty()
        total = len(uploaded_files)

        for idx, uploaded_file in enumerate(uploaded_files, start=1):
            file_path = _uploaded_file_to_temp_path(uploaded_file, input_dir)
            file_name = Path(uploaded_file.name).name
            status_box.info(f"Processing {idx} of {total}: {file_name}")

            theme = theme_sequence[(idx - 1) % len(theme_sequence)]
            cmap = cmap_sequence[(idx - 1) % len(cmap_sequence)]
            scatter_cmap = scatter_cmap_sequence[(idx - 1) % len(scatter_cmap_sequence)]
            correlation_cmap = correlation_cmap_sequence[(idx - 1) % len(correlation_cmap_sequence)]
            neutrality_cmap = scatter_cmap_sequence[(idx + 2) % len(scatter_cmap_sequence)]
            coa_cmap = scatter_cmap_sequence[(idx + 5) % len(scatter_cmap_sequence)]
            if not settings.get("auto_cycle_batch_styles", True):
                neutrality_cmap = settings["neutrality_cmap_name"]
                coa_cmap = settings["coa_cmap_name"]

            try:
                data = analyze_genbank_file(file_path, min_codons=settings["min_codons"], table_id=settings["table_id"])
                folder = unique_folder(batch_root, Path(file_name).stem)
                prefix = folder.name
                tables = save_analysis_package(
                    data,
                    output_dir=folder,
                    prefix_name=prefix,
                    formats=settings["formats"],
                    fig_width=settings["fig_width"],
                    fig_height=settings["fig_height"],
                    label_top=settings["label_top"],
                    palette_name=theme,
                    cmap_name=cmap,
                    scatter_cmap_name=scatter_cmap,
                    correlation_palette_name=theme if settings.get("auto_cycle_batch_styles", True) else settings["correlation_palette_name"],
                    correlation_cmap_name=correlation_cmap,
                    neutrality_palette_name=theme if settings.get("auto_cycle_batch_styles", True) else settings["neutrality_palette_name"],
                    neutrality_cmap_name=neutrality_cmap if settings.get("auto_cycle_batch_styles", True) else settings["neutrality_cmap_name"],
                    coa_palette_name=theme if settings.get("auto_cycle_batch_styles", True) else settings["coa_palette_name"],
                    coa_cmap_name=coa_cmap if settings.get("auto_cycle_batch_styles", True) else settings["coa_cmap_name"],
                    save_figures=settings["save_figures"],
                )

                genes = data.get("genome_stats", {}).get("Num_genes", "")
                msg = "Saved complete output package"
                row = {
                    "File": file_name,
                    "Input_path": file_name,
                    "Status": "Done",
                    "Genes": genes,
                    "Plot_theme": theme,
                    "Heatmap_cmap": cmap,
                    "Scatter_cmap": scatter_cmap,
                    "Correlation_cmap": correlation_cmap,
                    "Neutrality_cmap": neutrality_cmap,
                    "COA_cmap": coa_cmap,
                    "Output_folder": str(folder.relative_to(batch_root)),
                    "Message": msg,
                }
                batch_log_rows.append(row)

                metadata = {
                    "File": file_name,
                    "Input_path": file_name,
                    "Output_folder": str(folder.relative_to(batch_root)),
                    "Plot_theme": theme,
                    "Heatmap_cmap": cmap,
                    "Scatter_cmap": scatter_cmap,
                    "Correlation_cmap": correlation_cmap,
                    "Neutrality_cmap": neutrality_cmap,
                    "COA_cmap": coa_cmap,
                }
                for sheet_key, table_name in [
                    ("Batch_Genome_Stats", "summary"),
                    ("Batch_RSCU", "rscu"),
                    ("Batch_Amino_Acid_Usage", "amino_acid_usage"),
                    ("Batch_Stop_Codon_Usage", "stop_codon_usage"),
                    ("Batch_Optimal_Codons", "optimal_codons"),
                    ("Batch_QC_Log", "qc"),
                    ("Batch_Duplicate_Log", "duplicate_removed"),
                ]:
                    df = tables.get(table_name)
                    if df is None:
                        continue
                    df2 = df.copy()
                    for key, value in reversed(list(metadata.items())):
                        df2.insert(0, key, value)
                    combined[sheet_key].append(df2)

            except Exception as exc:
                error_msg = str(exc)
                row = {
                    "File": file_name,
                    "Input_path": file_name,
                    "Status": "Failed",
                    "Genes": "",
                    "Plot_theme": theme,
                    "Heatmap_cmap": cmap,
                    "Scatter_cmap": scatter_cmap,
                    "Correlation_cmap": correlation_cmap,
                    "Neutrality_cmap": neutrality_cmap,
                    "COA_cmap": coa_cmap,
                    "Output_folder": "",
                    "Message": error_msg,
                }
                batch_log_rows.append(row)
                combined["Batch_Error_Log"].append(pd.DataFrame([row]))

            progress.progress(idx / total if total else 1.0)

        log_df = pd.DataFrame(batch_log_rows)
        log_df.to_csv(batch_root / "Batch_Log.csv", index=False)
        workbook_tables = {"Batch_Run_Log": log_df}
        for sheet, parts in combined.items():
            if parts:
                workbook_tables[sheet] = pd.concat(parts, ignore_index=True)
        save_excel_workbook_generic(workbook_tables, batch_root / "Batch_Summary.xlsx")

        done = int((log_df["Status"] == "Done").sum()) if not log_df.empty else 0
        failed = int((log_df["Status"] == "Failed").sum()) if not log_df.empty else 0
        status_box.success(f"Batch complete: {done} succeeded, {failed} failed.")
        zip_bytes = _zip_directory_to_bytes(batch_root)
        return log_df, zip_bytes


def _sidebar_settings():
    """Collect analysis/export settings in the Streamlit left sidebar."""
    st.sidebar.markdown(
        """
        <div class="cc-sidebar-title-card">
            <div class="cc-sidebar-kicker">Workspace</div>
            <div class="cc-sidebar-title">Workspace Analysis Controls</div>
            <div class="cc-sidebar-subtitle">Analysis settings · figure styling · export options</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    with st.sidebar.expander("Analysis settings", expanded=True):
        min_codons = st.number_input("Minimum CDS length (codons)", min_value=1, max_value=10000, value=30, step=1)
        table_id = st.number_input("NCBI genetic code table", min_value=1, max_value=33, value=11, step=1)

    with st.sidebar.expander("Figure size & labels", expanded=True):
        fig_width = st.number_input("Figure width", min_value=3.0, max_value=30.0, value=8.0, step=0.5)
        fig_height = st.number_input("Figure height", min_value=2.5, max_value=30.0, value=6.0, step=0.5)
        label_top = st.number_input("Label top outlier points", min_value=0, max_value=100, value=0, step=1)

    theme_options = list(PLOT_PALETTES.keys())
    with st.sidebar.expander("General plot styling", expanded=True):
        palette_name = st.selectbox(
            "General plot theme",
            theme_options,
            index=theme_options.index("default"),
            help="Controls the background, text, grid, and line style for general figures.",
        )
        cmap_name = st.selectbox(
            "RSCU heatmap colormap",
            BATCH_HEATMAP_CMAPS,
            index=BATCH_HEATMAP_CMAPS.index("viridis"),
        )
        scatter_cmap_name = st.selectbox(
            "ENC/PR2 scatter colormap",
            BATCH_SCATTER_CMAPS,
            index=BATCH_SCATTER_CMAPS.index("viridis"),
        )

    with st.sidebar.expander("Special plot color controls", expanded=False):
        correlation_palette_name = st.selectbox(
            "Correlations theme",
            theme_options,
            index=theme_options.index("default"),
            help="Controls the Correlations plot background, labels, and colorbar styling.",
        )
        correlation_cmap_name = st.selectbox(
            "Correlations colormap",
            CORRELATION_CMAPS,
            index=CORRELATION_CMAPS.index("coolwarm"),
        )
        neutrality_palette_name = st.selectbox(
            "Neutrality theme",
            theme_options,
            index=theme_options.index("default"),
            help="Controls the Neutrality plot background, labels, grid, reference lines, and regression line.",
        )
        neutrality_cmap_name = st.selectbox(
            "Neutrality point colormap",
            BATCH_SCATTER_CMAPS,
            index=BATCH_SCATTER_CMAPS.index("plasma"),
        )
        coa_palette_name = st.selectbox(
            "COA theme",
            theme_options,
            index=theme_options.index("default"),
            help="Controls the COA plot background, labels, grid, and reference lines.",
        )
        coa_cmap_name = st.selectbox(
            "COA point colormap",
            BATCH_SCATTER_CMAPS,
            index=BATCH_SCATTER_CMAPS.index("turbo"),
        )

    with st.sidebar.expander("Export settings", expanded=True):
        formats = st.multiselect("Figure export formats", ["png", "pdf", "svg", "tiff"], default=["png"])
        save_figures = st.checkbox("Save figures in output package", value=True)
        auto_cycle_batch_styles = st.checkbox("Batch mode: auto-cycle themes/colormaps", value=True)

    return {
        "min_codons": int(min_codons),
        "table_id": int(table_id),
        "fig_width": float(fig_width),
        "fig_height": float(fig_height),
        "label_top": int(label_top),
        "palette_name": palette_name,
        "cmap_name": cmap_name,
        "scatter_cmap_name": scatter_cmap_name,
        "correlation_palette_name": correlation_palette_name,
        "correlation_cmap_name": correlation_cmap_name,
        "neutrality_palette_name": neutrality_palette_name,
        "neutrality_cmap_name": neutrality_cmap_name,
        "coa_palette_name": coa_palette_name,
        "coa_cmap_name": coa_cmap_name,
        "formats": formats,
        "save_figures": bool(save_figures),
        "auto_cycle_batch_styles": bool(auto_cycle_batch_styles),
    }

def _uploaded_file_signature(uploaded_file):
    """Stable-enough signature for clearing stale Streamlit results when a new upload is selected."""
    if uploaded_file is None:
        return None
    return f"{getattr(uploaded_file, 'name', '')}|{getattr(uploaded_file, 'size', '')}"


def _uploaded_files_signature(uploaded_files):
    """Signature for a batch upload list."""
    if not uploaded_files:
        return None
    return "||".join(_uploaded_file_signature(f) or "" for f in uploaded_files)


def _analysis_settings_signature(settings):
    """Signature for settings that affect the biological analysis output."""
    return f"min_codons={int(settings.get('min_codons', 30))}|table_id={int(settings.get('table_id', 11))}"


def _single_result_signature(uploaded_file, settings):
    """Full single-mode signature: uploaded file + analysis settings."""
    file_sig = _uploaded_file_signature(uploaded_file)
    if file_sig is None:
        return None
    return f"{file_sig}|{_analysis_settings_signature(settings)}"


def _batch_result_signature(uploaded_files, settings):
    """Full batch-mode signature: uploaded files + analysis/export settings."""
    files_sig = _uploaded_files_signature(uploaded_files)
    if files_sig is None:
        return None
    return f"{files_sig}|{_analysis_settings_signature(settings)}|{_export_settings_signature(settings)}|auto_cycle_batch_styles={bool(settings.get('auto_cycle_batch_styles', True))}"


def _export_settings_signature(settings):
    """Signature for settings that affect exported figures/packages."""
    export_keys = [
        "fig_width",
        "fig_height",
        "label_top",
        "palette_name",
        "cmap_name",
        "scatter_cmap_name",
        "correlation_palette_name",
        "correlation_cmap_name",
        "neutrality_palette_name",
        "neutrality_cmap_name",
        "coa_palette_name",
        "coa_cmap_name",
        "formats",
        "save_figures",
    ]
    parts = []
    for key in export_keys:
        value = settings.get(key)
        if isinstance(value, list):
            value = ",".join(map(str, value))
        parts.append(f"{key}={value}")
    return "|".join(parts)


def _clear_single_result_state():
    """Clear only the stored single-file analysis outputs."""
    for key in [
        "single_result_data",
        "single_result_tables",
        "single_result_uploaded_name",
        "single_result_signature",
        "single_result_file_signature",
        "single_result_error",
        "single_result_zip_bytes",
        "single_result_zip_signature",
        "single_result_section",
        "single_table_select",
        "single_plot_select",
    ]:
        st.session_state.pop(key, None)


def _clear_batch_result_state():
    """Clear only the stored batch-analysis outputs."""
    for key in [
        "batch_result_log_df",
        "batch_result_zip_bytes",
        "batch_result_signature",
        "batch_result_error",
    ]:
        st.session_state.pop(key, None)

def _render_single_mode(settings):
    st.markdown(
        """
        <div class="cc-panel">
            <div class="cc-section-title">Single-file processing</div>
            <div class="cc-muted">Upload one GenBank file, run the full codon-usage workflow, preview all outputs, and download a complete ZIP package.</div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    uploaded_file = st.file_uploader(
        "Upload one GenBank file",
        type=["gb", "gbk", "genbank"],
        accept_multiple_files=False,
        key="single_file_uploader",
    )

    if not uploaded_file:
        _clear_single_result_state()
        st.info("Upload a `.gb`, `.gbk`, or `.genbank` file to begin.")
        return

    current_file_signature = _uploaded_file_signature(uploaded_file)
    current_signature = _single_result_signature(uploaded_file, settings)

    previous_file_signature = st.session_state.get("single_result_file_signature")
    if previous_file_signature is not None and previous_file_signature != current_file_signature:
        _clear_single_result_state()

    run = st.button("Run single-file analysis", type="primary", key="run_single")

    stored_data_exists = st.session_state.get("single_result_data") is not None
    previous_result_signature = st.session_state.get("single_result_signature")
    analysis_settings_changed = (
        stored_data_exists
        and previous_result_signature is not None
        and previous_result_signature != current_signature
    )

    if analysis_settings_changed:
        st.info("Analysis settings changed. Re-running the analysis with the updated CDS length/genetic-code table...")

    if run or analysis_settings_changed:
        with st.spinner("Running codon usage analysis..."):
            try:
                with tempfile.TemporaryDirectory() as tmp_dir:
                    file_path = _uploaded_file_to_temp_path(uploaded_file, tmp_dir)
                    data = analyze_genbank_file(
                        file_path,
                        min_codons=settings["min_codons"],
                        table_id=settings["table_id"],
                    )
                tables = build_output_tables_from_data(data)
                st.session_state["single_result_data"] = data
                st.session_state["single_result_tables"] = tables
                st.session_state["single_result_uploaded_name"] = uploaded_file.name
                st.session_state["single_result_signature"] = current_signature
                st.session_state["single_result_file_signature"] = current_file_signature
                st.session_state["single_result_error"] = None
                st.session_state.pop("single_result_zip_bytes", None)
                st.session_state.pop("single_result_zip_signature", None)
            except Exception as exc:
                _clear_single_result_state()
                st.session_state["single_result_error"] = str(exc)

    if st.session_state.get("single_result_error"):
        st.error(f"Analysis failed: {st.session_state['single_result_error']}")
        return

    data = st.session_state.get("single_result_data")
    tables = st.session_state.get("single_result_tables")
    result_uploaded_name = st.session_state.get("single_result_uploaded_name", uploaded_file.name)

    if data is None or tables is None:
        st.info("File selected. Click **Run single-file analysis** to generate the results.")
        return

    st.markdown(
        f"""
        <div class="cc-status-good">
            ✅ Analysis complete. Results remain available while you switch tables, sections, and plot options.<br>
            <span style="font-weight:500; color:#627267;">Minimum CDS length: {settings['min_codons']} codons · NCBI genetic code table: {settings['table_id']}</span>
        </div>
        """,
        unsafe_allow_html=True,
    )
    stats = data.get("genome_stats", {})
    _render_metric_cards(stats)

    section_options = ["Tables", "Plots", "Download package"]
    if st.session_state.get("single_result_section") not in section_options:
        st.session_state["single_result_section"] = "Tables"

    selected_section = st.radio(
        "Result section",
        section_options,
        horizontal=True,
        key="single_result_section",
    )

    if selected_section == "Tables":
        st.markdown("<div class='cc-panel'><div class='cc-section-title'>Output tables</div>", unsafe_allow_html=True)
        _render_table_picker(tables, key_prefix="single")
        st.markdown("</div>", unsafe_allow_html=True)

    elif selected_section == "Plots":
        _render_plot_gallery(data, settings, key_prefix="single")

    elif selected_section == "Download package":
        st.markdown(
            """
            <div class="cc-panel">
                <div class="cc-section-title">Complete output package</div>
                <div class="cc-muted">The ZIP contains CSV tables, the Excel workbook, summary and methods text files, plus the selected figure files generated by the export workflow.</div>
            </div>
            """,
            unsafe_allow_html=True,
        )
        zip_signature = f"{st.session_state.get('single_result_signature')}|{_export_settings_signature(settings)}"
        if st.session_state.get("single_result_zip_signature") != zip_signature:
            with st.spinner("Preparing downloadable ZIP package..."):
                st.session_state["single_result_zip_bytes"] = _make_single_output_zip(data, result_uploaded_name, settings)
                st.session_state["single_result_zip_signature"] = zip_signature

        st.download_button(
            "Download complete single-file output ZIP",
            data=st.session_state.get("single_result_zip_bytes", b""),
            file_name=f"{sanitize_filename(Path(result_uploaded_name).stem)}_ChloroCodon_results.zip",
            mime="application/zip",
        )


def _render_batch_mode(settings):
    st.markdown(
        """
        <div class="cc-panel">
            <div class="cc-section-title">Batch processing</div>
            <div class="cc-muted">Upload multiple GenBank files and generate separated result folders plus a combined batch summary workbook.</div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    uploaded_files = st.file_uploader(
        f"Upload multiple GenBank files — maximum {BATCH_FILE_LIMIT} files",
        type=["gb", "gbk", "genbank"],
        accept_multiple_files=True,
        key="batch_file_uploader",
        help=f"Batch mode accepts up to {BATCH_FILE_LIMIT} GenBank files per run.",
    )

    if not uploaded_files:
        _clear_batch_result_state()
        st.info(f"Upload two or more GenBank files for batch processing. Maximum allowed: {BATCH_FILE_LIMIT} files per run. A single file also works, but single-file mode gives previews.")
        return

    if len(uploaded_files) > BATCH_FILE_LIMIT:
        _clear_batch_result_state()
        st.error(f"Batch upload limit exceeded: {len(uploaded_files)} files selected. Please select a maximum of {BATCH_FILE_LIMIT} GenBank files.")
        preview_df = pd.DataFrame({"File": [f.name for f in uploaded_files[:BATCH_FILE_LIMIT]], "Size_bytes": [f.size for f in uploaded_files[:BATCH_FILE_LIMIT]]})
        st.caption(f"Showing the first {BATCH_FILE_LIMIT} selected files only. Remove extra files before running batch analysis.")
        _display_dataframe_one_based(preview_df, use_container_width=True, height=180)
        return

    current_signature = _batch_result_signature(uploaded_files, settings)
    previous_signature = st.session_state.get("batch_result_signature")
    if previous_signature is not None and previous_signature != current_signature:
        _clear_batch_result_state()
        st.info("Batch input or settings changed. Click **Run batch analysis** again to generate updated results.")

    st.markdown(f"<div class='cc-panel cc-panel-tight'><b>Selected files:</b> {len(uploaded_files)} / {BATCH_FILE_LIMIT}</div>", unsafe_allow_html=True)
    file_df = pd.DataFrame({"File": [f.name for f in uploaded_files], "Size_bytes": [f.size for f in uploaded_files]})
    _display_dataframe_one_based(file_df, use_container_width=True, height=180)

    run = st.button("Run batch analysis", type="primary", key="run_batch")

    if run:
        if settings["save_figures"] and not settings["formats"]:
            st.warning("No figure format is selected. Tables, text files, and Excel outputs will still be saved.")

        try:
            log_df, zip_bytes = _run_batch_analysis_streamlit(uploaded_files, settings)
            st.session_state["batch_result_log_df"] = log_df
            st.session_state["batch_result_zip_bytes"] = zip_bytes
            st.session_state["batch_result_signature"] = current_signature
            st.session_state["batch_result_error"] = None
        except Exception as exc:
            _clear_batch_result_state()
            st.session_state["batch_result_error"] = str(exc)

    if st.session_state.get("batch_result_error"):
        st.error(f"Batch analysis failed: {st.session_state['batch_result_error']}")
        return

    log_df = st.session_state.get("batch_result_log_df")
    zip_bytes = st.session_state.get("batch_result_zip_bytes")

    if log_df is None or zip_bytes is None:
        st.info("Files selected. Click **Run batch analysis** to generate the batch package.")
        return

    st.markdown("<div class='cc-panel'><div class='cc-section-title'>Batch run log</div>", unsafe_allow_html=True)
    _display_dataframe_one_based(log_df, use_container_width=True, height=360)
    st.markdown("</div>", unsafe_allow_html=True)
    st.download_button(
        "Download complete batch output ZIP",
        data=zip_bytes,
        file_name=f"ChloroCodon_Batch_Results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip",
        mime="application/zip",
    )


def main():
    init_codon_maps(11)
    _render_global_css()
    settings = _sidebar_settings()
    _render_premium_hero()

    st.markdown("<div class='cc-processing-title'>Processing mode</div>", unsafe_allow_html=True)
    mode = st.radio(
        "Processing mode",
        ["Single file", "Batch"],
        horizontal=True,
        index=0,
        key="processing_mode",
        label_visibility="collapsed",
    )

    if mode == "Single file":
        _render_single_mode(settings)
    else:
        _render_batch_mode(settings)

    st.markdown(f"<div class='cc-footer'>{DEVELOPER_FOOTER}</div>", unsafe_allow_html=True)

if __name__ == "__main__":
    main()
