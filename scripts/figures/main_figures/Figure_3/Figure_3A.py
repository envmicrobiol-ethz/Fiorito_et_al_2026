#!/usr/bin/env python3

# DESCRIPTION
# Generates Figure 3A: a stacked histogram of predicted heme-binding motif
# density among Tier1 credible MHC candidates, grouped by the seven most
# abundant phyla.
#
# Heme-binding motif density is calculated as:
#   100 × total heme-binding motifs / protein length
#
# INPUT
# CSV table containing at least:
#   tier
#   sum_heme_motifs_member
#   seq_length_x
#   Phylum
#
# OUTPUT
# Figure 3A as HTML, SVG and PDF.
#
# USAGE
# python3 Figure_3A.py \
#   MHC_overview.csv \
#   figure_output_directory

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate Figure 3A: heme-binding motif density among "
            "Tier1 credible MHC candidates."
        )
    )
    parser.add_argument(
        "input_csv",
        type=Path,
        help="Input MHC overview CSV.",
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Directory for Figure 3A outputs.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=7,
        help="Number of most abundant phyla displayed separately. Default: 7.",
    )
    parser.add_argument(
        "--min-total-motifs",
        type=int,
        default=3,
        help="Minimum total number of heme-binding motifs. Default: 3.",
    )
    parser.add_argument(
        "--density-threshold",
        type=float,
        default=3.0,
        help="Vertical reference line in motifs per 100 aa. Default: 3.0.",
    )
    parser.add_argument(
        "--nbins",
        type=int,
        default=60,
        help="Number of histogram bins. Default: 60.",
    )
    return parser.parse_args()


def coerce_first_number(series: pd.Series) -> pd.Series:
    """Extract the first numerical value from each entry."""
    text = (
        series.astype("string")
        .str.replace(",", ".", regex=False)
    )

    number = text.str.extract(
        r"([0-9]+(?:\.[0-9]+)?)",
        expand=False,
    )

    return pd.to_numeric(
        number,
        errors="coerce",
    )


def strip_phylum_prefix(value: object) -> str:
    """Remove the GTDB p__ prefix from a phylum label."""
    label = str(value)

    if label.startswith("p__"):
        return label[3:]

    return label


def main() -> None:
    args = parse_args()

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_html = (
        args.output_dir
        / "Figure_3A_heme_binding_motif_density.html"
    )
    output_svg = (
        args.output_dir
        / "Figure_3A_heme_binding_motif_density.svg"
    )
    output_pdf = (
        args.output_dir
        / "Figure_3A_heme_binding_motif_density.pdf"
    )

    df = pd.read_csv(
        args.input_csv,
        low_memory=False,
        dtype=str,
    )

    plot_df = df.loc[
        df["tier"].astype(str).eq("Tier1_credible")
    ].copy()

    plot_df["sum_heme_motifs_member"] = (
        coerce_first_number(
            plot_df["sum_heme_motifs_member"]
        )
    )

    plot_df["seq_length_x"] = coerce_first_number(
        plot_df["seq_length_x"]
    )

    plot_df = plot_df.dropna(
        subset=[
            "sum_heme_motifs_member",
            "seq_length_x",
            "Phylum",
        ]
    ).copy()

    plot_df = plot_df.loc[
        (plot_df["seq_length_x"] > 0)
        & (
            plot_df["sum_heme_motifs_member"]
            >= args.min_total_motifs
        )
    ].copy()

    plot_df["heme_per_100aa"] = (
        100.0
        * plot_df["sum_heme_motifs_member"]
        / plot_df["seq_length_x"]
    )

    plot_df["Phylum"] = (
        plot_df["Phylum"]
        .astype("string")
        .fillna("p__Unknown")
    )

    phylum_counts = plot_df["Phylum"].value_counts()

    top_phyla = (
        phylum_counts
        .head(args.top_n)
        .index
        .tolist()
    )

    plot_df["Phylum_plot"] = np.where(
        plot_df["Phylum"].isin(top_phyla),
        plot_df["Phylum"],
        "Other",
    )

    plot_df["Phylum_plot_clean"] = (
        plot_df["Phylum_plot"]
        .map(strip_phylum_prefix)
    )

    top_phyla_clean = sorted(
        strip_phylum_prefix(phylum)
        for phylum in top_phyla
    )

    category_order = top_phyla_clean + ["Other"]

    phylum_colors = {
        "p__Acidobacteriota": "#B2182B",
        "p__Actinomycetota": "#FFF04F",
        "p__Bacteroidota": "#FE8CA1",
        "p__Chloroflexota": "#FF8C00",
        "p__Desulfobacterota": "#C77526",
        "p__Desulfobacterota_G": "#DD1C77",
        "p__FCPU426": "#EE82EE",
        "p__Planctomycetota": "#852F88",
        "p__Pseudomonadota": "#8C510A",
        "p__Verrucomicrobiota": "#D8B365",
        "p__Thermoproteota": "#08306B",
        "p__Halobacteriota": "#2171B5",
    }

    fallback_colors = [
        "#4DAF4A",
        "#377EB8",
        "#984EA3",
        "#A65628",
        "#999999",
        "#66C2A5",
        "#E6AB02",
    ]

    color_map = {}
    fallback_index = 0

    for phylum in top_phyla:
        clean_name = strip_phylum_prefix(phylum)

        if phylum in phylum_colors:
            color_map[clean_name] = phylum_colors[phylum]
        else:
            color_map[clean_name] = fallback_colors[
                fallback_index % len(fallback_colors)
            ]
            fallback_index += 1

    color_map["Other"] = "#BDBDBD"

    hover_columns = [
        column
        for column in (
            "Member",
            "SeqC_ID",
            "sum_heme_motifs_member",
            "seq_length_x",
            "heme_per_100aa",
        )
        if column in plot_df.columns
    ]

    figure = px.histogram(
        plot_df,
        x="heme_per_100aa",
        color="Phylum_plot_clean",
        nbins=args.nbins,
        opacity=1.0,
        color_discrete_map=color_map,
        category_orders={
            "Phylum_plot_clean": category_order
        },
        hover_data=hover_columns,
    )

    figure.update_traces(
        marker_line_width=0.45,
        marker_line_color="black",
        showlegend=False,
    )

    for name in category_order:
        figure.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                name=name,
                marker={
                    "size": 12,
                    "color": color_map[name],
                    "symbol": "square",
                    "line": {"width": 0},
                },
                showlegend=True,
                hoverinfo="skip",
            )
        )

    figure.update_layout(
        template="simple_white",
        xaxis_title=(
            "Predicted heme-binding motifs per 100 aa"
        ),
        yaxis_title="Number of proteins",
        barmode="stack",
        bargap=0.03,
        width=1100,
        height=600,
        margin={
            "l": 80,
            "r": 220,
            "t": 40,
            "b": 70,
        },
        font={
            "family": "Arial",
            "size": 14,
        },
        legend_title_text="Phylum",
        legend={
            "orientation": "v",
            "yanchor": "top",
            "y": 1.0,
            "xanchor": "left",
            "x": 1.02,
            "font": {"size": 12},
            "itemsizing": "constant",
            "traceorder": "normal",
        },
    )

    figure.update_xaxes(
        showgrid=True,
        gridwidth=0.5,
        ticks="outside",
        ticklen=6,
        tickwidth=1,
        zeroline=False,
        showline=True,
        linewidth=1,
        linecolor="black",
    )

    figure.update_yaxes(
        showgrid=True,
        gridwidth=0.5,
        ticks="outside",
        ticklen=6,
        tickwidth=1,
        zeroline=False,
        showline=True,
        linewidth=1,
        linecolor="black",
    )

    figure.add_vline(
        x=args.density_threshold,
        line_width=1,
        line_dash="dot",
        line_color="black",
    )

    figure.write_html(
        output_html,
        include_plotlyjs="cdn",
    )

    try:
        figure.write_image(output_svg)
        figure.write_image(output_pdf)
    except Exception as error:
        print(
            "[WARNING] SVG/PDF export was skipped. "
            "Install a Plotly-compatible Kaleido version to export "
            f"static images. Details: {error}"
        )

    print(f"Figure 3A HTML: {output_html}")
    print(f"Figure 3A SVG: {output_svg}")
    print(f"Figure 3A PDF: {output_pdf}")


if __name__ == "__main__":
    main()
