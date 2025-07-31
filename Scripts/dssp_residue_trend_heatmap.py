#!/usr/bin/env python3
"""
This script processes DSSP-derived data from multiple molecular dynamics replicas to calculate
per-residue secondary structure content averages and generate heatmaps.

It performs the following steps:
1. Loads DSSP summary data files for each system and replica.
2. Averages the secondary structure content across replicas.
3. Saves the averaged values to CSV files.
4. Plots and saves heatmaps of the average content per residue.

Input: DSSP summary `.dat` files (one per replica and system).
Output:
- CSV file with per-residue averages: `promedio_[system].csv`
- Heatmap image: `heatmap_[system].png`

Intended for structural analysis of residue-level trends across simulation replicas.
"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import warnings

warnings.filterwarnings("ignore")

# Input/output paths
input_folder = 'Input/Path/To/DSSP/Summaries'
output_folder = os.path.join(input_folder, 'promedios')
os.makedirs(output_folder, exist_ok=True)

# Systems and replicates
systems = ["S_KLVFF-NME", "S_ACE-KLVFF-NME"]
num_replicates = 3

# Residue names per monomer
residue_names = [
    "Ace", "Glu", "Val", "His", "His", "Gln", "Lys", "Leu", "Val", "Phe", "Fhe", "Ala",
    "Glu", "Asp", "Val", "Gly", "Ser", "Asn", "Lys", "Gly", "Ala", "Ile", "Ile", "Gly",
    "Leu", "Met", "Val", "Gly", "Gly", "Val", "Val", "Ile", "Ala"
]

# Custom color map
cmap = mcolors.LinearSegmentedColormap.from_list("custom_blues", [(1, 1, 1), "#005F8F"], N=256)

# Process each system
for system in systems:
    replica_files = [
        os.path.join(input_folder, f"sumout_{system}_{i+1}.dat") for i in range(num_replicates)
    ]
    replica_files = [f for f in replica_files if os.path.exists(f)]

    if len(replica_files) == num_replicates:
        try:
            dfs = [pd.read_csv(f, sep=None, engine="python") for f in replica_files]
            df_concat = pd.concat(dfs)

            avg = df_concat.mean(axis=0)
            avg.name = system
            avg.to_csv(os.path.join(output_folder, f"promedio_{system}.csv"))

            # Plot heatmap
            plt.figure(figsize=(12, 1.5))
            sns.heatmap(
                avg.values[np.newaxis, :],
                cmap=cmap,
                cbar_kws={"label": "Average"},
                xticklabels=residue_names[:len(avg)],
                yticklabels=[system]
            )
            plt.xticks(rotation=90)
            plt.title(f"Average per-residue content – {system}")
            plt.tight_layout()
            plt.savefig(os.path.join(output_folder, f"heatmap_{system}.png"), dpi=300)
            plt.close()

            print(f"Processed: {system}")

        except Exception as e:
            print(f"Error processing {system}: {e}")
    else:
        print(f"Missing replicas for {system}")

print("Processing completed.")
