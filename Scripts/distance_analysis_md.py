#!/usr/bin/env python3
"""
This script calculates geometric distances between specific residues
within and between monomers in a molecular dynamics PDB trajectory.
The calculations are based on alpha carbon (CA) coordinates.
Output is written in tabular `.dat` format for each PDB file.
"""

import os

# Configuration
OLIGOMER_LENGTH = 6      # Number of monomers
MONOMER_LENGTH = 33      # Number of residues per monomer
INPUT_FOLDER = 'Input/Path/To/PDB_Files'
OUTPUT_FOLDER = os.path.join(INPUT_FOLDER, "Distances")
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

pdb_files = [f for f in os.listdir(INPUT_FOLDER) if f.endswith(".pdb")]

print("Creating inter-monomer distance files:")

def process_distances(lines, model_indices, filename_suffix):
    distances = [[] for _ in range(OLIGOMER_LENGTH)]
    for i in range(len(model_indices) - 1):
        for j in range(OLIGOMER_LENGTH):
            idx1 = int(model_indices[i]) + 482 * j + j
            idx2 = int(model_indices[i]) + 482 * (j + 1) + (j + 1)
            x1, y1, z1 = map(float, [lines[idx1][30:38], lines[idx1][38:46], lines[idx1][46:54]])
            x2, y2, z2 = map(float, [lines[idx2][30:38], lines[idx2][38:46], lines[idx2][46:54]])
            dist = ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2) ** 0.5
            distances[j].append(dist)

    frames = list(range(1, len(distances[0]) + 1))
    filename = f"{os.path.splitext(os.path.basename(pdb_path))[0][:-5]}-{filename_suffix}.dat"
    with open(os.path.join(OUTPUT_FOLDER, filename), 'w') as f:
        f.write("Frames\t," + "\t,".join(f"M{i+1}-M{i+2}" for i in range(OLIGOMER_LENGTH - 1)) + "\n")
        for i in range(len(frames)):
            f.write("\t,".join([str(frames[i])] + [f"{distances[j][i]:.4f}" for j in range(OLIGOMER_LENGTH - 1)]) + "\n")

def calculate_inter_distances(pdb_path):
    try:
        with open(pdb_path, 'r') as f:
            lines = f.readlines()
        model_indices = [i for i, line in enumerate(lines) if "MODEL" in line]
        if model_indices:
            for offset, suffix, residue in [
                (132, 'DHC', 'Val18'),
                (241, 'DG', 'Gly25'),
                (312, 'DHS', 'Ile31'),
                (439, 'DCT', 'Val40')
            ]:
                indices = [idx + offset for idx in model_indices]
                process_distances(lines, indices, f"Inter_{suffix}_{residue}")
        else:
            print(f"No 'MODEL' lines found in {pdb_path}")
    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")

for pdb_file in pdb_files:
    pdb_path = os.path.join(INPUT_FOLDER, pdb_file)
    calculate_inter_distances(pdb_path)

print("Creating intra-monomer distance files:")

def process_intra_distances(lines, model_indices, filename_suffix, offset2):
    distances = [[] for _ in range(OLIGOMER_LENGTH)]
    for i in range(len(model_indices) - 1):
        for j in range(OLIGOMER_LENGTH):
            idx1 = int(model_indices[i]) + 482 * j + j
            idx2 = int(model_indices[i]) + offset2 + 482 * j + j
            x1, y1, z1 = map(float, [lines[idx1][30:38], lines[idx1][38:46], lines[idx1][46:54]])
            x2, y2, z2 = map(float, [lines[idx2][30:38], lines[idx2][38:46], lines[idx2][46:54]])
            dist = ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2) ** 0.5
            distances[j].append(dist)

    frames = list(range(1, len(distances[0]) + 1))
    filename = f"{os.path.splitext(os.path.basename(pdb_path))[0][:-5]}-{filename_suffix}.dat"
    with open(os.path.join(OUTPUT_FOLDER, filename), 'w') as f:
        f.write("Frames\t," + "\t,".join(f"M{i+1}" for i in range(OLIGOMER_LENGTH)) + "\n")
        for i in range(len(frames)):
            f.write("\t,".join([str(frames[i])] + [f"{distances[j][i]:.4f}" for j in range(OLIGOMER_LENGTH)]) + "\n")

def calculate_intra_distances(pdb_path):
    try:
        with open(pdb_path, 'r') as f:
            lines = f.readlines()
        model_indices = [i for i, line in enumerate(lines) if "MODEL" in line]
        if model_indices:
            for offset, suffix, offset2 in [
                (9, 'Glu10-Gly33', 341),
                (295, 'Gly29-Ala42', 179)
            ]:
                indices = [idx + offset for idx in model_indices]
                process_intra_distances(lines, indices, f"Intra_{suffix}", offset2)
        else:
            print(f"No 'MODEL' lines found in {pdb_path}")
    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")

for pdb_file in pdb_files:
    pdb_path = os.path.join(INPUT_FOLDER, pdb_file)
    calculate_intra_distances(pdb_path)
