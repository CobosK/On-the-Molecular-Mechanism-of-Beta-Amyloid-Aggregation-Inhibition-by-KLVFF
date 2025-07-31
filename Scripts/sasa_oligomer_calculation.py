#!/usr/bin/env python3
"""
This script calculates the Solvent Accessible Surface Area (SASA) of an oligomer
from AMBER molecular dynamics trajectories using MDAnalysis and FreeSASA.

For each simulation replicate, it:
1. Loads topology and trajectory files.
2. Selects only the oligomer (residues 1–198).
3. Computes the SASA for each frame using FreeSASA (Lee-Richards algorithm).
4. Saves the results in a tab-separated file with one row per frame.

Input: AMBER `.prmtop` and `.nc` files (three replicas).
Output: `*_SASA_oligomero.txt` per replica.

This is useful for analyzing the solvation profile of protein cores during dynamics.
"""

import warnings
import MDAnalysis as mda
import freesasa
import numpy as np
import tempfile
import os
import pandas as pd
from multiprocessing import Pool, cpu_count

# Simulation identifier (used in paths)
NAME = "S"

# FreeSASA probe radius (in Ångströms)
PROBE_RADIUS = 1.4
freesasa.setVerbosity(freesasa.silent)
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")

# FreeSASA configuration
sasa_params = freesasa.Parameters()
sasa_params.setAlgorithm(freesasa.LeeRichards)

# Simulation paths
base_path = f"Input/Path/To/Simulations/{NAME}"
simulations = [
    {"topology": f"{base_path}/{NAME}_1/{NAME}_1.prmtop", "trajectory": f"{base_path}/{NAME}_1/{NAME}_1_prod.nc", "id": "1"},
    {"topology": f"{base_path}/{NAME}_2/{NAME}_2.prmtop", "trajectory": f"{base_path}/{NAME}_2/{NAME}_2_prod.nc", "id": "2"},
    {"topology": f"{base_path}/{NAME}_3/{NAME}_3.prmtop", "trajectory": f"{base_path}/{NAME}_3/{NAME}_3_prod.nc", "id": "3"}
]

# Output directory
output_folder = os.path.join(base_path, 'SASA')
os.makedirs(output_folder, exist_ok=True)

def calculate_sasa(atomgroup):
    if len(atomgroup) == 0:
        return 0.0
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        atomgroup.write(tmp.name)
        structure = freesasa.Structure(tmp.name)
        result = freesasa.calc(structure, sasa_params)
        os.unlink(tmp.name)
    return round(result.totalArea(), 2)

def get_sequence(ag):
    return "-".join([res.resname for res in ag.residues])

def process_frame(frame, universe, oligomer):
    universe.trajectory[frame]
    sasa_value = calculate_sasa(oligomer)
    return [universe.trajectory.frame, sasa_value]

def process_simulation(sim):
    u = mda.Universe(sim['topology'], sim['trajectory'])
    oligomer = u.select_atoms("resid 1-198")
    sequence = get_sequence(oligomer)
    print(f"\n🔬 Oligomer sequence in simulation {sim['id']}:")
    print(sequence)
    print("=" * 80)

    output_name = os.path.splitext(os.path.basename(sim['topology']))[0]
    output_path = os.path.join(output_folder, f"{output_name}_SASA_oligomero.txt")

    frames = range(1, len(u.trajectory), 1)
    with Pool(processes=max(1, cpu_count() - 1)) as pool:
        results = pool.starmap(process_frame, [(f, u, oligomer) for f in frames])

    df = pd.DataFrame(results, columns=["Frame", "SASA_Oligomer"])
    df.to_csv(output_path, sep='\t', index=False, float_format="%.2f")
    print(f"✅ SASA calculated and saved: {output_path}")

# Main execution
for sim in simulations:
    process_simulation(sim)
