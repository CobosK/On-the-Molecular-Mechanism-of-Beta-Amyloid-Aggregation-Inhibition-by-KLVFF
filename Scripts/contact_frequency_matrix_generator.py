#!/usr/bin/env python3
"""
This script processes contact frequency data from molecular dynamics simulations.

It reads residue-level contact CSV files (starting with "contact_Inh-") containing contact information 
between peptide chains. For each input file, it creates an individual matrix with the number of frames 
in which each residue was involved in a contact (per chain). Then, it aggregates and normalizes all 
matrices across files to generate two summary files:

- `matriz-N-contacto.csv`: summed contact counts per residue and chain.
- `matriz-Normalizada-contacto.csv`: normalized frequencies per residue and chain (relative to total frames).

This allows for the analysis of contact hot spots across different chains and residues in the trajectory.

Input: folder with `contact_Inh-*.csv` files.
Output: one matrix per input + 2 global matrices.
"""

import os
import csv

# Input folder containing contact_Inh-*.csv files
folder = 'Input/Path/To/Contact_Files'

# Constants
TOTAL_FRAMES = 6000  # Total number of frames in the simulation
ELEMENTS = [
    "ACE10", "GLU11", "VAL12", "HIE13", "HIE14", "GLN15", "LYS16", "LEU17", "VAL18",
    "PHE19", "PHE20", "ALA21", "GLU22", "ASP23", "VAL24", "GLY25", "SER26", "ASN27",
    "LYS28", "GLY29", "ALA30", "ILE31", "ILE32", "GLY33", "LEU34", "MET35", "VAL36",
    "GLY37", "GLY38", "VAL39", "VAL40", "ILE41", "ALA42"
]
CADENAS = ["A", "B", "C", "D", "E", "F"]

def generar_archivo_salida(archivo_entrada, folder):
    """Generates an individual contact matrix file per input contact CSV."""
    try:
        directorio_entrada, nombre_archivo = os.path.split(archivo_entrada)
        nombre_archivo_sin_extension = os.path.splitext(nombre_archivo)[0]
        partes_nombre = nombre_archivo_sin_extension.split("_")

        if len(partes_nombre) < 2 or not nombre_archivo.startswith("contact_Inh"):
            print(f"File format not recognized: {archivo_entrada}")
            return

        cadena_residuo = partes_nombre[1]
        cadena, residuo_numero = cadena_residuo.split("-")
        nombre_salida = f"matriz-N-contacto-{cadena}_{residuo_numero}.csv"

        conteos_residuos = {c: {r: 0 for r in ELEMENTOS} for c in CADENAS}
        total_frames = 0

        with open(archivo_entrada, 'r') as archivo:
            lector_csv = csv.reader(archivo)
            next(lector_csv)
            for fila in lector_csv:
                residuos = [res.replace('"', '').replace("'", "") for res in fila[2].split(', ') if res]
                total_frames += 1
                for res in residuos:
                    letra_cadena = res[-1]
                    nombre_residuo = res[:-2]
                    if nombre_residuo in ELEMENTOS:
                        conteos_residuos[letra_cadena][nombre_residuo] += 1

        archivo_salida = os.path.join(folder, nombre_salida)
        with open(archivo_salida, 'w', newline='') as salida_csv:
            escritor = csv.writer(salida_csv)
            escritor.writerow(["Residue"] + CADENAS)
            for residuo in ELEMENTS:
                row = ["{:5d}".format(conteos_residuos[c][residuo]) for c in CADENAS]
                escritor.writerow([residuo] + row)

        print(f"Generated: {archivo_salida}")

    except FileNotFoundError:
        print(f"File not found: {archivo_entrada}")

def sumar_y_normalizar_archivos(folder):
    """Sums and normalizes all contact matrices in the folder."""
    sumas = {res: [0]*len(CADENAS) for res in ELEMENTS}

    for archivo in os.listdir(folder):
        if archivo.endswith(".csv") and archivo.startswith("matriz-N-contacto"):
            with open(os.path.join(folder, archivo), 'r') as f:
                lector = csv.reader(f)
                next(lector)
                for fila in lector:
                    if fila[0] in sumas:
                        for i, val in enumerate(fila[1:]):
                            sumas[fila[0]][i] += int(val)

    with open(os.path.join(folder, "matriz-N-contacto.csv"), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Residue"] + CADENAS)
        for res, vals in sumas.items():
            writer.writerow([res] + vals)

    print("Generated: matriz-N-contacto.csv")

    with open(os.path.join(folder, "matriz-Normalizada-contacto.csv"), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Residue"] + CADENAS)
        for res, vals in sumas.items():
            normalized = ["{:.6f}".format(val / TOTAL_FRAMES) for val in vals]
            writer.writerow([res] + normalized)

    print("Generated: matriz-Normalizada-contacto.csv")

def buscar_y_procesar_carpetas(folder):
    """Finds and processes all contact_Inh CSV files in the given folder."""
    for root, _, files in os.walk(folder):
        for archivo in files:
            if archivo.endswith(".csv") and archivo.startswith("contact_Inh"):
                generar_archivo_salida(os.path.join(root, archivo), folder)

    sumar_y_normalizar_archivos(folder)

buscar_y_procesar_carpetas(folder)
