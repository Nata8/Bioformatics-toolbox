import warnings
import freesasa

from pdbparser import Parser
from Bio import PDB
from sys import exit


# Constants
COLOR_ERROR = "\033[31m"
COLOR_BLUE = "\33[35m"
COLOR_GREEN = "\33[92m"
COLOR_ORANGE = "\33[33m"
COLOR_END = "\033[0m"

PDB_DIR = "pdb_files"

# List of apolar acids
APOLAR_ACIDS = [
    "GLY",
    "ALA",
    "VAL",
    "LEU",
    "ILE",
    "MET",
    "PRO",
    "PHE",
    "TRP",
]

AMINO_ACIDS = [
    "ARG",
    "LYS",
    "ASP",
    "GLU",
    "GLN",
    "ASN",
    "HIS",
    "SER",
    "THR",
    "TYR",
    "CYS",
]

def gcd(x, y):
    """Return the greatest common divisor"""

    while y != 0:
        (x, y) = (y, x % y)

    return x


def main():
    structure = None
    parser = None

    # Input
    while not structure or not parser:
        try:
            protein_name = input("Enter protein indentifier >>> ")
        except KeyboardInterrupt:
            exit()

        if protein_name:
            parser = Parser(protein_name)
            structure = freesasa.Structure(parser.structure_path)
        else:
            continue

    # Create the freesasa Result object and extracts the residues
    result = freesasa.calc(structure)
    amino_residues_dict = {residue_id: residue for residue_id, residue in dict(list(result.residueAreas().values())[0]).items()}
    amino_residues = [residue for residue in amino_residues_dict.values()]

    # Create a list for surface residues, threshold selected according to Chen, see sources
    surface_residue_list = []
    for amino_residue in amino_residues:
        if amino_residue.relativeTotal > 0.2:
            surface_residue_list.append(amino_residue)

    # Create a list for burried residues
    burried_residue_list = []
    for amino_residue in amino_residues:
        if amino_residue not in surface_residue_list:
            burried_residue_list.append(amino_residue)

# HISTOGRAM
    
    # Lists of both amino acid categories, surface and burried
    surface_residue_types = [residue.residueType for residue in surface_residue_list]
    burried_residue_types = [residue.residueType for residue in burried_residue_list]
    
    # Create a dictionary for histogram
    # Key -> Amino Acid Name Y; Value -> Count
    amino_acids_composition = {}
    amino_acids_composition["surface"] = {}
    amino_acids_composition["core"] = {}
    for amino_acid in (APOLAR_ACIDS + AMINO_ACIDS):
        amino_acids_composition["surface"][amino_acid] = surface_residue_types.count(amino_acid)
        amino_acids_composition["core"][amino_acid] = burried_residue_types.count(amino_acid)

    polar_burried_residues = [residue for residue in burried_residue_list if residue.residueType not in APOLAR_ACIDS]
    polar_surface_residues = [residue for residue in surface_residue_list if residue.residueType not in APOLAR_ACIDS]

    # The counts of individual residue list
    amino_residues_count = len(amino_residues)
    burried_residues_count = len(burried_residue_list)
    surface_residues_count = len(surface_residue_list)

    polar_burried_residues_count = len(polar_burried_residues)
    apolar_burried_residues_count = burried_residues_count - polar_burried_residues_count

    polar_surface_residues_count = len(polar_surface_residues)
    apolar_surface_residues_count = surface_residues_count - polar_surface_residues_count

    # The divisor for computing the ratio of surface and buried amino acids
    divisor = gcd(burried_residues_count, surface_residues_count)

# OUTPUT

    print("-" * 50)
    print(COLOR_GREEN + protein_name.upper().center(50, " ") + COLOR_END)
    print("-" * 50)
    print("➢  Diameter:\t\t\t", COLOR_BLUE, round(parser.get_structure_width(), 2), COLOR_END)
    print("➢  All amino acids:\t\t", COLOR_BLUE, amino_residues_count, COLOR_END)
    print("\n➢  Surface amino acids:\t\t", COLOR_BLUE, surface_residues_count, COLOR_END)
    print("   ○ Polar:\t\t\t", COLOR_BLUE, polar_surface_residues_count, COLOR_END)
    print("   ○ Apolar:\t\t\t", COLOR_BLUE, apolar_surface_residues_count, COLOR_END)
    print("\n➢  Burried amino acids:\t\t", COLOR_BLUE, burried_residues_count, COLOR_END)
    print("   ○ Polar:\t\t\t", COLOR_BLUE, polar_burried_residues_count, COLOR_END)
    print("   ○ Apolar:\t\t\t", COLOR_BLUE, apolar_burried_residues_count, COLOR_END)
    print("-" * 50)
    print(COLOR_GREEN + "Ratio of surface and buried amino acids:".center(50, " ") + COLOR_END)
    print("-" * 50)
    print("➢ ", int(surface_residues_count / divisor), ":", int(burried_residues_count / divisor))
    print("-" * 50)
    print(COLOR_GREEN + "Histogram of amino acids composition:".center(50, " ") + COLOR_END)
    print("-" * 50)

    # Print out the histogram of amino acids composition
    for category in amino_acids_composition:
        print(COLOR_ORANGE, f"\n{category.upper()}", COLOR_END)
        print("_" * 20)
        for residue_type, residue_count in amino_acids_composition[category].items():
            print(f"│ {residue_type} │", ("█" * residue_count))

    print("-" * 50)
    main()

if __name__ == "__main__":
    main()
