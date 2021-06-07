import warnings

from pdbparser import Parser
from Bio import PDB, BiopythonWarning
from sys import exit


COLOR_ERROR = "\033[31m"
COLOR_STRUCTURE = "\33[35m"
COLOR_RESIDUE = "\033[33m"
COLOR_ATOM = "\033[92m"
COLOR_END = "\033[0m"

DOWNLOAD_DIR = "pdb_files"

#warnings.simplefilter('ignore', BiopythonWarning)


def main():
    parser = None

    while not parser:
        protein_name = input("Enter protein indentifier >>> ")
        if protein_name:
            if protein_name == "exit":
                exit()
            else:
                # Create a parser object where the data are stored
                parser = Parser(protein_name)
        else:
            continue
    
    # Print the count of each member in the protein hierarchy
    print("-" * 50)
    print("\33[35m" + parser.structure.header["idcode"].center(50, " ") + "\33[0m")
    print("-" * 50)
    print("Structure width: \t\t%.3f\n" % parser.get_structure_width())
    print("Number of models:\t\t", parser.get_model_count())
    print("Number of chains:\t\t", parser.get_chain_count())
    print("Number of residues:\t\t", parser.get_residue_count())
    print("Number of atoms:\t\t", parser.get_atom_count())
    print("-" * 50, end="\n")

    #print(parser.GetHetResidues(0))

    # Stores a list of heterogen residues
    het_residues = parser.get_het_residues()
    
    # Radius and ligand input
    while True:
        try:
            residue_index = input("Ligand >>> ")
            if residue_index.replace(" ", "") == "exit":
                exit()
            residue = het_residues[int(residue_index) - 1]
            het_atoms = parser.get_het_atoms(residue)

            radius = float(input("Radius >>> "))

            if radius == "exit":
                exit()
            
            if radius < 0:
                parser.send_error("Radius must be positive")
                continue
        except ValueError:
            parser.send_error("Wrong radius type")
            continue
        except KeyboardInterrupt:
            exit()
        except:
            parser.send_error("No results...")
            continue
            
        # Searching object to search for the neighbors of an atom
        searcher = PDB.NeighborSearch([atom for atom in parser.structure.get_atoms()])

        for hetatm in het_atoms:
            # Store the coordinates a center position of the heterogen atom
            center = hetatm.get_vector().get_array()

            center_x = round(center[0], 3)
            center_y = round(center[1], 3)
            center_z = round(center[2], 3)

            # Store a list of all neighbors of an atom and a residue
            neighbor_atoms = [atom for atom in searcher.search(center, radius, level="A")]
            neighbor_residues = [residue for residue in searcher.search(center, radius, level="R")]

            # Print the neighbors
            print("-" * 50)
            print(f"{COLOR_ATOM}Atoms in radius {radius} with center coords X: {center_x} Y: {center_y} Z: {center_z}{COLOR_END}")
            print("-" * 50)
            for atm in neighbor_atoms:
                print(f"➢  Atom - {atm.get_id()}")

            print("-" * 50)
            print(f"{COLOR_RESIDUE}Residues in radius {radius} with center coords X: {center_x} Y: {center_y} Z: {center_z}{COLOR_END}")
            print("-" * 50)
            for res in neighbor_residues:
                print(f"➢  Residue - {res.resname}")
            
        print("-" * 50)


if __name__ == "__main__":
    main()