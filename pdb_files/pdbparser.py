from Bio import PDB, BiopythonWarning

import warnings
import sys
import os


# Constants
PDB_DIR = "pdb_files"

COLOR_ERROR = "\033[31m"
COLOR_STRUCTURE = "\33[35m"
COLOR_RESIDUE = "\033[33m"
COLOR_ATOM = "\033[92m"
COLOR_END = "\033[0m"

warnings.simplefilter('ignore', BiopythonWarning) # Filter out the irrelevant warnings 


class Parser(PDB.PDBParser):
    """The Parser for PDB files

        - protein: The protein name,
        - structure_path: The path to the pdb file,
        - structure: The structure object,
    """

    def __init__(self, protein, *args, **kwargs):
        super().__init__(args, kwargs)

        self.protein = protein
        self.structure_path = f"{PDB_DIR}/pdb{protein.lower()}.ent"

        # Load PDB file
        if os.path.exists(self.structure_path):
            self.load_pdb(self.structure_path)
        else:
            self.send_error("File not found!")

            try:
                self.download_pdb(protein_name=protein)
                self.structure = self.get_structure("x", self.structure_path)
            except:
                sys.exit()

    
    def get_model_count(self):
        """Return the count of models in the structure"""

        return len([model for model in self.structure.get_models()])
    
    def get_chain_count(self):
        """Return the count of chains in the structure"""

        return len([chain for chain in self.structure.get_chains()])
    
    def get_residue_count(self):
        """Return the count of residues in the structure"""

        return len([residue for residue in self.structure.get_residues()])
    
    def get_atom_count(self):
        """Return the count of atoms in the structure"""

        return len([atom for atom in self.structure.get_atoms()])
    
    def get_structure_width(self):
        """Return the width of the structure - maximum distance of any two atoms"""

        _max = 0
        atoms = self.structure.get_atoms()
        for atom in atoms:
            for _atom in atoms:
                if (atom - _atom) > _max:
                    _max = atom - _atom

        return _max
    
    def get_het_residues(self):
        """Return a list of all heterogen residues"""

        residues = self.structure[0].get_residues()
        het_residues = []
        
        for residue in residues:
            residue_id = residue.get_id()
            hetfield = residue_id[0]

            if hetfield[0] != " ":
                het_residues.append(residue)

        return het_residues
    
    def get_het_atoms(self, het_residues=None):
        """Return a list of all heterogen atoms"""

        het_residues = het_residues or self.get_het_residues()
        het_atoms = []

        for atom in het_residues.get_atoms():
            het_atoms.append(atom)
        
        return het_atoms
    
    def download_pdb(self, protein_name):
        """Download a PDB file to the dowload folder"""

        pdb_list = PDB.PDBList()
        pdb_list.retrieve_pdb_file(protein_name, file_format = 'pdb', pdir=PDB_DIR)
    
    def load_pdb(self, path):
        """Load the pdb file if it exists"""

        try:
            self.structure = self.get_structure("x", path)
        except:
            self.send_error("PDB failed to load")
    
    def send_error(self, message):
        """Send an error message to the console"""

        print(f"{COLOR_ERROR}[Error] {message}{COLOR_END}")