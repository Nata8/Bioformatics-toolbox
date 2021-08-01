# Bioformatics-toolbox

Implementations of basic methods from the field of sequence and structural
bioinformatics for mff course with id **NDBI044**.

# Getting Started

## Requirements
- [Python 3.7.5](https://www.python.org/),
- [BioPython](https://biopython.org/), `pip install bio`
- [FreeSASA](https://freesasa.github.io/python/) `pip install freesasa`

## Run this repository in virtual environment (copy and paste)


**Clone whole repository from github: <br />**
- `git clone https://github.com/Nata8/Bioformatics-toolbox.git` <br />

**Navigate to the repository: <br />**
- `cd Bioformatics-toolbox/` <br />

**Active your virtual environment: <br />**
- `source env/bin/activate` <br />

**Run python scripts: <br />**
- `python3 (proccesspdb.py or structpropts.py)` <br />

**Deactivate the virtual environment: <br />**
- `deactivate`  <br />

____________________________________________________________________________________________________________________________________________________________

# ABOUT PROJECT

## First Java Part

First five assignments have been written in Java (repository [Bioinformatics_toolbox_2nd](https://github.com/Nata8/Bioinformatics_toolbox_2nd)). 
 
## Second Python Part

Next two tasks have been written in Python (repository [Bioinformatics_toolbox](https://github.com/Nata8/Bioformatics-toolbox)). Output displays at the terminal.
`pdbparser.py` is the basic parser used in both tasks. Data are downloaded from PDB database, stored in `pdb_files` directory.

**6. Processing PDB files (`processpbd.py`)**

First, the user have to enter four-digit PDB identifier. If the file is not already downloaded,
the PDB file will be downloaded and stored. The implementation provides following functions as an output: <br />
- proteine structure width (maximum of distance of any two atoms) <br />
- number of models, chains, residues and atoms in the structure <br />

Subsequently, the user should enter the serial number of a ligand and radius (distance from given ligand).
Ligands can be found in PDB file - lines starting with HETATM. The output is a list of atoms and residues
being in given distance from given ligand and coordinates of ligand.

**7. Computing structure-related properties (`structpropts.py`)**

With help of the previous PDB parser, the implementation provides following functions as an output:

- the diameter of the protein and the ratio of surface and buried amino acids
- a histogram of amino acids composition of buried and exposed amino acids
- a portion of polar amino acids in the core and on the surface of the protein

Ratio of the surface and buried amino acids is computed using *FreeSASA* Pythom module.
*FreeSASA* values are calculated and residues are divided into two groups based on these values.
Buried group = *FreeSASA* value must be lower than `0.2` (threshold is set according to first source written by Chen). 


Some of the assumptions made in these tasks are not entirely realistic.

# Sources
1. Chen, H. (2005). Prediction of solvent accessibility and sites of deleterious mutations from protein sequence. Nucleic Acids Research, 33(10), 3193–3199.
2. http://siret.ms.mff.cuni.cz/hoksza/courses/bioinformatics
3. http://www.cse.chalmers.se/~kemp/teaching/UMF018/2010-2011/sequence2.pdf
4. https://freesasa.github.io/python/
5. https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
