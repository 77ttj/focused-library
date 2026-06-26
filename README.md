# Muulti-dimensional similarity-driven methods for molecular library construction
# 🧬 Overview
**In this work, we developed a comprehensive protein–ligand co-crystal database by systematically mining the Protein Data Bank (PDB). And on that basis, we proposed a set of multi-dimensional similarity-driven strategies for molecular library construction. The resulting libraries exhibit high novelty, potency, synthetic accessibility, and favorable drug-like properties, while covering diverse chemical space. These target-specific libraries provide rich candidate pools for hit identification and lead optimization.**

### Key Features
* Automated Data Retrieval: Maps UniProt IDs to PDB structures and handles batch downloading.
* Quality Control: Filters protein structures by crystallographic resolution.
* Intelligent Pocket Extraction: Automatically identifies co-crystallized ligands and extracts 6.0 Å binding pockets.
* Interaction Analysis: Calculates key physical interactions (Hydrogen bonds, Hydrophobic effects, Pi-Pi stacking, Pi-Cation, Metal coordination).
* Multi-Strategy Similarity Search: Compares targets using sequence alignment, Morgan/MACCS pocket fingerprints, motif overlap, and high-dimensional vector embeddings.
* De Novo Generation: Breaks down known ligands into fragment libraries and randomly recombines them, filtering outputs strictly by Lipinski's Rule of Five.

# 🛠️ Installation

### Prerequisites

Ensure you have Python 3.8+ installed. Due to the heavy reliance on cheminformatics and bioinformatics libraries, using a virtual environment (like conda or venv) is highly recommended.

#### Setup

1. Clone the repository:
```
Bash
git clone https://github.com/77ttj/focused-library
cd focused-library
```
2. Install the required dependencies:

```
Bash
pip install -r requirements.txt
```
Note: The requirements.txt includes the full computational environment. The core functionalities heavily rely on RDKit, Biopython, Pandas, SciPy, and Scikit-Learn.

# 🚀 Pipeline Workflow

**The repository is organized into sequential scripts. Run them in order to complete the full workflow.**

### Phase 1: Data Acquisition & Preprocessing

* 1_pdbid_search_download_base_uniprot_*.py: Queries the PDB API to find the best PDB structures for given UniProt IDs and downloads them locally (categorized by species: homo, mus, rattus, other).`chembl31_single_protein_list_deduplicated.csv` is the deduplicated target list.
* 2_downloaded_id_add_to_file.py: Verifies the successfully downloaded files and updates the metadata CSV.
* 3_move_processed_pdb_to_new_fold.py: Cleans up specific file formats (e.g., removing _1.pdb suffixes) and organizes them into a clean processing directory.
* 4_pdb_selection_by_resolution.py: Parses PDB headers and filters out low-quality structures (e.g., keeping only those with resolution $\le$ 4.5 Å).

### Phase 2: Pocket & Interaction Extraction

* 5_pdb_selection_ligand_constraint_extraction_pocket_1031.py: The core extraction engine. It identifies valid small-molecule ligands, extracts them, defines the binding pocket within 6.0 Å, and calculates specific non-covalent interactions to identify "key amino acids." Converts ligands to SMILES and saves pocket structures.

### Phase 3: Target Similarity & Molecular Generation

**These scripts correspond to four distinct similarity search strategies for molecular generation.**

* 6_compute_362base_pdb_sim.py: Treats entire binding pockets as molecules to calculate 3D Morgan and MACCS fingerprints, computing structural similarities against a baseline set.
* 7_generation_by_seq_similarity*.py: Uses local sequence alignment of the pocket amino acids to find similar targets. Extracts ligands, breaks them via RECAP, and generates up to 10,000 Lipinski-compliant novel molecules.
* 8_generation_by_keyaas_overlap*.py: Finds similar pockets based on the intersection/overlap ratio of critical interacting amino acids (or manually specified motifs like ['TYR', 'TYR', 'TYR']).
* 9_generation_by_pocket_similarity.py: Generates molecules based strictly on the 3D physicochemical fingerprint similarity (Dice similarity) of the binding pockets.
* 10_generation_by_compared_pocket_similarity.py: Embeds pockets into a high-dimensional vector space based on 362 reference structures and computes cosine similarity for downstream fragment generation.

# 🧪 Usage

### Reproduction

1. Obtain `uniprot_pdb_list_all_protein.csv`, which contains the corresponding PDB IDs of each target in the PDB database, and download the corresponding PDB files
```
run 1_pdbid_search_download_base_uniprot_all.py
```

```
PDB_all
|___O43451  # UniProt_id
|   |___2qly.pdb # pdb_id
|   |___...

|___O76074  # UniProt_id
|   |___1rkp.pdb # pdb_id
|   |___...
|
```
2. Write the successfully downloaded PDB entries to uniprot_pdb_list_protein_downloaded_new.csv.

```
run 2_downloaded_id_add_to_file.py
```
3. Cleans up specific file formats (e.g., removing _1.pdb suffixes) and organizes them into a clean processing directory.

```
run 3_move_processed_pdb_to_new_fold.py
```
4. Select high-resolution PDB structures(e.g., keeping only those with resolution $\le$ 4.5 Å).

```
run 4_pdb_selection_by_resolution.py
```
5. Identify valid small-molecule ligands, extracts them, defines the binding pocket within 6.0 Å, and calculates specific non-covalent interactions to identify "key amino acids."
```
run 5_pdb_selection_ligand_constraint_extraction_pocket_1031.py
```
After execution, the file `output_PDB_test_processed_select_ligandsplit_4_5_4.csv` is generated, which records pocket information, ligand information and key amino acid residues for each PDB structure.
The resulting file structure is as follows:
```
PDB_processed_select_ligandsplit_4_5_3
|___O43451  # UniProt_id
|   |___2qmj.pdb # pdb_id
|   |___2qmj_ligand_AC1.pdb
|   |___2qmj_ligand_GLC.pdb
|   |___2qmj_pocket_AC1.pdb
|   |___2qmj_pocket_GLC.pdb
|   |___...

|___O76074  # UniProt_id
|   |___1rkp.pdb # pdb_id
|   |___1rkp_ligand_IBM.pdb
|   |___1rkp_pocket_IBM.pdb
|
```
6. Computing structural similarities against a baseline set.
```
run 6_compute_362base_pdb_sim.py
```
After execution, the file `output_PDB_test_processed_select_ligandsplit_4_5_4_compute_362base_sim_refined.csv` is generated.
The file contains similarity scores between each PDB and the PDBs in the reference dataset.

7. Molecular Generation
```
run *_generation_by_*.py
```
* seq_similarity: Similarity Evaluation Based on Protein Sequences.
* keyaas_overlap_similarity: Similarity Evaluation Based on Key Amino Acid Overlap.
* pocket_similarity: Protein Pocket Similarity Evaluation.
* compared_pocket_similarity: Comparative Similarity Evaluation of Protein Pockets.

Running the four similarity strategy scripts separately will generate three types of files: 
`fragment_library_file_by_*_similarity.csv`, `fragments_file_by_*_similarity.csv`, and `new_molecules_file_by_*_similarity_filter_10000.csv`. 
These files store the ligands extracted from the top 1000 PDB structures ranked by similarity scores, the chemically reasonable fragments decomposed from these ligands via the RECAP (Retrosynthetic Combinatorial Analysis Procedure) algorithm, and the newly generated molecules, respectively.

