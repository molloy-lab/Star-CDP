KPTracer Data Preparation
-------------------------

**Step 1:** Download data from [this Zenodo data repository](https://zenodo.org/records/5847462) and select the following files from `KPTracer-Data/trees`.

**3513_NT_T1_Fam:**
* `3513_NT_T1_Fam_character_matrix.txt`
* `3513_NT_T1_Fam_priors.pkl`
* `3513_NT_T1_Fam_tree.nwk`
* `neighbor_joining/3513_NT_T1_Fam_tree_nj.processed.tree`

**3724_NT_All:**
* `3724_NT_All_character_matrix.txt`
* `3724_NT_All_priors.pkl`
* `3724_NT_All_tree.nwk`
* `neighbor_joining/3724_NT_All_tree_nj.processed.tree`

**Step 2:** Unpickle the mutation priors with this [script](tools/pickle2csv_priors.py). All priors are saved [here](kptracer-priors.tar.gz) because it requires careful installation of Cassiopeia and other Python modules.

**Step 3:** Reformat the character matrix for Startle and Star-CDP.
```
sed 's/\t/,/g' download/3513_NT_T1_Fam_character_matrix.txt | \
    sed 's/cellBC//g' | \
    sed 's/,-/,-1/g' | \
    > input/3513_NT_T1_Fam_character_matrix.csv

sed 's/\t/,/g' download/3724_NT_All_character_matrix.txt | \
    sed 's/cellBC//g' | \
    sed 's/,-/,-1/g' | \
    > input/3724_NT_All_character_matrix.csv
```

**Step 4:** Prune the matrix with [this script](https://github.com/raphael-group/startle/blob/main/scripts/prune.py); see [Sashittal, Schmidt, et al. (2023)](https://doi.org/10.1016/j.cels.2023.11.005).

**Step 5:** Run methods.


