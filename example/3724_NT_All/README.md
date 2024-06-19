Example for KPTracer on 3724_NT_All
--------------------------------------

To run Star-CDP on the data set `3724_NT_All`, we first need to go to the correct directory
```
cd Star-CDP/example/3724_NT_All/input
```

Then we need to perform heuristic search with [PAUP*](https://paup.phylosolutions.com). 

To run PAUP*, we need to create a binary version of character matrix, relabel the leaves to be suitable for PAUP*, and two add unedited cells, called `FAKEROOT` and `FAKEROOT2`.
```
python3 ../../tools/prepare_for_paup_and_starcdp.py \
    -i 3724_NT_All_pruned_character_matrix.csv \
    -w 3724_NT_All_priors.csv \
    -o 3724_NT_All_
```

This command produces the following files: 
* `3724_NT_All_pruned_character_matrix.csv-fakeroot` - version of character matrix with fake root
* `3724_NT_All_paup_binary.nex` - binary version of charater matrix
* `3724_NT_All_paup_leaf_map.csv` - leaf label name map
* `3724_NT_All_paup_camsok_hsearch_fast.nex` - commands for running heuristic search with PAUP*

Second, we download and execute PAUP* with the following commands:
```
./paup4a168_centos64 -n 3724_NT_All_paup_camsok_hsearch_fast.nex
```
if using Linux. **This can take hours; use our saved files**

This command conducts a heuristic search, saving the following files:
* `3724_NT_All_paup_all_saved.trees` - best 500 trees found
* `3724_NT_All_paup_all_score_saved.scores` - scores of best 500 trees found
* `3724_NT_All_paup_scon_high_score.tree` - strict consensus of high scoring trees found
* `3724_NT_All_paup_one_high_score.tree` - one of the high scoring trees found

However, these trees are not on the same label set as the original character matrix, so we need to relabel them.

```
cd ../output

python3 ../../tools/postprocess_from_paup.py \
    -i ../input/3724_NT_All_paup_scon_high_score.tree \
    -n ../input/3724_NT_All_paup_leaf_map.csv \
    -o rerun_paup_scon_high_score.tree

python3 ../../tools/postprocess_from_paup.py \
    -i ../input/3724_NT_All_paup_one_high_score.tree \
    -n ../input/3724_NT_All_paup_leaf_map.csv \
    -o rerun_paup_one_high_score.tree

python3 ../../tools/postprocess_from_paup_for_starcdp.py \
    -i ../input/3724_NT_All_paup_all_saved.trees \
    -n ../input/3724_NT_All_paup_leaf_map.csv \
    -o rerun_paup_all_saved_trees.trees
```

Now we can give these rooted trees as constraints to Star-CDP.
```
../../../src/star-cdp \
    -i ../input/3724_NT_All_pruned_character_matrix.csv-fakeroot \
    -m ../input/3724_NT_All_priors.csv \
    -t rerun_paup_all_saved.trees \
    -g FAKEROOT \
    -consensus \
    -o rerun_star_cdp
```

Lastly, we remove the `FAKEROOT` from the output of Star-CDP.
```
python3 ../../tools/postprocess_from_starcdp.py \
    -i rerun_star_cdp \
    -o nofakeroot
```

The parsimony score is listed in the output if you scroll up, but we can also check it using the following command:
```
../../../src/star-cdp \
    -q rerun_star_cdp_one_sol.tre-nofakeroot \
    -i ../input/3724_NT_All_pruned_character_matrix.csv \
    -m ../input/3724_NT_All_priors.csv | \
    grep "score"
```
returns `The star homoplasy  score: 4462.75`.

Now we compare this score to that of the tree computed with Startle-ILP tree; command
```
../../../src/star-cdp \
    -q startle_nni_python.tre \
    -i ../input/3724_NT_All_pruned_character_matrix.csv \
    -m ../input/3724_NT_All_priors.csv | \
    grep "score"
```
returns `The star homoplasy  score: 4990.52`.

For downstream analyses, we recommend using the **strict consensus tree**, stored in this file: `rerun_star_cdp_strict_consensus.tre-nofakeroot`.

For example, we can compare the Star-CDP-SCon tree to the Startle-ILP tree, *after contracting mutationless branches*, with the command:
```
python3 ../../tools/compare_two_rooted_trees_under_star.py \
    -t1 rerun_star_cdp_strict_consensus.tre-nofakeroot \
    -t2 startle_nni_python.tre \
    -c1 1 -c2 1 -ex1 -ex2 \
    -m ../input/3724_NT_All_pruned_character_matrix.csv 
```
which returns `1207,276,271,191,186,85,0.692029,0.686347,0.313653` for our analyses run on Linux.
* `1207` = the numbers indicate number of leaves/cells
* `276` = number of internal branches in t1
* `271` = number of internal branches in t2
* `191` = number of internal branches in t1 that are missing from t2
* `186` = number of internal branches in t2 that are missing from t1
* `86` = number of shared branches
* `0.692029` = 191 / 276
* `0.686347` = 186 / 271
* `0.313653` = 86 / 271

Importantly, this comparison assumes the input trees are rooted properly when identify branches with mutations.
We also included `-ex` flags, it will save the trees with mutationless branches contracted for later use. 

We can also compare Star-CDP-SCon to PAUP*-SCon, *after contracting mutationless branches*, with the command:
```
python3 ../../tools/compare_two_rooted_trees_under_star.py \
    -t1 rerun_star_cdp_strict_consensus.tre-nofakeroot \
    -t2 rerun_paup_scon_high_score.tree \
    -c1 1 -c2 1 -ex2 \
    -m ../input/3724_NT_All_pruned_character_matrix.csv 
```
which returns `1207,276,259,17,0,259,0.061594,0.000000,1.000000` for our analyses run on Linux.

We can also compare to the tree distributed with KPTracer (estimated with a special version of Cassiopeia), 
```
 python3 ../../tools/compare_two_rooted_trees_under_star.py \
     -t1 rerun_star_cdp_strict_consensus.tre-nofakeroot \
     -t2 ../download/3724_NT_All_tree.nwk  \
     -c1 1 -c2 1 -ex2 \
     -m ../input/3724_NT_All_pruned_character_matrix.csv 
```

Repeating this command for different trees shows us that Star-CDP-SCon is the most similar to KPTracer:
* Star-CDP-SCon vs. KPTracer `1207,276,283,112,119,164,0.405797,0.420495,0.579505`
* PAUP-SCon vs. KPTracer -   `1207,259,283,101,125,158,0.389961,0.441696,0.558304`
* Startle-NNI (Py) vs. KPTracer - `1207,271,283,196,208,75,0.723247,0.734982,0.265018` 
* Startle-NNI (C++) vs. KPTracer - `1207,235,283,185,233,50,0.787234,0.823322,0.176678`

**IMPORTANT** Note that PAUP* requires integer weights so we truncate at the second decimal place and multiple by 100. This means there can be differences in scoring and consensus by PAUP* and methods that use higher precision.

A copy of the PAUP* user manual is available [here](https://phylosolutions.com/paup-documentation/paupmanual.pdf); also see [https://rothlab.ucdavis.edu/genhelp/paupsearch.html](https://rothlab.ucdavis.edu/genhelp/paupsearch.html).
