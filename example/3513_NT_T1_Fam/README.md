Example for KPTracer on 3513_NT_T1_Fam
--------------------------------------

To run Star-CDP on the data set `3513_NT_T1_Fam`, we first need to go to the correct directory
```
cd Star-CDP/example/3513_NT_T1_Fam/input
```

Then we need to perform heuristic search with [PAUP*](https://paup.phylosolutions.com). 

To run PAUP*, we need to create a binary version of character matrix, relabel the leaves to be suitable for PAUP*, and add unedited cell, called `FAKEROOT`.
```
python3 ../../tools/prepare_for_paup_and_starcdp.py \
    -i 3513_NT_T1_Fam_pruned_character_matrix.csv \
    -w 3513_NT_T1_Fam_priors.csv \
    -o 3513_NT_T1_Fam_
```

This command produces the following files: 
* `3513_NT_T1_Fam_pruned_character_matrix.csv-fakeroot` - version of character matrix with fake root
* `3513_NT_T1_Fam_paup_binary.nex` - binary version of charater matrix
* `3513_NT_T1_Fam_paup_leaf_map.csv` - leaf label name map
* `3513_NT_T1_Fam_paup_camsok_hsearch_fast.nex` - commands for running heuristic search with PAUP*

Second, we download and execute PAUP* with the following commands:
```
./paup4a168_centos64 -n 3513_NT_T1_Fam_paup_camsok_hsearch_fast.nex
```
if using Linux.

This command conducts a heuristic search, saving the following files:
* `3513_NT_T1_Fam_paup_all_saved_trees.trees` - best 500 trees found
* `3513_NT_T1_Fam_paup_all_score_saved_trees.scores` - scores of best 500 trees found
* `3513_NT_T1_Fam_paup_high_saved_trees.trees` - high scoring trees found
* `3513_NT_T1_Fam_paup_high_score_saved_trees.scores` - scores of high scoring trees found
* `3513_NT_T1_Fam_paup_scon_high_score.tree` - strict consensus of high scoring trees found

However, these trees are not on the same label set as the original character matrix, so we need to relabel them.

```
cd ../output
python3 ../../tools/postprocess_from_paup.py \
    -i ../input/3513_NT_T1_Fam_paup_scon_high_score.tree \
    -n ../input/3513_NT_T1_Fam_paup_leaf_map.csv \
    -o rerun_paup_scon_high_score.tree
```

```
python3 ../../tools/postprocess_from_paup_for_starcdp.py \
    -i ../input/3513_NT_T1_Fam_paup_all_saved_trees.trees \
    -n ../input/3513_NT_T1_Fam_paup_leaf_map.csv \
    -o rerun_paup_all_saved_trees.trees
```

Now we can give these rooted trees as constraints to Star-CDP.
```
../../../src/star-cdp \
    -i ../input/3513_NT_T1_Fam_pruned_character_matrix.csv-fakeroot \
    -m ../input/3513_NT_T1_Fam_priors.csv \
    -t rerun_paup_all_saved_trees.trees \
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
    -i ../input/3513_NT_T1_Fam_pruned_character_matrix.csv \
    -m ../input/3513_NT_T1_Fam_priors.csv | \
    grep "score"
```
returns `The star homoplasy  score: 306.945`.

Now we compare this score to that of the tree computed with Startle-ILP tree; command
```
../../../src/star-cdp \
    -q startle_ilp.tre \
    -i ../input/3513_NT_T1_Fam_pruned_character_matrix.csv \
    -m ../input/3513_NT_T1_Fam_priors.csv | \
    grep "score"
```
returns `The star homoplasy  score: 315.507`.

For downstream analyses, we recommend using the **strict consensus tree**, stored in this file: `rerun_star_cdp_strict_consensus.tre-nofakeroot`.

For example, we can compare the Star-CDP-SCon tree to the Startle-ILP tree, *after contracting mutationless branches*, with the command:
```
python3 ../../tools/compare_two_rooted_trees_under_star.py \
    -t1 rerun_star_cdp_strict_consensus.tre-nofakeroot \
    -t2 startle_ilp.tre \
    -c1 1 -c2 1 -ex1 -ex2 \
    -m ../input/3513_NT_T1_Fam_pruned_character_matrix.csv 
```
which returns `86,13,15,8,10,5,0.615385,0.666667,0.333333` for our analyses run on Linux.
* 86 = the numbers indicate number of leaves/cells
* 13 = number of internal branches in t1
* 15 = number of internal branches in t2
* 8 = number of internal branches in t1 that are missing from t2
* 10 = number of internal branches in t2 that are missing from t1
* 5 = number of shared branches
* 0.615385 = 8 / 13
* 0.666667 = 10 / 15
* 0.333333 = 5 / 15

Importantly, this comparison assumes the input trees are rooted properly when identify branches with mutations.
We also included `-ex` flags, it will save the trees with mutationless branches contracted for later use. 

We can also compare Star-CDP-SCon to PAUP*-SCon, *after contracting mutationless branches*, with the command:
```
python3 ../../tools/compare_two_rooted_trees_under_star.py \
    -t1 rerun_star_cdp_strict_consensus.tre-nofakeroot \
    -t2 rerun_paup_scon_high_score.tree \
    -c1 1 -c2 1 -ex2 \
    -m ../input/3513_NT_T1_Fam_pruned_character_matrix.csv 
```
which returns `86,13,9,4,0,9,0.307692,0.000000,1.000000` for our analyses run on Linux.

A copy of the PAUP* user manual is available [here](https://phylosolutions.com/paup-documentation/paupmanual.pdf); also see [https://rothlab.ucdavis.edu/genhelp/paupsearch.html](https://rothlab.ucdavis.edu/genhelp/paupsearch.html).
