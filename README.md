Star-CDP
=========
Star-CDP is a method for estimating cell lineage trees from CRISPR/Cas9-induced mutations under the Star Homoplasy parsimony criterion score ([Sashittal et al., 2023](https://doi.org/10.1016/j.cels.2023.11.005)). 
The input characters have 0 representing the unedited state, positive integers  representing edited states, and -1 representing the missing or ambiguous state. 
Under the Star Homoplasy assumption, only mutations between the united state (0) and edited states are allowed.
Our approach is novel in that it is guaranteed to return an optimal solution to the problem that obeys the set of constraints (clades) given as input.
These constraints can be generated from trees recovered in prior analsyes or heuristic search, as described in [this example](example/README.md).

Usage
-----
To build, Dollo-CDP use commands:
```
git clone https://github.com/molloy-lab/Star-CDP.git
cd Star-CDP/src
make
```
Note: On Linux, we successfully compiled with gcc version 8.5.0 and version 9.3.0. On Mac OS X, we successfully compiled with Apple clang version 15.0.3; we were unable to compile with gcc installed via homebrew, unfortunately. The former requires Apple command line tools to be installed. This can be done with the following commands
```
# sudo rm -rf /Library/Developer/CommandLineTools
xcode-select --install
```
and then following the pop-up.

Alternatively, you could download binaries in a release. In either case, before running Star-CDP, you must download ASTRAL and extract the zip folder into the src directory:
```
git clone https://github.com/smirarab/ASTRAL.git
mv ASTRAL tmp-ASTRAL
unzip tmp-ASTRAL/Astral.*.zip
rm -rf tmp-ASTRAL
```

To run Star-CDP, we recommend working through [this example](example/README.md).

The usage options can be viewed with this command:
```
./star-cdp -h
```
The output should be
```
Star-CDP version 1.0.0
COMMAND: ./star-cdp 
===================================== Star-CDP =====================================
Star-CDP is a program that solves the large Star Homoplasy parsimony within a
clade-constrained version of tree space.

USAGE for large parsimony problem:
./star-cdp -i <input characters file>
           -m <mutations probability file>
           -t <trees from heuritic search or other sources>
           -g <outgroup or unedited cell label>
           -o <output file>

USAGE for small parsimony problem, i.e., computing score for given tree:
./star-cdp -i <input characters file> -q <input species tree>

OPTIONS:
[-h|--help]
        Prints this help message.
(-i|--input) <input characters file>
        Name of file containing input characters in CSV format
(-x|-g|--outgroup) <outgroup or unedited cell label>
        Comma separated list of outgroup cells (e.g. unedited cell) used to root
        solution space
(-m|--mutations) <mutations probability file>
        Name of file containing mutations probability
[(-t|--trees) <input trees file>]
        Name of file containing trees in newick format for constructing solution
        space with ASTRAL
[(-q) <input species file>]
        Name of file containing species trees in newick format
[(-o|--output) <output file>]
        Prefix of file for writing output species tree (default: stdout)
[(-nosupp)]
        Turn off the calculation of clade support, which is the fraction of optimal
        solutions clade appears in
[(-e|--equal)]
        Using equal weight for all mutations
[(-memory)]
        Amount of memory given to ASTRAL.(Defualt:16000M)

Contact: Post issue to Github (https://github.com/molloy-lab/Star-CDP)
         or email Junyan Dai (jdai1234@umd.edu) & Erin Molloy (ekmolloy@umd.edu)

If you use Star-CDP in your work, please cite:
  Dai and Molloy, 2024, Star-CDP, https://github.com/molloy-lab/Star-CDP/
====================================================================================
```
