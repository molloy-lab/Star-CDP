Star-CDP
=========
Star-CDP is a method for estimating cell lineage trees from CRISPR/Cas9-induced mutations under the Star Homoplasy parsimony criterion score ([Sashittal et al., 2023](https://doi.org/10.1016/j.cels.2023.11.005)). 
The input characters have 0 representing the unedited state, positive integers  representing edited states, and -1 representing the missing or ambiguous state. 
Under the Star Homoplasy assumption, only mutations between the unedited state (0) and edited states are allowed (note that this is equivalent to Camin-Sokal parsimony for binary characters, where 0 is the ancestral/unedited state and 1 is the only derived/edited state).
Our approach is novel in that it is guaranteed to return an optimal solution to the problem that obeys the set of constraints (clades aka subsets of cells) given as input.
In practice, the clade constraints are generated from candidate cell lineage trees recovered in prior analsyes or heuristic search, as described in [this example](example/README.md).

Usage
-----
To build, Star-CDP use commands:
```
git clone https://github.com/molloy-lab/Star-CDP.git
cd Star-CDP/src
make
```
Note: On Linux, we successfully compiled with gcc version 8.5.0. On Mac OS X, we successfully compiled with Apple clang version 15.0.0, which requires Apple command line tools to be installed. This can be done with the following commands
```
# sudo rm -rf /Library/Developer/CommandLineTools
xcode-select --install
```
and then following the instructions in the pop-up window.

Dependency
-----
1. C++ 17 or above
2. ASTRAL
3. Boost 1.80.0 or above

To install ASTRAL following the instructions
#In either case, before running Star-CDP, you must download ASTRAL and extract the zip folder into the src directory:
```
git clone https://github.com/smirarab/ASTRAL.git
mv ASTRAL tmp-ASTRAL
unzip tmp-ASTRAL/Astral.*.zip
rm -rf tmp-ASTRAL
```
To install boost
- On mac, it can be install via homebrew by the following instruction.
  ```
  brew install boost
  ```
  This installs up-to-date Boost(> 1.80.0) to ```/usr/local/```include and ```/usr/local/lib``` (or in /opt/homebrew on Apple Silicon).
- On Linux
  ```
  sudo apt update
  sudo apt install libboost-all-dev
  ```
  This installs Boost to ```/usr/local```

These defualt Boost path has been already set up in our [Makefile](https://github.com/molloy-lab/Star-CDP/blob/main/src/Makefile). However, if the Boost has been installed in a different path, the ```BOOST_INCLUDE_PATH``` and ```BOOST_LIB_PATH``` should be changed accordingly. 

Input
-----
Star-CDP requires the following three inputs files.
1. A file containing the character matrix, a comma-separated values (CSV) file that has rows representing cells and columns representing target sites. With same format as [Startle input character matrix](https://github.com/raphael-group/startle/blob/main/examples/n100_m30_d0.2_s0_p0.2_character_matrix.csv) Values of the character matrix must be either non-negative integers or '-1', with 0 indicating the unmutated state, other integers indicating mutated state, and '-1' as the missing data character.
2. A priors files containing the all mutations' probabilities, a comma-separated values (CSV) file with only three columns with the same format as [Startle Input priors csv file](https://github.com/raphael-group/startle/blob/main/examples/n100_m30_d0.2_s0_p0.2_mutation_prior.csv). The first column represents the site $x$, the second column represents a mutated state $y$ and the third column represents the probability that the mutation $0->y$ is on site x. 
3. A trees file containing trees of search space, a files containing lines of newick strings. This file could be constructed via heuristic search, refer to [this example](https://github.com/molloy-lab/Star-CDP/tree/main/example/3724_NT_All) for more details. 

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
./star-cdp -i <input character file>
           -m <input mutation probability file>
           -t <input trees from heuritic search or other sources>
           -g <label of outgroup (unedited / ancestor)>
           -o <output tree file>

USAGE for small parsimony problem, i.e., computing score for given tree:
./star-cdp -i <input character file> 
           -m <input mutation probability file>
           -q <input tree for scoring>

OPTIONS:
[-h|--help]
        Prints this help message.
(-i|--input) <input characters file>
        Name of file containing input characters in CSV format
(-m|--mutations) <mutations probability file>
        Name of file containing mutations probability
[(-t|--trees) <input trees file>]
        Name of file containing trees in newick format for constructing solution
        space with ASTRAL
(-x|-g|--outgroup) <outgroup or unedited cell label>
        Comma separated list of outgroup cells (e.g. unedited cell) used to root
        solution space
[(-q) <input species file>]
        Name of file containing tree for scoring in newick format
[(-o|--output) <output file>]
        Prefix of file for writing output tree (default: stdout)
[(-nosupp)]
        Turn off the calculation of clade support, which is the fraction of optimal
        solutions clade appears in
[(-consensus)]
        Compute greedy, majority, and strict consensus based on clade support
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
