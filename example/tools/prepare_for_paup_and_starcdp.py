"""
MIT License 
Copyright (c) 2024 Junyan Dai and Erin Molloy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
import argparse
import pandas
import numpy
import treeswift
import os
import sys


def binarize(df, priors=None):
    columns = list(df.columns)
    cell_list = df[columns[0]].values
    new_rows = {}
    new_cols = [columns[0]]
    for cell in cell_list:
        new_rows[cell] = [cell]
    weights = []

    for col in columns[1:]:
        character = df[col].values

        state_set = set(character)
        try:
            state_set.remove(0)
        except KeyError:
            pass
        try:
            state_set.remove(-1)
        except KeyError:
            pass

        states = sorted(list(state_set))
        states = numpy.array(states)
        nstates = states.size

        for state in states:
            xdf = priors[(priors["character"] == col) & \
                         (priors["state"] == state)]

            if xdf.shape[0] != 1:
                sys.exit("Bad mutation priors!\n")

            prob = xdf.probability.values[0]
            wght = -1.0 * numpy.log(prob)
            wght = int(numpy.round(numpy.round(wght, 3), 2) * 100)
            weights.append(str(wght))

        new_cols += [col + '_' + str(state) for state in states]

        for cell, x in zip(cell_list, character):
            if x == 0:
                new_rows[cell] += [0] * nstates
            elif x == -1:
                new_rows[cell] += [-1] * nstates
            else:
                data = [0] * nstates
                indx = numpy.where(states == x)[0][0]
                data[indx] = 1
                new_rows[cell] += data

    cols = new_cols
    rows = [row for row in new_rows.values()]
    return pandas.DataFrame(rows, columns=cols), weights


def write_to_nexus(df:pandas.DataFrame, weights:list, output:str):
    binary_file = output + "paup_binary.nex"
    leaf_map_file = output + "paup_leaf_map.csv"
    all_trees_file = output + "paup_all_saved.trees"
    all_scores_file = output + "paup_all_saved.scores"
    one_high_tree_file = output + "paup_one_high_score.tree"
    high_trees_file = output + "paup_high_score_saved.trees"
    high_scores_file = output + "paup_high_score_saved.scores"
    consensus_file = output + "paup_scon_high_score.tree"
    paup_file = output + "paup_camsok_hsearch_fast.nex"

    columns = list(df.columns)
    nchar = len(columns) - 1
    ntax = df.shape[0]

    with open(binary_file, 'w') as fp1, \
         open(leaf_map_file, 'w') as fp2:
        fp1.write("#NEXUS\n\n")
        fp1.write("Begin data;\n")
        fp1.write("\tDimensions ntax=%d nchar=%d;\n" % (ntax + 2, nchar))
        fp1.write("\tFormat datatype=standard gap=-;\n")
        fp1.write("\tMatrix\n")

        fake = ''.join(['0'] * nchar)
        fp1.write("FAKEROOT\n")
        fp1.write(fake + '\n')
        fp1.write("FAKEROOT2\n")
        fp1.write(fake + '\n')

        for index, row in df.iterrows():
            cell = row.values[0]
            new_label = str("LEAF%d" % index)
            fp1.write(new_label + '\n')
            fp2.write("%s,%s\n" % (cell, new_label))

            data = row.values[1:]
            data = [str(x) for x in data]
            data = ''.join(data)
            data = data.replace('-1', '-')
            fp1.write(data + '\n')

        fp1.write('\t;\n')
        fp1.write('End;\n')

    # Note - 
    # irreversible (Camin-Sokal); up means higher numbers are derived
    # you can use weight set of integers...
    # I can probably do round to second decimal and take base 100 weights
    with open(paup_file, 'w') as fp:
        fp.write("#NEXUS\n")
        fp.write("BEGIN PAUP;\n")
        fp.write("set maxtrees=510;\n")
        fp.write("set autoclose=yes warntree=no warnreset=no;\n")
        fp.write(f"execute {binary_file};\n")
        fp.write(f"outgroup FAKEROOT;\n")
        fp.write("typeset myctype = irrev.up:1-%d;\n" % nchar)
        fp.write("wtset mywtset vector = " + ' '.join(weights) + ";\n")
        fp.write("assume typeset=myctype wtset=mywtset;\n")
        fp.write("set criterion=parsimony;\n")
        fp.write("hsearch start=stepwise addSeq=random swap=None nreps=10 rseed=55555;\n")
        fp.write("hsearch start=1 swap=TBR nbest=500 rseed=12345;\n")
        fp.write("rootTrees;\n")
        fp.write("SortTrees;\n")
        fp.write(f"pscore all/scorefile={all_scores_file};\n")
        fp.write(f"savetrees File={all_trees_file} root=yes trees=all format=newick;\n")
        fp.write(f"savetrees File={one_high_tree_file} root=yes from=1 to=1 format=newick;\n")
        fp.write("filter best;\n")
        #fp.write(f"pscore all/scorefile={high_scores_file};\n")
        #fp.write(f"savetrees File={high_trees_file} root=yes trees=all format=newick;\n")
        fp.write(f"contree all/strict=yes treefile={consensus_file};\n")
        fp.write("END;\n")

def main(args):
    # Read inputs
    if args.priors is None:
        priors = None
    else:
        priors = pandas.read_csv(args.priors, sep=',')     

    cmat_df = pandas.read_csv(args.input, sep = ',')
    
    cmat_df.astype({col: int for col in cmat_df.columns[1:]})
    binary_df, weights = binarize(cmat_df, priors)

    write_to_nexus(binary_df, weights, args.output)

    # Add FAKEROOT to cmat_df
    columns = list(cmat_df.columns)
    nchar = len(columns) - 1
    rows = []
    row = {}
    row[columns[0]] = "FAKEROOT"
    for col in columns[1:]:
        row[col] = 0
    rows.append(row)
    new_df = pandas.DataFrame(rows, columns=columns)
    cmat_df = pandas.concat([cmat_df, new_df], ignore_index=True)
    cmat_df.to_csv(args.input + "-fakeroot", index=False)

    print("Finished writing paup nex file for " + args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str, default=None, required=True,
                        help="Input character matrix in TSV or CSV format")

    parser.add_argument("-w", "--priors", type=str, default=None, required=False,
                        help="Mutation priors in CSV format")

    parser.add_argument("-o", "--output", type=str, default=None, required=True,
                        help="Prefix for output files")

    main(parser.parse_args())
