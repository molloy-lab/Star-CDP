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


class CharacterMatrix:
    def __init__(self, inputf):
        self.inputf = inputf
        self.is_csv = self.find_character_matrix_format()
        self.is_miss_dash = self.find_missing_symbol()
        self.columns = self.find_column_names()
        self.cell2data = None
        self.eqclass = None
        self.eqclass_upto_miss = None
        self.df = None
        self.df_pruned = None       # Prune cells with same sequence data
        self.df_pruned_upto_miss = None  # Also allow for missing data

    def find_character_matrix_format(self):
        with open(self.inputf, 'r') as fin:
            header = fin.readline()
            col_csv = header.split(',')
            col_tsv = header.split('\t')
            if len(col_csv) > 1:
                return True
            elif len(col_tsv) > 1:
                return False
            else:
                sys.exit("Unable to read character matrix!\n")

    def find_missing_symbol(self):
        with open(self.inputf, 'r') as fin:
            header = fin.readline()
            for line in fin:
                end = len(line)
                if line[-2:-1] == '-':
                    return True
                else:
                    ind = line.rfind('-')
                    miss = line[ind:ind+2]
                    if self.is_csv:
                        if miss == '-,':
                            return True
                    else:
                        if miss == '-\t':
                            return True
        return False

    def find_column_names(self):
        with open(self.inputf, 'r') as fin:
            header = fin.readline().splitlines()[0]
            if self.is_csv:
                return header.split(',')
            else:
                return header.split('\t')

    def read_to_pandas_df(self):
        if self.df is not None:
            return

        if self.is_csv:
            self.df = pandas.read_csv(self.inputf, sep=',')
        else:
            self.df = pandas.read_csv(self.inputf, sep='\t')

        if self.is_miss_dash:
            self.df = self.df.replace('-', '-1')
    
        for ind, col in enumerate(self.columns):
            if ind > 0:
                self.df[col] = self.df[col].astype(int)

    def read_to_dict(self):
        if self.cell2data is not None:
            return

        cell2data = {}
        with open(self.inputf, 'r') as fin:
            header = fin.readline()
            for line in fin:
                if self.is_csv:
                    [cell, data] = line.split(',', 1)
                else:
                    [cell, data] = line.split('\t', 1)

                data = data.splitlines()[0]
                if self.is_miss_dash:
                    data = data.replace('-', '-1')

                cell2data[cell] = data

        self.cell2data = cell2data

    def find_cells_with_same_seq(self):
        if self.eqclass is not None:
            return

        if self.cell2data is None:
            self.read_to_dict()

        cells = list(self.cell2data.keys())
        ncell = len(cells)
        cell_set = set(cells)
        eqclass = {}
        for i in range(0, ncell):
            ci = cells[i]
            if ci in cell_set:
                eqclass[ci] = set([])
                cell_set.remove(ci)
                lookat = list(cell_set)
                for cj in lookat:
                    if self.cell2data[ci] == self.cell2data[cj]:
                        eqclass[ci].add(cj)
                        cell_set.remove(cj)

        self.eqclass = eqclass

    def sort_cells_by_missingness(self, cell_list):
        if self.cell2data is None:
            self.read_to_dict()

        cell_list = sorted(cell_list)  # First sort alphabetically

        ncell = len(cell_list)
        nmiss = [0] * ncell

        for ind, cell in enumerate(cell_list):
            if self.is_csv:
                data = self.cell2data[cell].split(',')
            else:
                data = self.cell2data[cell].split()
            for x in data:
                if x == '-1':
                    nmiss[ind] += 1

        inds = numpy.argsort(nmiss)
        sorted_cell_list = [cell_list[ind] for ind in inds]

        return sorted_cell_list

    def find_cells_with_same_seq_upto_miss_helper(self):
        if self.cell2data is None:
            self.read_to_dict()

        if self.eqclass is None:
            self.find_cells_with_same_seq()

        cells_with_same_seq = list(self.eqclass.keys())
        order = self.sort_cells_by_missingness(cells_with_same_seq)

        cells = list(order)
        ncell = len(cells)
        cell_set = set(cells)
        eqclass = {}
        # found = 0
        for i in range(0, ncell):
            ci = cells[i]
            if self.is_csv:
                di = self.cell2data[ci].split(',')
            else:
                di = self.cell2data[ci].split()
            if ci in cell_set:
                eqclass[ci] = set([])
                cell_set.remove(ci)
                lookat = list(cell_set)
                for cj in lookat:
                    if self.is_csv:
                        dj = self.cell2data[cj].split(',')
                    else:
                        dj = self.cell2data[cj].split()
                    issame = True
                    for xi, xj in zip(di, dj):
                        if (xi != xj):
                            if xj == '-1':
                                pass
                            else:
                                issame = False
                    if issame:
                        eqclass[ci].add(cj)
                        cell_set.remove(cj)
        
        return eqclass

    def find_cells_with_same_seq_upto_miss(self):
        if self.eqclass_upto_miss is not None:
            return

        if self.eqclass is None:
            self.find_cells_with_same_seq()

        cells_same = self.eqclass
        cells_same_upto_miss = self.find_cells_with_same_seq_upto_miss_helper()

        eqclass = {}
        keep_list = list(cells_same_upto_miss.keys())
        for celli in keep_list:
            prune_set = cells_same_upto_miss[celli]
            prune_list = list(prune_set)
            for cellj in prune_list:
                prune_set = prune_set.union(cells_same[cellj])
            prune_set = prune_set.union(cells_same[celli])
            eqclass[celli] = prune_set

        self.eqclass_upto_miss = eqclass

    def create_pruned_df(self, isuptomiss=False):
        if isuptomiss:
            if self.df_pruned_upto_miss is not None:
                return
        else:
            if self.df_pruned is not None:
                return

        if self.cell2data is None:
            self.read_to_dict()

        if isuptomiss:
            if self.eqclass_upto_miss is None:
                self.find_cells_with_same_seq_upto_miss()
            eqclass = self.eqclass_upto_miss
        else:
            if self.eqclass is None:
                self.find_cells_with_same_seq()
            eqclass = self.eqclass

        keep_set = set(eqclass.keys())

        rows = []
        for cell in list(keep_set):
            if self.is_csv:
                row = [cell] + [int(x) for x in self.cell2data[cell].split(',')]
            else:
                row = [cell] + [int(x) for x in self.cell2data[cell].split('\t')]
            rows.append(row)

        df = pandas.DataFrame(rows, columns=self.columns)

        if isuptomiss:
            self.df_pruned_upto_miss = df
        else:
            self.df_pruned = df

    def write_pruned_matrix(self, keep_list, outputf):
        keep_set = set(keep_list)
        with open(outputf, 'w') as fout:
            fout.write(','.join(self.columns) + '\n')
            #if self.is_csv:
            #    fout.write(','.join(self.columns) + '\n')
            #else:
            #    fout.write('\t'.join(self.columns) + '\n')
            for cell in keep_list:
                data = self.cell2data[cell]
                if cell in keep_set:
                    if self.is_csv:
                        data = data.split(',')
                        #fout.write(cell + ',' + data + '\n')
                    else:
                        #fout.write(cell + '\t' + data + '\n')
                        data = data.split('\t')
                    fout.write(cell + ',' + ','.join(data) + '\n')

    def write_eqclass(self, outputf, isuptomiss=False):
        if isuptomiss:
            if self.eqclass_upto_miss is None:
                self.find_cells_with_same_seq_upto_miss()
            eqclass = self.eqclass_upto_miss
        else:
            if self.eqclass is None:
                self.find_cells_with_same_seq()
            eqclass = self.eqclass

        cell_list = list(eqclass.keys())

        with open(outputf, 'w') as fout:
            for key in cell_list:
                try:
                    vals = list(eqclass[key])
                    if len(vals) == 0:
                        fout.write(key + '\n')
                    else:
                        fout.write(key + ',' + ','.join(vals) + '\n')
                except KeyError:
                    pass

def infer_states_under_star_at_internal_helper(states_below):
    # Assume ambiguous state = -1
    #        unedited state = 0
    #        edited states = 1, 2, ...

    states_below_list = sorted(list(states_below))

    nbelow = len(states_below_list)
    if nbelow == 1:
        state = states_below_list[0]
    elif nbelow == 2:
        if -1 in states_below:
            state = states_below_list[1]
        else:
            state = 0
    else:
        state = 0

    return state


def infer_states_under_star_at_leaf(node, cmat):
    tmp = cmat.df.loc[cmat.df[cmat.columns[0]] == node.label].values[0][1:]

    states_at_leaf = []
    for s in tmp:
        if s == '-':
            states_at_leaf.append(-1)
        else:
            states_at_leaf.append(s)

    node.states = []
    node.states_below = []
    for state in states_at_leaf:
        state_set = set([state])
        node.states_below.append(state_set)
        node.states += list(state_set)


def infer_states_under_star_at_internal(node, children):
    nchar = len(children[0].states_below)
    node.states = []
    node.states_below = []

    for i in range(nchar):
        # Get states below
        states_below = children[0].states_below[i]
        for child in children[1:]:
            states_below = states_below.union(child.states_below[i])
        node.states_below.append(states_below)

        # Infer ancestral states
        state = infer_states_under_star_at_internal_helper(states_below)
        node.states += [state]

    for child in children:
        del child.states_below


def get_mutations_on_edge(head, tail):
    nchar = len(head.states)
    muts = []
    for i in range(nchar):
        hs = head.states[i]
        ts = tail.states[i]
        if ts == -1:
            pass
        elif hs == ts:
            pass
        else:
            if hs != 0:
                sys.exit("Error in star labeling!")
            muts.append((str("r%d" % i), ts))
    return muts


def get_mutations_above_root(root):
    nchar = len(root.states)
    muts = [] 
    for i in range(nchar):
        rs = root.states[i]
        if rs == -1:
            pass
        elif rs == 0:
            pass
        else:
            muts.append((str("r%d" % i), rs))
    return muts


def infer_muts_under_star_model(tree, cmat):
    for node in tree.traverse_postorder():
        if node.is_leaf():
            infer_states_under_star_at_leaf(node, cmat)
        else:
            # Check if node is binary
            children = node.child_nodes()

            infer_states_under_star_at_internal(node, children)

            # Annotate child edges with mutations
            for child in children:
                child.muts = get_mutations_on_edge(node, child)
                child.set_edge_length(len(child.muts))

    # Lastly process the root!
    root = tree.root
    root.muts = get_mutations_above_root(root)
    root.set_edge_length(len(root.muts))

    return tree


def contract_zero_length_branches(tree):
    # Set incoming edges incident to root and leaves
    # to be greater than 0 so they are not contracted
    rlen = tree.root.get_edge_length()
    if rlen == 0:
        tree.root.set_edge_length(0.5)

    for node in tree.traverse_leaves():
        if node.get_edge_length() == 0:
            node.set_edge_length(0.5)

    # Contract zero length edges
    for node in tree.traverse_postorder():
        if node.get_edge_length() == 0:
            node.contract()

    # Re-set edge lengths associated with root and leaves
    for node in tree.traverse_leaves():
        if node.get_edge_length() == 0.5:
            node.set_edge_length(0.0)

    if rlen == 0:
        tree.root.set_edge_length(0.0)


def score_tree_under_star_model(tree, priors=None):
    mutcounts = {}
    for node in tree.traverse_postorder():
        muts = node.muts
        for mut in muts:
            char, state = mut
            try:
                mutcounts[char][state] = 0
            except KeyError:
                mutcounts[char] = {}
                mutcounts[char][state] = 0

    score = 0.0
    wscore = 0.0
    for node in tree.traverse_postorder():
        muts = node.muts
        for mut in muts:
            char, state = mut
            mutcounts[char][state] += 1
            score += 1.0
            if priors is not None:
                xdf = priors[(priors["character"] == char) & \
                             (priors["state"] == state)]
                if xdf.shape[0] == 0:
                    prob = 0.0
                elif xdf.shape[0] == 1:
                    prob = xdf.probability.values[0]
                else:
                    sys.exit("Bad mutation priors!\n")
                wscore += -1.0 * numpy.log(prob)

    nmuts = 0            # Number of chars/states
    nmuts_homoplasy = 0  # Number of chars/states with homoplasy
    chars = list(mutcounts.keys())
    for char in chars:
        states = list(mutcounts[char].keys())
        for state in states:
            nmuts += 1
            if mutcounts[char][state] > 1:
                nmuts_homoplasy += 1

    if priors is None:
        scores = [nmuts, nmuts_homoplasy, score]
    else:
        scores = [nmuts, nmuts_homoplasy, score, wscore]

    return scores, mutcounts


def write_mutation_counts(mutcounts, cmat, priors, outputf):
    # Get pruned version of the data frame
    cmat.create_pruned_df()
    cmat.create_pruned_df(isuptomiss=True)

    with open(outputf + "_mutation_counts.csv", 'w') as fout:
        fout.write("character,state,nmuts_on_tree,ncell,ncell_eqclass,ncell_eqclass_upto_miss")
        if priors is not None:
            fout.write(",prior")
        fout.write('\n')

        chars = list(mutcounts.keys())
        for char in chars:
            states = list(mutcounts[char].keys())
            for state in states:
                fout.write("%s,%d" % (char, state))
                fout.write(",%d" % mutcounts[char][state])

                char_data = cmat.df[char].values
                ncell = len(numpy.where(char_data == state)[0])
                fout.write(",%d" % ncell)

                char_data = cmat.df_pruned[char].values
                ncell = len(numpy.where(char_data == state)[0])
                fout.write(",%d" % ncell)

                char_data = cmat.df_pruned_upto_miss[char].values
                ncell = len(numpy.where(char_data == state)[0])
                fout.write(",%d" % ncell)

                if priors is not None:
                    xdf = priors[(priors["character"] == char) & \
                                 (priors["state"] == state)]
                    if xdf.shape[0] == 0:
                        prob = 0.0
                        fout.write(",0.0")
                    elif xdf.shape[0] == 1:
                        prob = xdf.probability.values[0]
                        fout.write(",%f" % prob)
                    else:
                        sys.exit("Bad mutation priors!\n")

                fout.write('\n')


def main(args):
    "Works except when file starts with comma"

    # Read files
    cmat = CharacterMatrix(args.chars)
    cmat.read_to_pandas_df()

    tree = treeswift.read_tree(args.tree, schema="newick")
    if args.priors is None:
        priors = None
    else:
        priors = pandas.read_csv(args.priors, sep=',')

    ## Infer mutations on branches under the star model
    infer_muts_under_star_model(tree, cmat)

    ## Compute the star score
    scores, mutcounts = score_tree_under_star_model(tree, priors)
    if args.output:
        # Write information about mutations
        write_mutation_counts(mutcounts, cmat, priors, args.output)

        # Write tree after contracting branches with no mutations
        contract_zero_length_branches(tree)
        tree.write_tree_newick(args.output + "_contracted_tree.nwk")

    # Write score to standard out
    if priors is None:
        [nmuts, nmuts_homoplasy, uscore] = scores
        sys.stdout.write("%d,%d,%d\n" % (nmuts, nmuts_homoplasy, uscore))
    else:
        [nmuts, nmuts_homoplasy, uscore, wscore] = scores
        sys.stdout.write("%d,%d,%d,%f\n" % (nmuts, nmuts_homoplasy, uscore, wscore))

    sys.stdout.flush()
    os._exit(0)  # CRITICAL ON BLUE WATERS LOGIN NODE


if __name__ == "__main__":
    sys.exit("Not well tested... don't use")

    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--tree", type=str, default=None, required=True,
                        help="Input tree in newick format")

    parser.add_argument("-m", "--chars", type=str, default=None, required=True,
                        help="Input character matrix in TSV or CSV format")

    parser.add_argument("-w", "--priors", type=str, default=None, required=False,
                        help="Mutation priors in CSV format")

    parser.add_argument("-o", "--output", type=str, default=None, required=False,
                        help="Output name")

    main(parser.parse_args())
