"""
This python script infers mutations under star homoplasy and then contracts mutationless branches

"""
import argparse
import pandas
import treeswift
import numpy
import os
import sys


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
    states_at_leaf = cmat.loc[node.label].values.tolist()
    
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
            muts.append((i, ts))
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
            muts.append((i, rs))
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


def main(args):
    # Read character matrix
    chars = pandas.read_csv(args.chars, index_col=0)

    # Read tree
    tree = treeswift.read_tree(args.tree, schema="newick")

    # Infer mutations under star homoplasy model
    infer_muts_under_star_model(tree, chars)

    # Contract branches with no mutations
    contract_zero_length_branches(tree)
    
    # Write tree with mutationless branches contracted
    if args.output is None:
        sys.stdout.write(tree.newick() + '\n')
    else:
        with open(args.output, 'w') as file:
            file.write(tree.newick() + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--tree", type=str, default=None,
                        help="File containing tree in newick format",
                        required=True)

    parser.add_argument("-m", "--chars", type=str, default=None,
                        help="Character data in CSV format",
                        required=True)

    parser.add_argument("-o", "--output", type=str, default=None,
                        help="Output file",
                        required=False)

    main(parser.parse_args())
