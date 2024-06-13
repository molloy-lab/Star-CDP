import argparse
import treeswift
import sys

def get_leaf_set(tree):
    return set([leaf.label for leaf in tree.traverse_leaves()])


def get_branch_set(tree):
    branch_set = set([])

    leaves = sorted([leaf.label for leaf in tree.traverse_leaves()])
    nl = len(leaves)
    root = leaves[0]

    leaf_set = set(leaves)

    for node in tree.traverse_postorder():
        if node.is_leaf():
            pass
        elif node.is_root():
            pass
        else:
            leaves_below = sorted([leaf.label for leaf in node.traverse_leaves()])
            if leaves_below[0] == root:
                leaf_set_below = leaf_set.difference(set(leaves_below))
                leaves_below = sorted(list(leaf_set_below))

            nb = len(leaves_below)
            if (nb > 1) and (nb < nl - 1):
                branch = ','.join(leaves_below)
                branch_set.add(branch)

    return branch_set


def compare_two_trees(tree1, tree2):
    # Suppress unifurcations
    tree1.suppress_unifurcations()
    tree2.suppress_unifurcations()

    # Restrict trees to same leaf set
    leafset1 = get_leaf_set(tree1)
    leafset2 = get_leaf_set(tree2)

    keepset = leafset1.intersection(leafset2)

    tree1 = tree1.extract_tree_with(keepset)
    tree2 = tree2.extract_tree_with(keepset)

    # Get branch set for trees
    brset1 = get_branch_set(tree1)  # let t1 be true
    brset2 = get_branch_set(tree2)  # let t2 be esti

    # Compare tree topologies
    n1 = len(leafset1)
    n2 = len(leafset2)
    nl = len(keepset)
    i1 = len(brset1)
    i2 = len(brset2)
    fn = 0
    fp = 0

    for branch in brset1:
        if not branch in brset2:
            fn += 1

    for branch in brset2:
        if not branch in brset1:
            fp += 1

    tp = i1 - fn
    if tp != i2 - fp:
        sys.exit("Something bad happened!")

    return [n1, n2, nl, i1, i2, fn, fp, tp]


def main(args):
    # Get prefix of CSV line
    if args.prefix is None:
        prefix = ""
    else:
        prefix = str(args.prefix + ',')

    # Read trees
    with open(args.tree1, 'r') as fin:
        line = fin.readline()
        #while len(line.strip()) == 0:
        #    line = fin.readline()
        tree1 = treeswift.read_tree(line, "newick")

    with open(args.tree2, 'r') as fin:
        line = fin.readline()
        #while len(line.strip()) == 0:
        #    line = fin.readline()
        tree2 = treeswift.read_tree(line, "newick")

    # Compare two trees
    info = compare_two_trees(tree1, tree2)
    [n1, n2, nl, i1, i2, fn, fp, tp] = info

    # Write CSV to standard output
    sys.stdout.write("%s%d,%d,%d,%d,%d,%d,%d,%d\n" \
        % (prefix, n1, n2, nl, i1, i2, fn, fp, tp))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t1", "--tree1", type=str,
                        help="Input containing tree 1",
                        required=True)
    parser.add_argument("-t2", "--tree2", type=str,
                        help="Input containing tree 2",
                        required=True)
    parser.add_argument("-p", "--prefix", type=str,
                        help="Add prefix to CSV output",
                        required=False)

    main(parser.parse_args())
