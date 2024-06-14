import argparse
import treeswift
from typing import Dict


def read_name_map(input:str, right2left=True) -> Dict[str,str]:
    nmap = {}
    with open(input, 'r') as fin:
        for line in fin:
            tmp = line.strip()
            [left, right] = tmp.split(',')
            if right2left:
                nmap[right] = left
            else:
                nmap[left] = right
    return nmap


def relabel_and_remove_outgroup(tre:treeswift.Tree, nmap:Dict[str, str]) -> treeswift.Tree:
    for leaf in tre.traverse_leaves():
        lab = leaf.label
        try:
            leaf.label = nmap[leaf.label]
        except KeyError:
            pass
    return tre.extract_tree_without(["FAKEROOT", "FAKEROOT2"])


def main(args):
    nmap = read_name_map(args.name_map)

    if args.output is None:
        outf = sys.stdout
    else:
        outf = open(args.output, 'w')

    with open(args.input, 'r') as inf:
        for line in inf:
            start = line.find("[&R]")
            if start > -1:
                tre = treeswift.read_tree_newick(line[start+4:])
                tre = relabel_and_remove_outgroup(tre, nmap)
                nwk = tre.newick().replace('[&R]', '') + '\n'
                outf.write(nwk)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str, default=None, required=True,
                        help="Input character matrix in TSV or CSV format")

    parser.add_argument("-n", "--name_map", type=str, default=None, required=True,
                        help="Leaf label map in CSV format")

    parser.add_argument("-o", "--output", type=str, default=None, required=True,
                        help="Output")

    main(parser.parse_args())
