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
    return tre.extract_tree_without(["FAKEROOT0", "FAKEROOT2"])

def main(args):
    nmap = read_name_map(args.name_map)
    with open(args.input, 'r') as inf, \
         open(args.output, 'w') as outf:
        for line in inf:
            tre = treeswift.read_tree_newick(line)
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
