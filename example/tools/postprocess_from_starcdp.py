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

def main(args):
    suffixes = ["_one_sol.tre",
                "_greedy_consensus.tre",
                "_majority_consensus.tre",
                "_strict_consensus.tre"]
    for suffix in  suffixes:
        with open(args.input + suffix, 'r') as inf, \
             open(args.input + suffix + "-" + args.output, 'w') as outf:
            for line in inf:
                tre = treeswift.read_tree_newick(line)
                tre = tre.extract_tree_without(["FAKEROOT"])
                nwk = tre.newick().replace('[&R]', '') 
                outf.write(nwk.strip() + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str, default=None, required=True,
                        help="Input Star-CDP prefix")

    parser.add_argument("-o", "--output", type=str, default=None, required=True,
                        help="Output suffix")

    main(parser.parse_args())
