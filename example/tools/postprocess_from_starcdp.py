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
