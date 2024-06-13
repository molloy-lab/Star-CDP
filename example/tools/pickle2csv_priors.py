'''' 
Written by Junyan Dai (jdai1234@terpmail.umd.edu) in October 2023.
'''

import cassiopeia
import pickle
import argparse
import csv
import sys


def pickle2csv(input_file, output_file):
    with open(input_file, 'rb') as fp:
        dic = pickle.load(fp)

    with open(output_file, 'w') as fp:
        fp.write("character,state,probability\n")

        #states_set = set()
        #for character, inner_dict in dic.items():
        #    for state, prob in inner_dict.items():
        #        states_set.add(state)
        #states_set = sorted(states_set)

        for char, inner_dict in  dic.items():
            #for state in states_set:
            #    prob = 0.0
            #    if state in inner_dict.keys():
            #        prob = inner_dict[state]

            states_set = list(inner_dict.keys())
            states_set = sorted(states_set)

            for state in states_set:
                prob = inner_dict[state]
                fp.write("r%s,%s,%f\n" % (char, state, prob))


def main(args):
    if args.output:
        output = args.output
    else:
        output = args.input[0:-3] + 'csv'

    pickle2csv(args.input, output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str, help="Input pickle file", required=True)
    parser.add_argument("-o", "--output", type=str, help="Output CSV file", required=False)

    main(parser.parse_args())

