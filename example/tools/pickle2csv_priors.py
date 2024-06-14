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

