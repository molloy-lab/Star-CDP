/*

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

*/
#ifndef UTILITY
#define UTILITY
#include "bipartition.hpp"
#include<string>
#include<vector>
#include<algorithm>
#include<string.h>
#include<set>
#include<sstream>
#include<cstdint>
#include<boost/unordered_map.hpp>
#include<unordered_set>
#include<cmath>
#include<unordered_map>
#endif

long double weight_function(std::vector<std::unordered_map<int,long double>> &mut_char_by_state, int char_id, int state, bool equal_weight) {
	if (equal_weight) {
		return 1;
	}
	long double prob = mut_char_by_state[char_id][state];
	if (prob == 0) {std::cout << "Warning! want to read a 0 probability-> char_id: " << "c" << char_id << " state: " << state << std::endl;}
	
	return -log(prob);
}


// NOTE: did not store ancestral state 0 and missing state -1 in states set and
// mut_char_by_state.
long double score(std::vector<int> &l, std::vector<int>
&l1, std::vector<int> &l2, std::vector<std::unordered_map<int,long double>>
&mut_char_by_state, int m, bool equal_weight)
{

    long double R = 0.0;
    for (int i = 0; i < m; i++) {
        if (l[i] == 0) {
            if (l1[i] != 0 && l1[i] != -1) R += weight_function(mut_char_by_state, i, l1[i], equal_weight);
            if (l2[i] != 0 && l2[i] != -1) R += weight_function(mut_char_by_state, i, l2[i], equal_weight);
        }
    }
    return R;
}

std::vector<int> get_state(Bipartition A, unsigned int n, unsigned int m, std::vector<std::vector<int>>
&charbytaxa) {
    std::vector<int> R;
    for (int i = 0; i < m; i++) {
        std::unordered_set<int> memo;
        memo.clear();
        for (int j = 0; j < n; j++) {
            if (A.contain_index(j) && charbytaxa[i][j] != -1) memo.insert(charbytaxa[i][j]);
        }

        if (memo.size() > 1) {
            R.push_back(0);
        } else if (memo.size() == 1){
            R.push_back(*(memo.begin()));
        } else {
            R.push_back(-1);
        }

        /*
        bool all_missing = true;
        bool all_same = true;
        bool flag = false;
        std::string hold = "-1";
        for (int j = 0; j < n; j++) {
            if (A.contain_index(j) && charbytaxa[i][j] == "-1") {
                    continue;
                } else if (A.contain_index(j) && !flag) {
                    hold = charbytaxa[i][j];
                    all_missing = false;
                    flag = true;
                    } else if (A.contain_index(j) && charbytaxa[i][j] != hold) {
                        R.push_back("0");
                        //all_missing = false;
                        //all_same = false;
                        continue;
                }
        }
        if (all_missing) R.push_back("-1");
        if (all_same) R.push_back(hold);
    */
    }
    return R;
}


