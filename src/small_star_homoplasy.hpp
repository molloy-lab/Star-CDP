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
#ifndef SMALL_STAR_HOMOPLASY
#define SMALL_STAR_HOMOPLASY
#include "tree-lib.hpp"
#include "reader.hpp"
#include<queue>
#include<algorithm>
#include<string>
#include<cstdint>
#include<vector>
#include<fstream>
#include<list>
#include<array>
#include<unordered_set>
#include<utility>
#endif
// Assume in the input tree file there is only contain one newick string.
std::string get_newick(std::string input_tree) {
  std::ifstream file(input_tree);

  if (file.fail()) {
    std::cout << "Failed to open the input tree file" << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string line;

  std::getline(file, line);

  file.close();
  return line;
}



long double one_step(int i, std::vector<std::vector<int>>
&charbytaxa,boost::unordered_map<std::string, unsigned int> &label2index, Node*
r,std::vector<std::unordered_map<int, long double>>
&mut_char_by_state,  bool equal_weight) {
    long double res = 0.0;
    for (auto j = Traverse::PostOrder(r); j != j.end(); j++) {
        
        if ((*j)->is_leaf()) {
            
            std::string leaf_label = (*j)->label;

            (*j)->state.push_back(charbytaxa[i][label2index[leaf_label]]);

        } else {
            
			std::list<Node*> children = (*j)->get_children();
            
            std::unordered_set<int> states_set;
           	
			states_set.clear();

			for (auto child : children) {
					states_set.insert(child->state[i]);                
            }
			states_set.erase(-1);
			
			/*
			std::cout << " states_set: ";
			for (auto x : states_set) 
					std::cout << x << " ";
			std::cout << std::endl;
			*/

			if ((*j)->is_root()) {
				(*j)->state.push_back(0);
			} else if (states_set.size() > 1) {
                (*j)->state.push_back(0);
            } else if (states_set.size() == 1) {
                (*j)->state.push_back(*(states_set.begin()));
            } else if (states_set.size() == 0) {
                (*j)->state.push_back(-1);
            }

			//std::cout << " state we added in " << (*j)->state[i] << std::endl;
			
			if ((*j)->state[i] == 0) {
				for (auto &child: children) {
						//std::cout << "i: " << i << " child->state.size(): " << child->state.size() << std::endl;
					if (child->state[i] != -1 && child->state[i] != 0) {
							long double added = weight_function(mut_char_by_state, i, child->state[i], equal_weight);
							res += added;

					//if (child->is_leaf()) {
					//		std::cout << "This is a leaf child state: " << child->state[i] << " adding weight: "<< added << std::endl;
					//} else {
					//	std::cout << "This is not a leaf child state: " << child->state[i] << " adding weight: "<< added << std::endl;
					//}
					}
					if (child->state[i] != 0 && child->state[i] != -1) child->len += 1;
				}
			}

        }

    }
    return res;
}



std::pair<long double, Tree*> star_homoplasy_score(std::vector<std::vector<int>>
&charbytaxa, unsigned int m, std::string input_tree,
std::vector<std::unordered_map<int,long double>>&mut_char_by_state,
boost::unordered_map<std::string, unsigned int> &label2index, bool equal_weight) {
    
    long double tot = 0.0;

    std::string newick_str = get_newick(input_tree);

    Tree *T = new Tree(newick_str);
    for (int i = 0; i < m; i++) {
        tot += one_step(i, charbytaxa,
        label2index, T->get_root(), mut_char_by_state, equal_weight);
        }
    
    return std::make_pair(tot, T);

}

std::pair<long double, Tree*> star_homoplasy_score_from_tree(std::vector<std::vector<int>>
&charbytaxa, unsigned int m, Tree* T,
std::vector<std::unordered_map<int,long double>>&mut_char_by_state,
boost::unordered_map<std::string, unsigned int> &label2index, bool equal_weight) {

    long double tot = 0.0;
    std::vector<std::vector<std::string>> labelings;
	
    for (int i = 0; i < m; i++) {
        tot += one_step(i, charbytaxa,
        label2index, T->get_root(), mut_char_by_state, equal_weight);
        }
	
    return std::make_pair(tot, T);

}

Tree* contract_mutless_edges(std::string input_tree, std::vector<std::vector<int>> &charbytaxa, unsigned int m, boost::unordered_map<std::string, unsigned int> &label2index) {
		std::string nwk = get_newick(input_tree);
		Tree *T = new Tree(nwk);
		Node *r = T->get_root();
		for (int i = 0; i < m; i++) {	
			for (auto j = Traverse::PostOrder(r); j != j.end(); j++) {
						
						if ((*j)->is_leaf()) {
							std::string leaf_label = (*j)->label;
							(*j)->state.push_back(charbytaxa[i][label2index[leaf_label]]);
						} else {
							std::list<Node*> children = (*j)->get_children();
							std::unordered_set<int> states_set;
							states_set.clear();
							for (auto child : children) {
								if (child->state.back() == -1) {
										continue;
								} else {
									states_set.insert(child->state.back());
								}
							}
							if ((*j)->is_root()) {
								(*j)->state.push_back(0);
							} else if (states_set.size() > 1) {
								(*j)->state.push_back(0);
							} else if (states_set.size() == 1) {
								(*j)->state.push_back(*(states_set.begin()));
							} else {
								(*j)->state.push_back(-1);
							}
						}
			}
		}
	
		for (auto j = Traverse::PostOrder(r); j != j.end(); j++) {
				std::list<Node*> children = (*j)->get_children();
				std::vector<Node*> want2contract;
				want2contract.clear();
				for (auto child : children) {
					bool flag = true;
					for (int k = 0; k < m; k++) {
						if (child->state[k] != -1 && child->state[k] != (*j)->state[k]) {
							flag = false;
							break;
						}
					}
						if (flag && !child->is_leaf()) {want2contract.push_back(child);}
				}

				for (int i = 0; i < want2contract.size(); i++)
						want2contract[i]->contract();
		}
		return T;
}
