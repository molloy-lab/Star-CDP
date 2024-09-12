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
#ifndef READER
#define READER
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<boost/unordered_map.hpp>
#include<cstdint>
#include<cstring>
#include<limits>
#endif

const long double INF = std::numeric_limits<long double>::infinity();

// https://www.geeksforgeeks.org/cpp-string-to-vector-using-delimiter/
std::vector<std::string> split(std::string str, std::string delimiter) {
    std::vector<std::string> v;
    if (!str.empty()) {
        int start = 0;
        do {
            // Find the index of occurrence
            int idx = str.find(delimiter, start);
            if (idx == std::string::npos) {
                break;
            }

            // If found add the substring till that
            // occurrence in the vector
            int length = idx - start;
            v.push_back(str.substr(start, length));
            start += (length + delimiter.size());
        } while (true);
        v.push_back(str.substr(start));
    }

    return v;
}


void read_characters_matrix(const std::string& filename, unsigned int &n, unsigned int &m, 
boost::unordered_map<std::string, unsigned int> &label2index,
std::vector<std::string> &labels, std::vector<std::vector<int>>
&index2_leaf_labeling, std::vector<std::vector<int>> &charbytaxa, 
std::unordered_set<std::string>outgroup_set) {
    
	std::ifstream file;
	
    file.open(filename);
    

	if (file.fail()) {
        std::cout << "Failed to open the character matrix file" << std::endl;
        exit(EXIT_FAILURE);

    }
    
	unsigned int cur_num = 0;
    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> row = split(line, ",");
        if (cur_num == 0) {
            m = row.size() - 1;

        } else {
            std::string taxa = row[0];
            label2index[taxa] = cur_num - 1;
            labels.push_back(taxa);
            
			if (outgroup_set.find(taxa) != outgroup_set.end()) {
				outgroup_set.erase(taxa);
			}

			row.erase(row.begin());
			
			std::vector<int> states;
			states.clear();
			for (auto &x : row) {
				states.push_back(std::stoi(x));
			}
            index2_leaf_labeling.push_back(states);
        }
        cur_num++;
    }
   
   std::vector<int> unedited_states(m, 0);
   
   for (auto& outg : outgroup_set) {
	label2index[outg] = cur_num - 1;
	cur_num++;
	labels.push_back(outg);
	index2_leaf_labeling.push_back(unedited_states);
   }
   
   n = index2_leaf_labeling.size();
   charbytaxa.resize(m);
   
   for (int i = 0; i < m; i++) {
   
       charbytaxa[i].resize(n);
        for (int j = 0; j < n; j++) {
            charbytaxa[i][j] = index2_leaf_labeling[j][i];
            
        }
   }
   
}

unsigned int read_mutation_prob(const std::string& filename,
std::vector<std::unordered_map<int,long double>> &mut_char_by_state, unsigned int m) {

    std::ifstream file(filename);
    
    unsigned int r = 0;
    
    if (file.fail()) {
        std::cout << "Failed to open the mutation probability file" <<
        std::endl;
        exit(EXIT_FAILURE);
    }

    unsigned int cur_num = 0;
    //mut_char_by_state.resize(m);
	
	std::unordered_set<std::string> states_set;
	std::string line;
	
	while(std::getline(file,line)) {
		if (cur_num > 0) {
				std::vector<std::string> row = split(line, ",");
				std::string char_name = row[0];
				char_name.erase(0,1);
				int char_id = std::stoi(char_name);
				states_set.insert(row[1]);

				mut_char_by_state[char_id][std::stoi(row[1])] = std::stold(row[2]);
		}
		cur_num++;
	}

	return states_set.size();

}

std::unordered_set<Bipartition> get_binary_clades(unsigned int n, unsigned int m, unsigned int r, std::vector<std::vector<int>> &charbytaxa, std::vector<std::string> &labels, boost::unordered_map<std::string, unsigned int> &label2index, std::vector<std::vector<Bipartition>> &char_trees) {
		
	std::unordered_set<Bipartition> clades;
	for (unsigned int i = 0; i < m; i++) {
			std::unordered_map<int,std::vector<std::string>> cluster;
			
			for (unsigned int j = 0; j < n; j++) {
					if (charbytaxa[i][j] != -1) {
							if (cluster.find(charbytaxa[i][j]) == cluster.end()) {
									cluster[charbytaxa[i][j]] = std::vector<std::string>();
							}
							cluster[charbytaxa[i][j]].push_back(labels[j]);
					}
			}
			// idnex 0 : state 0 index 1 : state -1
			std::vector<Bipartition> cur_char_tre;
			cur_char_tre.clear();
			for (auto& pair : cluster) {
				std::vector<std::string> group = pair.second;
				if (group.size() > 1) {
						boost::dynamic_bitset<> bs(label2index.size());
					for (std::string lab : group) {
						bs.set(label2index[lab]);
					}
					Bipartition cla(bs);
					if (pair.first != 0) 
							clades.insert(cla);
					if (pair.first != 0) 
							cur_char_tre.push_back(cla);
				}
			}
			
			// adding backgound state 0 taxa
			std::vector<std::string> background = cluster[0];
			
			for (auto &lab : background) {
					boost::dynamic_bitset<> bs(label2index.size());
					bs.set(label2index[lab]);
					Bipartition cla(bs);
					cur_char_tre.push_back(bs);
			}

			char_trees.push_back(cur_char_tre);
			
	}
	return clades;
}


void read_anatomical_labes(
	std::string &primary_tumor,
	std::unordered_set<std::string> &outgroup_set,
	std::vector<std::string> &labels,
	boost::unordered_map<std::string, unsigned int> label2index,
	std::unordered_map<Bipartition, std::string> &taxon2anatomical, 
	std::unordered_set<std::string> &anatomical_labels,
	std::ifstream &labels_file) {
	std::string line;
    // Read the file line by line
    while (std::getline(labels_file, line)) {
        std::istringstream iss(line);
        std::string cellId, anatomicalLabel;

        // Split the line into cellId and anatomicalLabel
        if (iss >> cellId >> anatomicalLabel) {
			boost::dynamic_bitset<> bs(labels.size());
			bs.set(label2index[cellId]);
            taxon2anatomical[Bipartition(bs)] = anatomicalLabel;
			anatomical_labels.insert(anatomicalLabel);
        }
    }

	for (auto &outg : outgroup_set) {
		boost::dynamic_bitset<> bs(labels.size());
		bs.set(label2index[outg]);
		taxon2anatomical[Bipartition(bs)] = primary_tumor;
		anatomical_labels.insert(primary_tumor);
	}
}