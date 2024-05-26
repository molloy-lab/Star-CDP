#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<boost/unordered_map.hpp>
#include<cstdint>
#include<cstring>
#include<limits>
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
&index2_leaf_labeling, std::vector<std::vector<int>> &charbytaxa) {
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
            
			row.erase(row.begin());
			
			std::vector<int> states;
			states.clear();
			for (auto &x : row) {
				states.push_back(std::stoi(x));
			}
            index2_leaf_labeling.push_back(states);

			//std::cout << " taxa: " << taxa << " ";
			//for (auto &x : row) {
			//	std::cout << x << " ";
			//}
			//std::cout << std::endl;
        }
        cur_num++;
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
	/*
    std::vector<std::vector<std::string>> table;
    std::string line;
    while (std::getline(file,line)) {
        std::vector<std::string> row = split(line, ",");
        
        if (cur_num > 0) {
            table.push_back(row);
            //std::cout << line << std::endl;
        }
        cur_num++;
    }
    std::string first_char = table[0][0];
       
    states.push_back("-1");
    states.push_back("0");
    state2index["0"] = 0;
    state2index["-1"] = 1;
    while (table[r][0] == first_char) { 
        states.push_back(table[r][1]);
        state2index[table[r][1]] = r + 2;
        r++;
    }

    cur_num = 0;
    
    for (int i = 0; i < m; i++) {
        mut_char_by_state[i].resize(r + 2);
        for (int j = 0; j < r + 2; j++) {
            if (j < 2) {
                mut_char_by_state[i][j] = 0;
            } else {

				if (table[cur_num][2] != "0") {
                	mut_char_by_state[i][j] = equal_weight ? 1 : -log(std::stold(table[cur_num][2]));
				} else {
					mut_char_by_state[i][j] = equal_weight ? 1 : INF;
				}
				cur_num++;
            }
        }

    }

    return r;*/
}
/*
std::vector<std::vector<Bipartition>> extract_clades_from_matrix(unsigned int n, unsigned int m, unsigned int r, boost::unordered_map<std::string, unsigned int> &state2index, std::vector<std::vector<std::string>> &charbytaxa, std::vector<std::string> &labels, boost::unordered_map<std::string, unsigned int> &label2index) {
		std::vector<std::vector<Bipartition>> collection;
		std::vector<Bipartition> clades;
	for (unsigned int i = 0; i < m; i++) {
		std::unordered_map<std::string,std::vector<std::string>> cluster;

		for (unsigned int j = 0; j < n; j++) {
			if (charbytaxa[i][j] != "-1") {
				cluster[charbytaxa[i][j]].push_back(labels[j]);
			}
		}
		// idnex 0 : state 0 index 1 : state -1

		for (unsigned int t  = 2; t < cluster.size(); t++) {
				if (cluster[t].size() > 1) {
						boost::dynamic_bitset<> bs(label2index.size());
					for (std::string lab : cluster[t]) {
						bs.set(label2index[lab]);
					}
					Bipartition cla(bs);
					clades.push_back(cla);
				}
		}
	}

	for (int i = 0; i < clades.size(); i++) {
		std::vector<Bipartition> x_y_collection;
		for (int j = i + 1; j < clades.size(); j++) {
			boost::dynamic_bitset<> x_bs = clades[i].get_bitset();
			boost::dynamic_bitset<> y_bs = clades[j].get_bitset();
			boost::dynamic_bitset<> intersect = x_bs & y_bs;
			if (!x_bs.is_subset_of(y_bs) && !y_bs.is_subset_of(x_bs) && intersect.count() > 0) {
				Bipartition inter(intersect);
				Bipartition x_minus_inter(x_bs - intersect);
				Bipartition y_minus_inter(y_bs - intersect);
				x_y_collection.push_back(intersect);
				x_y_collection.push_back(x_minus_inter);
				x_y_collection.push_back(y_minus_inter);
			}
		}
		if (x_y_collection.size() > 0) collection.push_back(x_y_collection);
	}
	return collection;

}
*/
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
			/*
			for (unsigned int t  = 2; t < cluster.size(); t++) {
					if (cluster[t].size() > 1) {
							boost::dynamic_bitset<> bs(label2index.size());
							for (std::string lab : cluster[t]) {
									bs.set(label2index[lab]);
							}
							Bipartition cla(bs);
							clades.insert(cla);
							cur_char_tre.push_back(cla);
					}
			}
			char_trees.push_back(cur_char_tre);*/
	}
	return clades;
}


// void write_newick_from(std::ostream &os, unsigned int n, unsigned int m,
// unsigned int r, std::vector<std::vector<std::string>> &charbytaxa, std::vector<std::string> &labels, boost::unordered_map<std::string, unsigned int> &label2index) {
		
// 		std::vector<Bipartition> clades;
		
// 		for (unsigned int i = 0; i < m; i++) {
        
//         	std::unordered_map<std::string,std::vector<std::string>> cluster;

//         	// cluster.resize(r + 2);
// 			for (unsigned int j = 0; j < n; j++) {
            
            
//             cluster[state2index[charbytaxa[i][j]]].push_back(labels[j]);
            
			
//             // if (charbytaxa[i][j] == "-1") {
//             //    std::cout << "visited -1" << std::endl;
//             //}
//             	if (charbytaxa[i][j] != "-1") {
//                 	cluster[charbytaxa[i][j]].push_back(labels[j]);
//                 	}
                
//         	}
//         	int count = 0;

//         	std::string str0 = "(";

//         	unsigned int first_non_zero_num = 0;

//         	while(cluster[first_non_zero_num].size() < 1) first_non_zero_num++;
//        		for (unsigned int k = 0; k < cluster[first_non_zero_num].size(); k++) {
//         	if (k == 0) {
//             	str0 += cluster[first_non_zero_num][k];

//         	} else {
//             	str0 += "," + cluster[first_non_zero_num][k];

//         	}
//        	}
    
//        	bool flag = false;
//        	for (unsigned int t = first_non_zero_num + 1; t < cluster.size(); t++) {
//             	std::string str1 = "(";

// 				if (cluster[t].size() <= 1) continue;
            	
// 				boost::dynamic_bitset<> bs(label2index.size());
				
// 				for (std::string lab : cluster[t]) {
// 					bs.set(label2index[lab]);
// 				}
				
// 				Bipartition cla(bs);

// 				clades.push_back(cla);		
				
// 				for (unsigned int k = 0; k < cluster[t].size(); k++) {
						
// 						if (k == 0) {
//                     		str1 += cluster[t][k];

//                 		} else {
//                     		str1 += "," + cluster[t][k];
//                 		}
//             	}
            
//             	if (str1 != "(") {
//                 	str0 += "," + str1 + ")";
//                 	if (!flag) flag = true;
//             	}

//        		}
    
//        		if (flag) os << str0 << ");" << std::endl;
       
//     	}
// 		/*	
// 		//one tree one clade
// 		std::cout << "The number of clades: " << clades.size() << std::endl;
// 		std::vector<Bipartition> clades_comp;
// 		for(int i = 0; i < clades.size(); i++) {
// 			for (int j = i + 1; j < clades.size(); j++) {
// 				boost::dynamic_bitset<> x_bs = clades[i].get_bitset();
// 				boost::dynamic_bitset<> y_bs = clades[j].get_bitset();
// 				boost::dynamic_bitset<> intersect = x_bs & y_bs;
				
// 				if (!x_bs.is_subset_of(y_bs) && !y_bs.is_subset_of(x_bs) && intersect.count() > 0) {
// 						std::cout << "find one conflict! " << std::endl;
// 					Bipartition inter(intersect);
// 					Bipartition x_minus_inter(x_bs - intersect);
// 					Bipartition y_minus_inter(y_bs - intersect);
// 					clades_comp.push_back(intersect);
// 					clades_comp.push_back(x_minus_inter);
// 					clades_comp.push_back(y_minus_inter);
// 				}
// 			}
// 		}
// 		std::cout << "Number of added extrac clades: " << clades_comp.size() << std::endl; 
// 		for(auto &clade : clades_comp) {
// 			std::string str0 = "(";
// 			int zero_count = 0;
// 			int one_count = 0;
// 			std::string str1 = "(";
			
// 			for (std::string lab : labels) {
// 				if (!clade.get_bitset().test(label2index[lab])) {
// 					if (zero_count == 0) {
// 						str0 += lab;
// 					} else {
// 						str0 += "," + lab;
// 					}
// 					zero_count += 1;
// 				} else {
// 					if (one_count == 0) {
// 						str1 += lab;
// 					} else {
// 						str1 += "," + lab;
// 					}
// 					one_count += 1;
// 				}
// 			}
// 			str0 += "," + str1 + ")";
// 			os << str0 << ");" << std::endl;
// 		}
		
// 		//one tree one clade 
// 		*/
// }

