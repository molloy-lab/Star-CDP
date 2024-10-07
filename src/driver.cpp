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

#include<cstdio>
#include<iostream>
#include<fstream>
#include<vector>
#include <cstdint>
#include<cstdlib>
#include<unistd.h>
#include<iostream>
#include<unordered_set>
#include "migration.hpp"

using namespace std;
const std::string help =
"===================================== Star-CDP =====================================\n"
"Star-CDP is a program that solves the large Star Homoplasy parsimony within a\n"
"clade-constrained version of tree space.\n\n"
"USAGE for large parsimony problem:\n"
"./star-cdp -i <input character file>\n"
"           -m <input mutation probability file>\n"
"           -t <input trees from heuritic search or other sources>\n"
"           -g <label of outgroup (unedited / ancestor)>\n"
"           -o <output tree file>\n\n"
"USAGE for small parsimony problem, i.e., computing score for given tree:\n"
"./star-cdp -i <input character file> \n"
"           -m <input mutation probability file>\n"
"           -q <input tree for scoring>\n\n"
"OPTIONS:\n"
"[-h|--help]\n"
"        Prints this help message.\n"
"(-i|--input) <input characters file>\n"
"        Name of file containing input characters in CSV format\n"
"(-m|--mutations) <mutations probability file>\n"
"        Name of file containing mutations probability\n"
"[(-t|--trees) <input trees file>]\n"
"        Name of file containing trees in newick format for constructing solution\n"
"        space with ASTRAL\n"
"(-x|-g|--outgroup) <outgroup or unedited cell label>\n"
"        Comma separated list of outgroup cells (e.g. unedited cell) used to root\n"
"        solution space\n"
"[(-q) <input species file>]\n"
"        Name of file containing tree for scoring in newick format\n"
//"[(-k)]\n"
//"        Write characters as rooted trees with one internal branch and exit\n"
"[(-o|--output) <output file>]\n"
"        Prefix of file for writing output tree (default: stdout)\n"
"[(-nosupp)]\n"
"        Turn off the calculation of clade support, which is the fraction of optimal\n"
"        solutions clade appears in\n"
"[(-XOUTG)]\n"
"        Output trees induced with only ingroup taxa\n"
"[(-consensus)]\n"
"        Compute greedy, majority, and strict consensus based on clade support\n"
"[(-e|--equal)]\n"
"        Using equal weight for all mutations\n"
"[(-memory)]\n"
"        Amount of memory given to ASTRAL.(Defualt:16000M)\n\n"
"Contact: Post issue to Github (https://github.com/molloy-lab/Star-CDP)\n"
"         or email Junyan Dai (jdai1234@umd.edu) & Erin Molloy (ekmolloy@umd.edu)\n\n"
"If you use Star-CDP in your work, please cite:\n"
"  Dai and Molloy, 2024, Star-CDP, https://github.com/molloy-lab/Star-CDP/\n"
"====================================================================================\n\n";

std::vector<std::string> get_outgroup(std::string oname) {
  std::vector<std::string> lat;
  boost::split(lat, oname, boost::is_any_of(","));
  return lat;
}


int main(int argc, char** argv) {
    std::cout << "Star-CDP version 1.0.0\nCOMMAND: ";

    for (int j = 0; j < argc; j++) std::cout << argv[j] << " ";
    std::cout << std::endl;
    
    if (argc == 1) {std::cout << help; return 0;}

    auto start = std::chrono::high_resolution_clock::now();
    
    string filename1 = ""; // input character matirx
    string filename2 = ""; // mutation probablity
    string filename3 = ""; // output file
    string filename4 = ""; // user defined search space
    string mig_file = ""; // anatomitical sites labels
    string input_tree = "";

    string outname;
	string guided_tree_filename = "";
    bool large_star = true;
    bool write_bptrees_and_exit = false;
    bool user_defined_search_space = false;
    bool nosupp = false;    
    bool heuristic = false;
    bool faster = false;
    bool equal_weight = false;
	bool greedy_mode = false;
	bool contract_mode = false;
    bool consensus = false;
    std::string leaf_map_file = "";
	string memory = "-Xmx16000M";
    bool no_outgroup = false;
    std::string primary_tumor = "";
    int num_sols = 0;
    for (int i = 0; i < argc; i++) {
    string opt(argv[i]);
    if (opt == "-h" || opt == "--help") {std::cout << help; return 0;}
    if (opt == "-i" || opt == "--input" && i < argc - 1) filename1 = argv[++ i];
    if (opt == "-q" && i < argc - 1) {large_star = false; input_tree = argv[++ i];}
    if (opt == "-k") {write_bptrees_and_exit = true;}
    if (opt == "-m" || opt == "--mutations" && i < argc - 1) {filename2 = argv[++ i];}
    if (opt == "-x" || opt == "--outgroup" || opt == "-g" && i < argc - 1) outname = argv[++ i];
    if (opt == "-o" || opt == "--output"&& i < argc - 1) filename3 = argv[++ i];
    if (opt == "-nosupp") nosupp = true;
    if (opt == "-consensus") consensus = true;
    if (opt == "-e" || opt == "--equal") equal_weight = true;
    if (opt == "-t" || opt == "--tress") {filename4 = argv[++ i];}
	if (opt == "-j") {guided_tree_filename = argv[++ i];}
	if (opt == "-contract" && i < argc - 1) {contract_mode = true; input_tree = argv[++ i]; large_star = false;}
	if (opt == "-memory" && i < argc - 1) {memory = "-Xmx" + std::string(argv[++ i]);}
    if (opt == "-XOUTG") no_outgroup = true;
    if (opt == "-leafmap" && i < argc - 1) {leaf_map_file = argv[++ i];}
    if (opt == "-MIG" && i < argc - 1) {mig_file = argv[++ i];}
    if (opt == "-p" && i < argc - 1) {primary_tumor = argv[++ i];}
    if (opt == "-trees" && i < argc - 1) {num_sols = std::stoi(argv[++ i]);}
	}
	
	if (argc < 3) {
			std::cout << "Not enough arguments given!\n\n";
			std::cout << help;
			return 1;
	}
	
	if (contract_mode) {
		std::cout << "Perform contract mode" << std::endl;
	}

  if (filename1 == "") {
    std::cout << "Need to provide file name for input characters!\n\n";
    std::cout << help;
    return 1;
  }

  if (filename1[0] == '-') {
     std::cout << "Warning: May not have correctly specified file name for -i option!\n\n";
  }

  if (outname == "" && large_star) {
    std::cout << "Need to provide outgroup name!\n\n";
    std::cout << help;
    return 1;
  }

   if (input_tree == "" && !large_star) {
    std::cout << "Need to provide file name for input species tree when using -q option!\n\n";
    std::cout << help;
    return 1;
  }

  if (input_tree[0] == '-') {
     std::cout << "Warning: May not have correctly specified file name for -q option!\n\n";
  }
    vector<string> outgroup = get_outgroup(outname);
	std::unordered_set<std::string> outgroup_set(outgroup.begin(), outgroup.end());
    unsigned int n;
    unsigned int m;
    unsigned int r;
    boost::unordered_map<std::string, unsigned int> label2index;

    std::vector<std::string> labels;
    std::unordered_set<std::string> ingroup;

    std::vector<std::vector<int>> index2_leaf_labeling;

    std::vector<std::vector<int>> charbytaxa;

    std::vector<std::string> states;

    std::unordered_map<Clade, std::string> taxon2anatomical;

    std::unordered_set<std::string>anatomical_labels;


    std::unordered_map<Clade, std::vector<int>> clade2state;

    std::unordered_map<std::string,std::vector<std::string>> leaves_eq_map;
    
    std::unordered_map<std::string, std::string> cell2anatomical;


    read_characters_matrix(filename1, n, m, label2index, labels,
    index2_leaf_labeling, charbytaxa, outgroup_set);
    
    std::cout << "Finish reading characters matrix" << std::endl; 
    std::vector<std::unordered_map<int,long double>> mut_char_by_state(m);
    
    for (const auto& taxa : labels) {
        if (outgroup_set.find(taxa) == outgroup_set.end()) {
            ingroup.insert(taxa);
        } 
    }

    if (mig_file != "") {
        std::ifstream anatomical_file_stream(mig_file);
        // Check if the file opened successfully
        if (!anatomical_file_stream.is_open()) {
            std::cerr << "Error opening sites labeling file!" << std::endl;
            return 1;
        } else {
            std::cout << "Reading anatomical sites" << std::endl;
            read_anatomical_labels(primary_tumor, outgroup_set,labels, label2index,taxon2anatomical, anatomical_labels, cell2anatomical, anatomical_file_stream);

        }

        anatomical_file_stream.close();

        if (leaf_map_file != "") {
            std::ifstream eqifs(leaf_map_file);

            if (!eqifs.is_open()) {
                std::cerr << "Error opeing leaves map file!" << std::endl;
            } else {
                std::cout << "Reading leaves map " << std::endl;
                load_eqclass(eqifs,leaves_eq_map);
                std::cout << "leaves_eq_map k#: " << leaves_eq_map.size() << std::endl;
                // for (const auto& [key, values] : leaves_eq_map) {cout << "Key: " << key << ", Values: ";for (const auto& val : values) {cout << val << " ";}cout << endl;}
                // std::cout << leaves_eq_map.size() << std::endl;
            }
            eqifs.close();
        }

        

    }


	//std::cout << "m: " << m << std::endl;
    //std::cout << "n: " << n << std::endl;
    
	if (!contract_mode) {
			r = read_mutation_prob(filename2, mut_char_by_state, m);
	}
	if (!contract_mode) std::cout << "Finish reading mutation probabilities" << std::endl;

    auto read_characters_end = std::chrono::high_resolution_clock::now();
	
    if (!large_star && !contract_mode) {
        auto R = star_homoplasy_score(charbytaxa, m, input_tree, mut_char_by_state, label2index, equal_weight);
        std::cout << "The star homoplasy  score: " << R.first << std::endl;
        if (nosupp) {
            std::cout << "The newick string with len: " << R.second->newick(false, true, false) << std::endl;
      } else {
            std::cout << "The newick string with len: " <<
            R.second->newick(false, true, true) << std::endl;
          }
    auto small_dollo_end = std::chrono::high_resolution_clock::now();
    auto duration_of_small_dollo =  std::chrono::duration_cast<std::chrono::milliseconds>(small_dollo_end - start);
    std::cout << "execution time of small star homoplasy: " << duration_of_small_dollo.count()  << "ms" << std::endl;
    return 0;

    }

	if (contract_mode) {
			Tree *tre = contract_mutless_edges(input_tree, charbytaxa, m, label2index);
			std::cout << "After contracting all mutationless edges: " << tre->newick(false, false, false) << std::endl;
			auto small_dollo_end = std::chrono::high_resolution_clock::now();
			auto duration_of_small_dollo =  std::chrono::duration_cast<std::chrono::milliseconds>(small_dollo_end - start);
			std::cout << "execution time of contraction of no mutation edges: " << duration_of_small_dollo.count()  << "ms" << std::endl;
			return 0;
	}

    if (filename4 == "") {
        filename4 = filename3 + "-auto-generated.bptrees";
		
        std::ofstream generated_clades_file(filename4);
    
    if (!generated_clades_file.is_open()) {
      std::cerr << "Error: could not create generated_caldes_file " << std::endl;
      return 1;
    }
        std::unordered_set<Clade> all_clades_set;
		std::vector<std::vector<Clade>> trees_clades;
		std::vector<std::vector<Clade>> char_trees_clades;
        std::unordered_set<Clade> binary_clades = get_binary_clades(n,m,r, charbytaxa, labels, label2index, trees_clades);
        std::vector<Clade> binary_clades_vec(binary_clades.begin(), binary_clades.end());
        // for (int i = 0; i < binary_clades_vec.size(); i++) {
        //     Clade c1 = binary_clades_vec[i];
        //     std::vector<Clade> tree1;
        //     tree1.push_back(c1);
        //     trees_clades.push_back(tree1);

        // }
        for (int i = 0; i < trees_clades.size(); i++) {
			Tree tre = build_full_leaf_tree_from_uncomplement_clades(trees_clades[i], labels, label2index);
			std::string nwk = tre.newick(false, false,false);
			generated_clades_file << nwk << std::endl;
		}

        if (write_bptrees_and_exit) return 0;

    }


    
	
    
    Clades_Set Sigma = read_search_space(filename4, label2index, labels, outgroup, memory);
    
    auto search_space_end = std::chrono::high_resolution_clock::now();
    std::tuple<long double, SIESTA> sol_pair =  cclshp(Sigma, charbytaxa, mut_char_by_state, index2_leaf_labeling, labels, label2index, equal_weight,
    clade2state);

    SIESTA I = std::get<1>(sol_pair);

    
    boost::dynamic_bitset<> tbs(labels.size());

    tbs.flip();

    Bipartition S(tbs);

     
    std::unordered_map<Clade, cpp_rational> fre;
    
    cpp_rational num_opt = compute_frequency(I, Sigma, S, fre, labels);
    std::cout << "The number of optimal solutions: " << num_opt << std::endl;
    Tree one_sol = one_solution(I,labels, fre); 

    if (no_outgroup) {
        one_sol = *one_sol.get_induced_subtree_copy(ingroup);
        cout << "NO OUT GROUP!!!" << std::endl;

    }

    auto duration_of_characters =  std::chrono::duration_cast<std::chrono::milliseconds>(read_characters_end - start);
  std::cout << "execution time of reading characters matrix: " << duration_of_characters.count()  << "ms" << std::endl;
  auto duration_of_search_space =  std::chrono::duration_cast<std::chrono::milliseconds>(search_space_end - read_characters_end);
  std::cout << "execution time of computing search space: " << duration_of_search_space.count() << "ms" << std::endl;
    

    bool write2file = filename3 != "";
    
    std::string one_sol_name = filename3 + "_one_sol.tre";
    std::string strict_name =  filename3 + "_strict_consensus.tre"; 
    std::string majority_name = filename3 + "_majority_consensus.tre";
    std::string greedy_name =  filename3 + "_greedy_consensus.tre";
    std::string mig_name = filename3 + "_migration.tre";
    std::string num_of_sol_name = filename3 + "_number_of_sol.csv";
    cout << "One optimal solution: " << endl;    
    
    if (write2file) {
        std::ofstream num_opt_out(num_of_sol_name);
        num_opt_out << num_opt << endl;
    }
    


    if (nosupp) {
        cout << one_sol.newick(false, false, false) << endl;
        if (write2file) {
            std::ofstream one_sol_out(one_sol_name);
            one_sol_out << one_sol.newick(false, false, false) << endl;
            }
    } else {
        cout << one_sol.newick(false, false, true) << endl;
        if (write2file) {
            std::ofstream one_sol_out(one_sol_name);
            one_sol_out << one_sol.newick(false, false, true) <<endl;
            
            }
    }
	if (consensus) {
            Tree strict = strict_consensus(fre, Sigma, labels);

        if (no_outgroup) {
            strict = *strict.get_induced_subtree_copy(ingroup);
        }
        cout << "The strict consensus tree: " << endl;
    
        if (nosupp) {
            cout << strict.newick(false, false, false) << endl;
            if (write2file) {
                std::ofstream strict_out(strict_name);
                strict_out << strict.newick(false, false, false) << endl;}

        } else {
            cout << strict.newick(false, false, true) << endl;
            if (write2file) {
            
                std::ofstream strict_out(strict_name);
                strict_out << strict.newick(false, false, true) << endl; }
        }

        Tree majority = majority_consensus(fre, Sigma, labels);
        
        if (no_outgroup) {
            majority = *majority.get_induced_subtree_copy(ingroup);
        }

        cout << "The majority consensus tree: " << endl;
        if (nosupp) {
            cout << majority.newick(false, false, false) << endl;
            if (write2file) {
                std::ofstream majority_out(majority_name);
                majority_out << majority.newick(false, false, false) << endl;
            }
        } else {
            cout << majority.newick(false, false, true) << endl;
            if (write2file) {
                std::ofstream majority_out(majority_name);
                majority_out << majority.newick(false, false, true) << endl;
                }
        }
        Tree greedy = greedy_consensus(fre, Sigma, labels, label2index);

        if (no_outgroup) {
            greedy = *greedy.get_induced_subtree_copy(ingroup);
        }

        cout << "The greedy consensus tree: " << endl;
        if (nosupp) {
            cout << greedy.newick(false, false, false) << endl;
            if (write2file) {
                std::ofstream greedy_out(greedy_name);
                greedy_out << greedy.newick(false, false, false) << endl;
            }
        } else {
            cout << greedy.newick(false, false, true) << endl;
            if (write2file) {
                std::ofstream greedy_out(greedy_name);
                greedy_out << greedy.newick(false, false, true) << endl;
            }
        }

    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "execution time: " << duration.count() << "ms" << std::endl;

    if (mig_file != "") {

        std::vector<std::string> anatomical_labels_vec(anatomical_labels.begin(), anatomical_labels.end());

        // computing leaves weight

        std::unordered_map<Bipartition, std::unordered_map<std::string, double>> taxon_weight_meta;
        std::unordered_map<Bipartition, std::unordered_map<std::string, double>> taxon_weight_reseeding;

        bool is_weight = leaf_map_file != "";

        // std::cout << " 475 pruned_cell2anatomical size: " << pruned_cell2anatomical.size() << std::endl;
        if (is_weight) {
            load_weights_for_leaves(taxon_weight_meta,
            taxon_weight_reseeding,
            primary_tumor,
            labels,
            anatomical_labels_vec,
            label2index,
            leaves_eq_map,
            cell2anatomical, 
            taxon2anatomical);
        } else {
            load_weights_for_leaves(taxon_weight_meta,taxon_weight_reseeding,labels,anatomical_labels_vec,label2index,taxon2anatomical);
        }
        

        // std::cout << "pruned_cell2anatomical size: " << pruned_cell2anatomical.size() << std::endl;
        // for (const auto& [k,v] : taxon_weight_meta) {
        //      std::cout << v.at(primary_tumor) << std::endl;

        // }


        std::pair<std::pair<double,double>, std::unordered_map<std::pair<Clade, std::string>,  std::pair<std::pair<Clade, Clade>, std::pair<std::string, std::string>>,PairHash, PairEqual>> mig_res = 
        mig_dp(primary_tumor, taxon2anatomical,anatomical_labels_vec,Sigma,S,clade2state,I, taxon_weight_meta, taxon_weight_reseeding);
        
        double final_meta = mig_res.first.first;
        double final_reseeding = mig_res.first.second;
        Tree mig_tree = binary_mig_solution(mig_res.second,fre, S,primary_tumor, labels, 
        taxon_weight_meta,
        taxon_weight_reseeding);

        if (no_outgroup) {
            mig_tree = *mig_tree.get_induced_subtree_copy(ingroup);

        }

        cout << " The most parsimonious refinement of migration tree: " << endl;

        if (nosupp) {
            cout << mig_tree.newick(false, false, false) << endl;
            if (write2file) {
                std::ofstream mig_out(mig_name);
                mig_out << mig_tree.newick(false, false, false) << endl;
            }
        } else {
            cout << mig_tree.newick(false, false, true) << endl;

            if (write2file) {
                std::ofstream mig_out(mig_name);
                mig_out << mig_tree.newick(false, false, true) << endl;
            }
        }

    }


    return 0;

}






