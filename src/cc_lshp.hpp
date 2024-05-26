#ifndef CC_LSHP_HPP
#define CC_LSHP_HPP
#include"read_search_space.hpp"
#include"small_star_homoplasy.hpp"
#include<string>
#include<vector>
#include<algorithm>
#include<unordered_set>
#include<cstdint>
//#include<limits>
#include <chrono>
#include <iomanip>
#include <unordered_map>
#include <utility>
#include<queue>
#include<tuple>
#include <boost/multiprecision/cpp_int.hpp>
#include<ctime>
#include <random>

typedef Bipartition Clade;
typedef std::pair<Clade, Clade> CladePair;
typedef std::unordered_set<Bipartition> Clades_Set;
typedef std::unordered_map<Clade, std::vector<CladePair>> SIESTA;
typedef std::unordered_map<Clade, long double> DP_Table;

#endif
//const long double INF = std::numeric_limits<long double>::infinity();

using namespace boost::multiprecision;


// may or may not need in the future
struct PairHash {
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T, U>& pair) const {
        // Combine the hash values of the pair's elements
        return std::hash<T>()(pair.first) ^ std::hash<U>()(pair.second);
    }
};

struct PairEqual {
    template <typename T, typename U>
    bool operator()(const std::pair<T, U>& p1, const std::pair<T, U>& p2) const {
        // Compare the pair's elements for equality
        return p1.first == p2.first && p1.second == p2.second;
    }
};


std::vector<Clade> add2tree(std::vector<Clade> &want2add, std::vector<Clade> &tree_clades) {
		std::vector<Clade> res(tree_clades.begin(), tree_clades.end());
		for (auto &cla1 : want2add) {
			boost::dynamic_bitset<> cla1_bs = cla1.get_bitset();
			bool flag = true;
			for (auto &cla2 : res) {
				boost::dynamic_bitset<> cla2_bs = cla2.get_bitset();
				if (!cla2_bs.is_subset_of(cla1_bs) && !cla1_bs.is_subset_of(cla2_bs) && cla2_bs.intersects(cla1_bs)) {
						flag = false;
						break;
				}
			}
			if (flag) res.push_back(cla1);
		}
		return res;
}

Tree build_tree_from_compat_clades(std::vector<Clade> &original_clades,
std::vector<std::string> &labels, std::unordered_map<Clade, cpp_rational>
&frequency, bool annotate) {
	
		std::unordered_set<Clade> remove_duplicated(original_clades.begin(), original_clades.end());
    	
		std::vector<Clade> clades(remove_duplicated.begin(), remove_duplicated.end());
		
		Node *root = new Node();
    
    if (annotate) root->fre = 1.0;
    
    //std::cout << "start building a tree" << std::endl;
    //std::cout << "use the following clade" << std::endl;
     
    //for (Clade &cla : clades) {
    //    std::cout << cla.to_labels(labels) << std::endl;
    //}
    
    std::unordered_map<Clade, Node*> cal2node;
    
    for (auto &e : clades) {

        if (e.count() == 1) {
            std::string lab = e.to_labels(labels);
            Node *chil = new Node(lab);
            //std::cout << "leaf lab: "<< lab << std::endl;
            cal2node[e] = chil;
            root->add_child(chil);
            if (annotate) chil->fre = 1.0;
        }
    }

    
    sort(clades.begin(), clades.end(), [](Clade &s1, Clade &s2) {
        return s1.count() > s2.count();
    });
	
	//std::cout << "builded a star tree! " << std::endl;

    for (auto &e : clades) {
        
        if (e.count() == labels.size()) continue;
        if (e.count() == 1) break;


        std::vector<Clade> taxons;
        taxons.clear();

        for (int i = 0; i < labels.size(); i++) {
            if (e.contain_index(i)) {
                boost::dynamic_bitset<> bs(labels.size());
                bs.set(i);
                Clade tax(bs);
                taxons.push_back(tax);
            }
        }

        Node* new_parent = new Node();
        if (annotate) new_parent->fre = frequency[e];

        auto it = cal2node.find(taxons[0]);

        if (it == cal2node.end()) std::cout << "did not store in cal2node "  <<
            std::endl;

        Node* old_parent = cal2node[taxons[0]]->get_parent();


        old_parent->add_child(new_parent);

        for (auto &tax : taxons) {
            Node *chil = cal2node[tax];

            if (chil->get_parent() != old_parent) {
                std::cout << "Not Compatiable!" << std::endl;
                exit(EXIT_FAILURE);
            }
            /* assume compatiable */

            old_parent->remove_child(chil);

           new_parent->add_child(chil);
        }
    }
    Tree R = Tree(root);
    return R;
    }

Tree build_tree_from_compat_clades(std::vector<Clade> &clades, std::vector<std::string> &labels){
		std::unordered_map<Clade, cpp_rational> dummy_fre;
		return build_tree_from_compat_clades(clades, labels, dummy_fre, false);
}

Tree build_full_leaf_tree_from_uncomplement_clades(std::vector<Clade> &clades, std::vector<std::string> &labels, boost::unordered_map<std::string, unsigned int> &label2index) {
	std::vector<Clade> completed_clades(clades.begin(), clades.end());
	for (std::string lab : labels) {
		boost::dynamic_bitset<> bs(labels.size());
		bs.set(label2index[lab]);
		Bipartition cla(bs);
		completed_clades.push_back(cla);
	}
	return build_tree_from_compat_clades(completed_clades, labels); 
}

// Tree build_partial_leaf_tree_from_uncomplement_clades(std::vector<Clade> &clades, std::vector<std::string> &labels, boost::unordered_map<std::string, unsigned int> &label2index) {
		
// 		std::vector<Clade> completed_clades(clades.begin(), clades.end());
// 		boost::dynamic_bitset<> leaf_bs(labels.size());
// 		boost::dynamic_bitset<> leaves_bs;
// 		for (auto &cla : clades) {
// 			boost::dynamic_bitset<> cla_bs = cla.get_bitset();
// 			leaf_bs |= cla_bs;
// 		}

// 		for (int i = 0; i < labels.size(); i++) {
// 			if (leaf_bs.test(i)) {
// 				boost::dynamic_bitset<> leaf_bs(labels.size());
// 				leaf_bs.set(i);
// 				Bipartition leaf_bp(leaf_bs);
// 				completed_clades.push_back(leaf_bp);
// 			}
// 		}
// 		return build_tree_from_compat_clades(completed_clades, labels);
// }


long double dp(Clades_Set &Sigma, std::vector<std::vector<int>> charbytaxa,
DP_Table &f, SIESTA &I, Clade &S, std::vector<std::unordered_map<int,long double>>
&mut_char_by_state, std::vector<std::vector<int>> &index2_leaf_labeling,
std::vector<std::string> &labels, boost::unordered_map<std::string, unsigned
int> &label2index, bool equal_weight) {
    auto start = std::chrono::high_resolution_clock::now();
    
    unsigned int m = charbytaxa.size();
    unsigned int n = charbytaxa[0].size();

    std::vector<Bipartition> Sigma_vec(Sigma.begin(), Sigma.end());
    
    sort(Sigma_vec.begin(), Sigma_vec.end(), [](Clade &s1, Clade &s2) {
        return s1.count() < s2.count();
        });

  std::cout << std::endl;

  std::cout << "The size of the search space(subproblems): " << Sigma_vec.size() << std::endl;

  unsigned int subp_num = 1;
  
  std::unordered_map<Clade, std::vector<int>> clade2state; 
  
  for (auto &u : Sigma_vec) {

    I[u] = std::vector<CladePair>();
    

    if (u.count() == 1) {
        f[u] = 0.0;
        
        //std::string u_str = u.to_string();
        
        std::string taxon = u.to_labels(labels);
        
        //std::cout << "base case: " << taxon << std::endl;
        
        clade2state[u] = index2_leaf_labeling[label2index[taxon]];
       // std::cout << "states size: " << clade2state[u].size() << std::endl;
        //std::cout << "state: " << clade2state[u][0] << std::endl;
        
        if(subp_num % 1000 == 0) std::cout << std::setw(12) << " subproblems computed." << std::endl;
        std::cout << "\rComputing subproblem:" << std::setw(12) << subp_num++;
        std::cout << std::flush;
    } else if (u.count() >= 2) {
        Clades_Set memo;
        memo.clear();
        clade2state[u] = get_state(u, n, m, charbytaxa);   
       
       long double tmp = INF;
        
        for (auto &a : Sigma_vec) {
            
            if (a.count() >= u.count()) {
                break;
            }

            if (a.get_bitset().is_subset_of(u.get_bitset())) {
                        
                Clade a_comp = a.complement(u);

                auto ptr_to_memo = memo.find(a_comp);

                if (ptr_to_memo != memo.end()) {
                    continue;
                } else {
                    memo.insert(a);
                
                }

                auto it1 = Sigma.find(a_comp);
                
                if (it1 != Sigma.end()) {
                    
                    long double v = f[a] + f[a_comp] + score(clade2state[u],
                    clade2state[a], clade2state[a_comp], mut_char_by_state, m, equal_weight);
                    
                    if (tmp > v) {
                        I[u].clear();
                        tmp = v;

                        I[u].push_back(std::make_pair(a, a_comp));
                        
                    } else if (tmp == v) {
                        
                       I[u].push_back(std::make_pair(a, a_comp));
                    }

                }

            }
        }

        f[u] = tmp;

        if (subp_num % 1000 == 0) std::cout << std::setw(12) << subp_num <<
            "subsubproblems computed." << std::endl;
            std::cout << "\rComputing subsubproblem:" << std::setw(12) << subp_num++;
            std::cout.flush();
    }
  }

    std::cout << std::endl;
    std::cout << "The star homoplasy score: " << f[S] << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "execution time for DP part: " << duration.count() << "ms" << std::endl;
     return f[S];
}

std::tuple<long double,SIESTA> cclshp(Clades_Set &Sigma,
std::vector<std::vector<int>> charbytaxa, std::vector<std::unordered_map<int,long double>>
&mut_char_by_state, std::vector<std::vector<int>> &index2_leaf_labeling,
std::vector<std::string> &labels, boost::unordered_map<std::string, unsigned
int> &label2index, bool equal_weight) {
    DP_Table f;
    SIESTA I;
    /* compute full set */
    
    boost::dynamic_bitset<> tbs(labels.size());

    tbs.flip();

    Bipartition S(tbs);

    long double star_score = dp(Sigma, charbytaxa, f, I, S, mut_char_by_state, index2_leaf_labeling, labels, label2index, equal_weight);
    
    return std::tuple<long double, SIESTA>{star_score, I};
}

Tree one_solution(SIESTA &I, std::vector<std::string> &labels,
std::unordered_map<Clade, cpp_rational> &freq) {
    
    boost::dynamic_bitset<> tbs(labels.size());
    
    std::vector<Clade> clades;
    
    tbs.flip();

    Clade S(tbs);
    
    clades.push_back(S);
    
    std::queue<Clade> que;

    que.push(S);
    
    while (!que.empty()) {
        Clade t = que.front();
        que.pop();
        clades.push_back(t);
        if (I[t].size() > 0) {
            Clade x = I[t][0].first;
            Clade y = I[t][0].second;
        
            que.push(x);
            que.push(y);
        }
    }
    
    std::cout << "Compute all clades in one sol "<< std::endl;
    return build_tree_from_compat_clades(clades, labels, freq, true);

}

std::unordered_map<Clade, cpp_rational> compute_below(SIESTA &I, Clades_Set
&Sigma,std::vector<std::string> &labels) {
    std::vector<Clade> Sigma_vec(Sigma.begin(), Sigma.end());
    
    std::unordered_map<Clade, cpp_rational> Below;
    std::sort(Sigma_vec.begin(), Sigma_vec.end(), [](Clade &s1, Clade &s2){
    return s1.count() < s2.count();
    });
    
    for(int i = 0; i < Sigma_vec.size(); i++) {
        Clade A = Sigma_vec[i];
        Below[A] = 0;
        if (Sigma_vec[i].count() == 1) {
            Below[A] = 1;
        } else {
            for (int j = 0; j < I[A].size(); j++) {
                Clade x = I[A][j].first;
                Clade y = I[A][j].second;
                
                Below[A] += Below[x]*Below[y];

            }

        }
    }
    return Below;

}


std::unordered_map<Clade, cpp_rational> compute_above(SIESTA &J, Clades_Set &Sigma,
std::vector<std::string> &labels, std::unordered_map<Clade, cpp_rational> &Below)
{   
    std::unordered_map<Clade, cpp_rational> Above;
    std::vector<Clade> Sigma_vec(Sigma.begin(), Sigma.end());
    sort(Sigma_vec.begin(), Sigma_vec.end(), [](Clade &s1, Clade &s2) {
            return s1.count() > s2.count();
        });
    for (int i = 0; i < Sigma_vec.size(); i++) {
        
        Clade cla = Sigma_vec[i];
        
        Above[cla] = 0;

        if (cla.count() == labels.size()) {
            Above[cla] = 1;        
        } else {
            for (int j = 0; j < J[cla].size(); j++) {
                Clade x = J[cla][j].first;
                Clade y = J[cla][j].second;
                if (!y.get_bitset().is_subset_of(x.get_bitset())) {
                    std::cout << "error in compute above" << std::endl;
                }
                Above[cla] += Above[x]*Below[y];
            }

        }

    }
    return Above;
}



cpp_rational compute_frequency(SIESTA &I, Clades_Set &Sigma, Clade &S,
std::unordered_map<Clade, cpp_rational> &fre, std::vector<std::string> &labels)
{   


    SIESTA J;
    
    std::unordered_map<Clade, bool> visited;
    for (auto &u : Sigma) {
        J[u] = std::vector<CladePair>();
        
        visited[u] = false;
    }
    
    std::unordered_map<Clade, cpp_rational> Below = compute_below(I, Sigma, labels);
    
    std::queue<Clade> que;
    que.push(S);
    
    visited[S] = true;
     


    while(!que.empty()) {
        Clade t = que.front();
        que.pop();
        
        for (int i = 0; i < I[t].size(); i++) {
            
            
            Clade x = I[t][i].first;
            Clade y = I[t][i].second;
                
                if (!y.get_bitset().is_subset_of(t.get_bitset())) {
                
                    std::cout << "error in compute J" << std::endl;
                    exit(0);   
                }
                J[x].push_back(std::make_pair(t, y));
                J[y].push_back(std::make_pair(t, x));    
                
                if (!visited[x]) {
                    que.push(x);
                    visited[x] = true;
                }
                if (!visited[y]) {
                    que.push(y);
                    visited[y] = true;
                }
            }
        
    }
    
    
    
    
    
    std::unordered_map<Clade, cpp_rational> Above = compute_above(J, Sigma, labels, Below);
       for (auto &e : Sigma) { 
            cpp_rational abv = (Above[e]);
            cpp_rational bels = (Below[S]);
            cpp_rational bele = (Below[e]);
            fre[e] = abv / bels * bele;
            if (fre[e] > 1) {
               // std::cout << "Error frequency > 1" << std::endl;
            }
        }
    // std::cout << "The number of optimal solutions: " << Below[S] << std::endl;
    return Below[S];
}



std::vector<Clade> alpha_clades(cpp_rational alpha, Clades_Set &Sigma,
std::unordered_map<Clade, cpp_rational> &fre) {
        
        std::vector<Clade> R;
        
        for (auto &e : Sigma) {
            if (fre[e] > alpha) {
                R.push_back(e);
            }
        }
        
        return R;
}


std::vector<Clade> get_compatiable(std::vector<Clade> &clades) {
    
    std::vector<Clade> compat;

    for (Clade &cla : clades) {
        boost::dynamic_bitset<> cla_bs = cla.get_bitset();
        
        bool good = true;

        for (Clade& other :  compat) {
            boost::dynamic_bitset<> other_bs = other.get_bitset();
            
            if (!other_bs.is_subset_of(cla_bs) && !cla_bs.is_subset_of(other_bs)
                && other_bs.intersects(cla_bs)) {
                good = false;
                break;
            }
            
        }
        
        if (good) compat.push_back(cla);

    }

    return compat;
}


std::vector<Clade> get_clades(Node *r, std::vector<std::string> &labels,
boost::unordered_map<std::string, unsigned int> &label2index) {
        std::vector<Clade> R;;
        for (auto j = Traverse::PostOrder(r); j != j.end(); j++) {
            if ((*j)->is_leaf()) {
                boost::dynamic_bitset<> bs(labels.size());
                (*j)->update_label_list((*j)->label);
                bs.set(label2index[(*j)->label]);
                Clade cla(bs);
                R.push_back(cla);
                //std::cout << "leaf node: " << (*j)->label << std::endl;
            } else {
                std::list<std::string> labels_list = (*j)->get_label_list();
                
                boost::dynamic_bitset<> bs(labels.size());
                
                for (std::string label: labels_list) {
                    bs.set(label2index[label]);
                }
                Clade cla(bs);
                R.push_back(cla);
                /*
                std::cout << "non leaf node: " << std::endl;
                for (const std::string & lab : labels_list) {
                    std::cout << lab << " ";
                }
                std::cout << std::endl;
                */
            }
            if (!(*j)->is_root()) {
                Node *par = (*j)->get_parent();
                par->update_label_list((*j)->get_label_list());
            }
        }
        return R;
    }

/*
std::vector<Clade> get_clades(Node *r, std::vector<std::string>
&labels,boost::unordered_map<std::string, unsigned int> &label2index) {
    std::vector<Clade> R;
    if (r->is_leaf()) {
        r->update_label_list(r->label);
    } else {
        std::list <Node*> children = r->get_children();
        for (auto child : children) {
            std::vector<Clade> t = get_clades(child, labels, label2index);
            R.insert(R.end(), t.begin(), t.end());
            r->update_label_list(child->label_list);
        }
    }
    boost::dynamic_bitset<> bs(labels.size());
    for (std::string lab : r->label_list) {
        bs.set(label2index[lab]);
    }
    Clade cla(bs);
    R.push_back(cla);
    return R;
}
*/
bool is_good_clade(Clade &cla, std::vector<Clade> &clades, Clade &Le) {
	bool good = true;
	boost::dynamic_bitset<> cla_bs = cla.get_bitset();
	boost::dynamic_bitset<> Le_bs = Le.get_bitset();
	
    if (!cla_bs.is_subset_of(Le_bs)) return false;
	
    for (Clade &x: clades) {
		boost::dynamic_bitset<> x_bs = x.get_bitset();
		if (!x_bs.is_subset_of(cla_bs) && !cla_bs.is_subset_of(x_bs) && x_bs.intersects(cla_bs)) {
			good = false;
			break;
		}

	}
	return good;
}

std::vector<Clade> get_good_clades(std::vector<Clade> clades_S,
std::vector<Clade> clades_T, Clade Le, std::vector<std::string> &labels) {
		std::unordered_set<Clade> R;
		for (Clade &cla : clades_S) {
			if (is_good_clade(cla, clades_T, Le)) {
				R.insert(cla);
			}
		}
		for (Clade &cla : clades_T) {
            R.insert(cla);
		}

        for (unsigned int i = 0; i < labels.size(); i++) {
            if (Le.contain_index(i)) {
                boost::dynamic_bitset<> bs(labels.size());
                bs.set(i);
                Clade cla(bs);
                R.insert(cla);
            }
        }
		std::vector<Clade>res(R.begin(), R.end());
		//std::cout << "# good clades: " << res.size() << std::endl;
        return res;
}
/*
void resolve_polytomies(Node *r) {
    if (r->num_children() > 2) {
        std::list<Node*> children = r->get_children();
        std::vector<Node*> vec(children.begin(), children.end());
        std::random_shuffle(vec.begin(), vec.end());
        r->children.assign(vec.begin(), vec.end());
        Node* left = r->children.back();
        r->remove_child(left);
        Node* right = new Node();
        right->add_child(left);
        right->add_child(r->children.back());
        r->remove_child(children.back());
        r->add_child(right);
        resolve_polytomies(right);

    }
    for (Node* child: r->get_children()) {
        resolve_polytomies(child);
    }
}
*/
Tree refine_by_guided(std::vector<Clade> &Clade_S, Tree *T, std::vector<std::string> &labels,
boost::unordered_map<std::string, unsigned int> &label2index) {
    std::vector<Clade> Clade_T = get_clades(T->get_root(), labels, label2index);
    Clade Le_T = Clade_T.back();
    
    std::vector<Clade> good_clades = get_good_clades(Clade_S, Clade_T, Le_T,
    labels);
    std::unordered_map<Clade, cpp_rational> fre;
    Tree R = build_tree_from_compat_clades(good_clades, labels, fre, false);
    //resolve_polytomies(R.get_root());
    return R;
}


Tree strict_consensus(std::unordered_map<Clade, cpp_rational> &fre, Clades_Set &Sigma,
std::vector<std::string> &labels) {
    
    std::vector<Clade> clades;

    for (auto &e : Sigma) {
        if (fre[e] >= 1) clades.push_back(e);
    }
    
    Tree R =  build_tree_from_compat_clades(clades, labels, fre, true);
    Node* root = R.get_root();
    return R;
}

Tree majority_consensus(std::unordered_map<Clade, cpp_rational> &fre, Clades_Set
&Sigma, std::vector<std::string> &labels) {
    
    std::vector<Clade> clades = alpha_clades(0.5, Sigma, fre);

    return build_tree_from_compat_clades(clades, labels, fre, true);

}


Tree greedy_consensus(std::unordered_map<Clade, cpp_rational> &freq,  Clades_Set
&Sigma, std::vector<std::string> &labels, boost::unordered_map<std::string,
unsigned int> &label2index) {
    
    std::vector<Clade> clades = alpha_clades(0.0, Sigma, freq);
    
    std::sort(clades.begin(), clades.end(), [&freq](const Clade &s1, const Clade  &s2){
        return freq[s1] > freq[s2];
});

    std::vector<Clade> compat = get_compatiable(clades);
    


    return build_tree_from_compat_clades(compat, labels, freq, true);

}





    
    

