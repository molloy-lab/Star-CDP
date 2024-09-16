#ifndef HERITABLE
#define HERITABLE
#include "migration.hpp"

#endif



long double compute_hertiable_cost(std::vector<int> &l, std::vector<int>
&l1, std::vector<int> &l2, std::vector<std::unordered_map<int,long double>>
&mut_char_by_state, int m, bool equal_weight)
{

    long double R = 0.0;
    for (int i = 0; i < m; i++) {
        if (l[i] != -1) {
            if (l1[i] == -1) R += 1;
            if (l2[i] == -1) R += 1;
        }
    }
    return R;
}



std::pair<long double, long double> dp(Clades_Set &Sigma, std::vector<std::vector<int>> charbytaxa,
DP_Table &f, SIESTA &I, Clade &S, std::vector<std::unordered_map<int,long double>>
&mut_char_by_state, std::vector<std::vector<int>> &index2_leaf_labeling,
std::vector<std::string> &labels, boost::unordered_map<std::string, unsigned
int> &label2index, bool equal_weight, 
std::unordered_map<Clade, std::vector<int>> &clade2state,
std::unordered_map<Clade, long double> &hertiable_cost) {
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
  
//   std::unordered_map<Clade, std::vector<int>> clade2state; 
  std::unordered_map<Clade, int>edges_num;

  for (auto &u : Sigma_vec) {

    I[u] = std::vector<CladePair>();
    

    if (u.count() == 1) {
        f[u] = 0.0;
        hertiable_cost[u] == 0;
        // edges_num[u] = 0;
        
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
    //    int best_edges_num = 0;
       long double tmp = INF;
       long double lowest_hertiable_cost = INF;
    // long double tmp = 0;
        
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
                    
                    long double w = hertiable_cost[a] + hertiable_cost[a_comp] + compute_hertiable_cost(clade2state[u],
                    clade2state[a], clade2state[a_comp], mut_char_by_state, m, equal_weight);
                    


                    if (tmp + lowest_hertiable_cost > v + w) {
                        I[u].clear();
                        tmp = v;
                        lowest_hertiable_cost = w;
                        // best_edges_num = cur_edges_num;
                        I[u].push_back(std::make_pair(a, a_comp));
                        
                    } else if (tmp + lowest_hertiable_cost == v + w) {
                        
                       I[u].push_back(std::make_pair(a, a_comp));
                    }

                }

            }
        }

        f[u] = tmp;
        hertiable_cost[u] = lowest_hertiable_cost;
        // edges_num[u] = best_edges_num;

        // std::cout << "Clade size: " << u.count() << " edges # " << best_edges_num << std::endl;
        if (subp_num % 1000 == 0) std::cout << std::setw(12) << subp_num <<
            "subsubproblems computed." << std::endl;
            std::cout << "\rComputing subsubproblem:" << std::setw(12) << subp_num++;
            std::cout.flush();
    }
  }

    long double total_parsimony_score = f[S];
    long double total_hertiable_score = hertiable_cost[S];
    std::cout << std::endl;
    std::cout << "The star homoplasy score: " << total_parsimony_score << std::endl;
    std::cout << "The hertiable parsimony score: " << total_hertiable_score << std::endl;
    std::cout << "The mixed parsimony score: " << total_parsimony_score + total_hertiable_score << std::endl; 
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "execution time for DP part: " << duration.count() << "ms" << std::endl;
     return std::make_pair(f[S],hertiable_cost[S]);
}








std::tuple<long double,long double,SIESTA> cclshp(Clades_Set &Sigma,
std::vector<std::vector<int>> charbytaxa, std::vector<std::unordered_map<int,long double>>
&mut_char_by_state, std::vector<std::vector<int>> &index2_leaf_labeling,
std::vector<std::string> &labels, boost::unordered_map<std::string, unsigned
int> &label2index, bool equal_weight, 
std::unordered_map<Clade, std::vector<int>> &clade2state,
std::unordered_map<Clade, long double> &hertiable_cost) {
    DP_Table f;
    SIESTA I;
    /* compute full set */
    
    boost::dynamic_bitset<> tbs(labels.size());

    tbs.flip();

    Bipartition S(tbs);

    std::pair<long double,long double> result_pair = dp(Sigma, 
    charbytaxa, f, I, S, mut_char_by_state, 
    index2_leaf_labeling, labels, label2index, 
    equal_weight, clade2state, hertiable_cost);
    
    long double total_star_score = result_pair.first;
    long double total_hertiable_cost = result_pair.second;

    return std::tuple<long double, long double,SIESTA>{total_star_score, total_hertiable_cost, I};
}


