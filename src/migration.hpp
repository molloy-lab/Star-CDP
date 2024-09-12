#ifndef MIGRATION
#define MIGRATION
#include "cc_lshp.hpp"

#endif 

typedef std::tuple<std::string, std::string> pp;

int INT_INF = std::numeric_limits<int>::max();



// struct TupleHash {
//     template <typename... Args>
//     std::size_t operator()(const std::tuple<Args...>& tuple) const {
//         std::size_t seed = 0;
//         boost::hash_combine(seed, tuple);
//         return seed;
//     }
// };

// // Define a custom equality function for tuples
// struct TupleEqual {
//     template <typename... Args>
//     bool operator()(const std::tuple<Args...>& lhs, const std::tuple<Args...>& rhs) const {
//         return lhs == rhs;
//     }
// };









std::pair<int, std::unordered_map<std::pair<Clade, std::string>,  std::pair<std::pair<Clade, Clade>, std::pair<std::string, std::string>>,PairHash, PairEqual>> mig_dp(
    std::string &primary_tumor,
    std::unordered_map<Clade, std::string> &taxon2anatomical,
    std::unordered_set<std::string> &anatomical_labels,
    Clades_Set &Sigma,
    Clade &S,
    std::unordered_map<Clade,std::vector<int>> &clade2state,
    SIESTA &I
) {
    std::unordered_map<std::pair<Clade, std::string>, std::pair<int, int>,PairHash, PairEqual> mig_dp;

    std::unordered_map<std::pair<Clade, std::string>,  std::pair<std::pair<Clade, Clade>, std::pair<std::string, std::string>>,
                PairHash, PairEqual> trace_back;


    std::vector<Bipartition> Sigma_vec(Sigma.begin(), Sigma.end());
    
    sort(Sigma_vec.begin(), Sigma_vec.end(), [](Clade &s1, Clade &s2) {
        return s1.count() < s2.count();
        });


    for (auto &u : Sigma_vec) {
        for (auto &site_lab : anatomical_labels) {
            if (u.count() == 1) {
                if (site_lab == taxon2anatomical[u]) {
                    mig_dp[std::make_pair(u, site_lab)] = std::make_pair(0,0);
                } else {
                    mig_dp[std::make_pair(u, site_lab)] = std::make_pair(0x3f3f3f3f,0);
                }
            } else {
                
                int best_mig_score = 0x3f3f3f3f;
                int best_edges_num = 0;
                // std::vector<CladePair> candids;
                
                CladePair best_candids;

                std::pair<std::string, std::string> best_candids_labs;

                for (auto &children_tax : I[u]) {

                    for (auto &site_lab_l : anatomical_labels) {

                        for (auto &site_lab_r : anatomical_labels) {
                            int s1 = mig_dp[std::make_pair(children_tax.first,site_lab_l)].first;
                            int edges_num = mig_dp[std::make_pair(children_tax.first, site_lab_l)].second;
                            int s2 = mig_dp[std::make_pair(children_tax.second, site_lab_r)].first;
                            edges_num += mig_dp[std::make_pair(children_tax.second, site_lab_r)].second;

                            if (clade2state[children_tax.first] != clade2state[u]) {
                                edges_num += 1;
                                
                                }

                            if (clade2state[children_tax.second] != clade2state[u]) {
                                edges_num += 1;

                                }

                                int score = s1 + s2;

                                if (site_lab_l != site_lab) {
                                    score += 1;
                                }

                                if (site_lab_r != site_lab) {
                                    score += 1;
                                }

                        

                         if (score < best_mig_score) {
                            best_mig_score = score;
                            best_edges_num = edges_num;
                            best_candids = std::make_pair(children_tax.first,children_tax.second);
                            best_candids_labs = std::make_pair(site_lab_l, site_lab_r);
                        } 


                        // if ((2*u.count()-1)/(edges_num+1) + score <  (2*u.count()-1)/(best_edges_num+1) + best_mig_score) {
                        //     best_mig_score = score;
                        //     best_edges_num = edges_num;
                        //     best_candids = std::make_pair(children_tax.first,children_tax.second);
                        //     best_candids_labs = std::make_pair(site_lab_l, site_lab_r);
                        // } 
                        // else if (score == best_mig_score) {
                        //     if (best_edges_num < edges_num) {
                        //         best_edges_num = edges_num;
                        //         best_candids = std::make_pair(children_tax.first,children_tax.second);
                        //         best_candids_labs = std::make_pair(site_lab_l, site_lab_r);
                        //     }
                        // }

                        }
                    }
                    
                }

                trace_back[std::make_pair(u, site_lab)] = std::make_pair(best_candids, best_candids_labs);
                mig_dp[std::make_pair(u, site_lab)] = std::make_pair(best_mig_score, best_edges_num);
            } 
        }
    
    }

    int final_score = INT_INF;
    int final_edges_num = 0;
    
        std::cout << " site_lab: " << primary_tumor << " " << mig_dp[std::make_pair(S,primary_tumor)].first;
    

    final_score = mig_dp[std::make_pair(S, primary_tumor)].first;
    return std::make_pair(final_score, trace_back);

}


Tree binary_mig_solution(std::unordered_map<std::pair<Clade, std::string>,  std::pair<std::pair<Clade, Clade>, std::pair<std::string, std::string>>,
                PairHash, PairEqual> &trace_back,
                std::unordered_map<Clade, cpp_rational> &freq,
                Clade &S,
                std::string &primary_tumor,
                std::vector<std::string>labels
                ) {
                    std::vector<Clade> clades;
                    clades.push_back(S);
                    std::queue<std::pair<Clade, std::string>> que;
                    que.push(std::make_pair(S, primary_tumor));
                    while (!que.empty()) {
                        Clade t = que.front().first;
                        std::string lab = que.front().second;
                        que.pop();
                        clades.push_back(t);
                        std::pair<Clade, std::string> parameter(t,lab);

                        if (trace_back.find(parameter) != trace_back.end()) {
                            Clade x = trace_back[parameter].first.first;
                            Clade y = trace_back[parameter].first.second;
                            std::string x_lab = trace_back[parameter].second.first;
                            std::string y_lab = trace_back[parameter].second.second;
                            que.push(std::make_pair(x,x_lab));
                            que.push(std::make_pair(y,y_lab));
                            }
                        }

                        std::cout << " # Clades = " << clades.size() << std::endl;
                        return build_tree_from_compat_clades(clades, labels, freq, true);
                }








