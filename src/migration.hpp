#ifndef MIGRATION
#define MIGRATION
#include "cc_lshp.hpp"

#endif 

typedef std::tuple<std::string, std::string> pp;



std::pair<std::pair<double,double>, std::unordered_map<std::pair<Clade, std::string>,  std::pair<std::pair<Clade, Clade>, std::pair<std::string, std::string>>,PairHash, PairEqual>> mig_dp(
    std::string &primary_tumor,
    std::unordered_map<Clade, std::string> &taxon2anatomical,
    std::vector<std::string> &anatomical_labels_vec,
    Clades_Set &Sigma,
    Clade &S,
    std::unordered_map<Clade,std::vector<int>> &clade2state,
    SIESTA &I,
    std::unordered_map<Bipartition, std::unordered_map<std::string, double>> &taxon_weight_meta,
    std::unordered_map<Bipartition, std::unordered_map<std::string, double>> &taxon_weight_reseeding

) {
    std::unordered_map<std::pair<Clade, std::string>, std::pair<int, int>,PairHash, PairEqual> mig_dp;

    std::unordered_map<std::pair<Clade, std::string>,  std::pair<std::pair<Clade, Clade>, std::pair<std::string, std::string>>,
                PairHash, PairEqual> trace_back;


    std::vector<Bipartition> Sigma_vec(Sigma.begin(), Sigma.end());
    
    sort(Sigma_vec.begin(), Sigma_vec.end(), [](Clade &s1, Clade &s2) {
        return s1.count() < s2.count();
        });


    for (auto &u : Sigma_vec) {
        for (auto &site_lab : anatomical_labels_vec) {
            
            if (u.count() == 1) {
                
                mig_dp[std::make_pair(u, site_lab)] = std::make_pair(taxon_weight_meta[u][site_lab], taxon_weight_reseeding[u][site_lab]);

            } else {
                
                    double best_meta = INF;
                    double best_reseeding = INF;
                    CladePair best_candids;
                    std::pair<std::string, std::string> best_candids_labs;
                    for (auto &children_tax : I[u]) {

                        for (auto &site_lab_l : anatomical_labels_vec) {

                            for (auto &site_lab_r : anatomical_labels_vec) {
                                double meta_1 = mig_dp[std::make_pair(children_tax.first,site_lab_l)].first;
                                double reseeding_1 = mig_dp[std::make_pair(children_tax.first, site_lab_l)].second;
                                double meta_2 = mig_dp[std::make_pair(children_tax.second, site_lab_r)].first;
                                double reseeding_2 = mig_dp[std::make_pair(children_tax.second, site_lab_r)].second;


                                double meta = meta_1 + meta_2;
                                

                                double reseeding = reseeding_1 + reseeding_2;

                                if (site_lab_l != site_lab) {
                                    if (site_lab_l == primary_tumor) {
                                        reseeding += 1;
                                    } else {
                                        meta += 1;
                                    }
                                }

                                if (site_lab_r != site_lab) {
                                    if (site_lab_r == primary_tumor) {
                                        reseeding += 1;
                                    } else {
                                        meta += 1;
                                    }
                                }

                                // 
                                if ((meta + reseeding < best_meta + best_reseeding) 
                                || (meta + reseeding == best_meta + best_reseeding && reseeding < best_reseeding) 
                                ) {
                                
                                best_meta = meta;
                                best_reseeding = reseeding;
                                
                                best_candids = std::make_pair(children_tax.first,children_tax.second);
                                best_candids_labs = std::make_pair(site_lab_l, site_lab_r);
                                }
                            }
                        }
                    
                    }

                    // std::cout << " cur best reseeding " << best_reseeding << std::endl;
                    trace_back[std::make_pair(u, site_lab)] = std::make_pair(best_candids, best_candids_labs);
                    mig_dp[std::make_pair(u, site_lab)] = std::make_pair(best_meta, best_reseeding);
                    
            } 
        }
    }

    double final_meta = mig_dp[std::make_pair(S,primary_tumor)].first;
    double final_reseeding =  mig_dp[std::make_pair(S,primary_tumor)].second;
    double final_score = final_meta + final_reseeding;
    
        std::cout << " site_lab: " << primary_tumor << " total: "  << final_score << " reseeding: " << mig_dp[std::make_pair(S,primary_tumor)].second << std::endl;
    
    return std::make_pair(std::make_pair(final_meta, final_reseeding), trace_back);

}


Tree binary_mig_solution(std::unordered_map<std::pair<Clade, std::string>,  std::pair<std::pair<Clade, Clade>, std::pair<std::string, std::string>>,
                PairHash, PairEqual> &trace_back,
                std::unordered_map<Clade, cpp_rational> &freq,
                Clade &S,
                std::string &primary_tumor,
                std::vector<std::string>labels,
                std::unordered_map<Bipartition, std::unordered_map<std::string, double>> &taxon_weight_meta,
                std::unordered_map<Bipartition, std::unordered_map<std::string, double>> &taxon_weight_reseeding
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








