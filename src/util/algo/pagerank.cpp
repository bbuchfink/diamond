#include <vector>
#include <unordered_map>
//#include <numeric>
#include <iostream>
#include "algo.h"
//#define _REENTRANT
//#include "lib/ips4o/ips4o.hpp"
#include "util/log_stream.h"
#include "basic/config.h"
using namespace std;

// #define LOGGING

namespace Util { namespace Algo {

template<typename Int>
void sweep_cut_min_conductance(
    FlatArray<Int>& neighbors,
    const Int rep,
    const vector<double>& scores,
    vector<Int>& reps
) {
    int n = neighbors.size(), unclustered = 0;
    vector<int> deg(n, 0);
    vector<int> order;
    for (int i = 0; i < n; ++i) {
        if (reps[i] == -1 && scores[i] > 0.0) {
            deg[i] = (int)neighbors.count(i);
            order.push_back(i);
            ++unclustered;
        }
    }
    if (unclustered <= 1) {
        reps[rep] = rep;
        return;
    }
    
    sort(order.begin(), order.end(), [&](int a, int b) {
        double sa = (deg[a] > 0) ? scores[a] / deg[a] : 0.0;
        double sb = (deg[b] > 0) ? scores[b] / deg[b] : 0.0;
        if (sa != sb) return sa > sb;
        return a < b;
        });

    vector<char> inS(n, false);
    long long volS = 0;
    long long totalVol = 0;
    for (int d : deg) totalVol += d;
    long long cutS = 0;

    double best_phi = 1.0;
    int best_k = 0;

    for (int idx = 0; idx < unclustered; ++idx) {
        int u = order[idx];
        inS[u] = true;
        volS += deg[u];
        const auto end = neighbors.cend(u);
        for (auto i = neighbors.cbegin(u); i < end; ++i) {
            if (reps[*i] != -1 || scores[*i] <= 0.0)
                continue;
            if (inS[*i]) {
                cutS -= 1;
            }
            else {
                cutS += 1;
            }
        }
        
        if(volS > 0) {
            //double phi = (double)cutS / (double)min(volS, totalVol - volS);
            double phi = (double)cutS / (double)volS;
            if (phi <= best_phi) {
                best_phi = phi;
                best_k = idx + 1;
            }
#ifdef LOGGING
            cout << u << ' ' << phi << endl;
#endif
        }
    }

    for (int i = 0; i < best_k; ++i) {
        reps[order[i]] = rep;
    }
	reps[rep] = rep;
#ifdef LOGGING
	int rep_n = neighbors.count(rep);
    cout << "best_k=" << best_k << " " << " cc=" << unclustered << " rep_n=" << rep_n;
    if (best_k < unclustered)
        cout << " uneven";
    if (rep_n < unclustered)
        cout << " indirect";
    cout << endl;
#endif
}

template<typename Int>
void build_cluster(
    FlatArray<Int>& neighbors,
    const Int rep,
    const vector<double>& scores,
    vector<Int>& reps
) {
    double cutoff = scores[rep]; // *config.ppr_multiplier;
    vector<double> nb;
    auto end = neighbors.cend(rep);
    for (auto i = neighbors.cbegin(rep); i < end; ++i) {
        if (reps[*i] != -1)
            continue;
        nb.push_back(scores[*i]);
    }
    if (nb.empty()) {
        reps[rep] = rep;
        return;
    }
    sort(nb.begin(), nb.end(), greater<double>());
    Int cut = 1; // std::max((Int)round((double)nb.size() * config.ppr_neighbor_cover), (Int)1);
    cutoff = std::min(cutoff, nb[cut - 1]);
    //cutoff *= config.ppr_multiplier;
    for (Int i = 0; i < scores.size(); ++i) {
        if (reps[i] == -1 && scores[i] >= cutoff) {
            reps[i] = rep;
        }
    }
    reps[rep] = rep;
}

template<typename Int>
void personalized_pagerank(
    FlatArray<Int>& neighbors,
    Int rep,
    vector<Int>& reps,
    double alpha = 0.85,
    double eps = 1e-8,
    int max_iters = 200
) {
	int n = (int)neighbors.size(), unclustered = 0;
    
    vector<int> deg(n, 0);
    for (int i = 0; i < n; ++i) {
        if (reps[i] == -1) {
            deg[i] = (int)neighbors.count(i);
			++unclustered;
        }
    }

    if(unclustered <= 1) {
		reps[rep] = rep;
        return;
	}

    vector<double> r(n, 0.0);
    r[rep] = 1.0;
    vector<double> rnext(n);

    for (int iter = 0; iter < max_iters; ++iter) {
        fill(rnext.begin(), rnext.end(), 0.0);

        double dangling_mass = 0.0;
        for (int j = 0; j < n; ++j) {
            if (reps[j] != -1)
                continue;
            if (deg[j] == 0) dangling_mass += r[j];
            else {
                double contrib = alpha * r[j] / deg[j];
                auto end = neighbors.cend(j);
                for (auto i = neighbors.cbegin(j); i < end; ++i)
                    if (reps[*i] == -1)
                        rnext[*i] += contrib;
            }
        }

        double teleport_factor = (1.0 - alpha);
        rnext[rep] += teleport_factor;
        rnext[rep] += alpha * dangling_mass;
        
        double diff = 0.0;
        for (int i = 0; i < n; ++i) diff += fabs(rnext[i] - r[i]);

        r.swap(rnext);

        if (diff < eps) {
#ifdef LOGGING
            cerr << "Converged in " << iter+1 << " iterations, diff=" << diff << " rep=" << rep << "\n";
#endif
            break;
        }
    }

#ifdef LOGGING
    vector<pair<double, Int>> v;
    v.reserve(n);
    for (Int i = 0; i < n; ++i)
        if (r[i] > 0.0)
            v.emplace_back(r[i], i);
    std::sort(v.begin(), v.end());
    for (auto& i : v)
        cout << i.first << ' ' << i.second << endl;
#endif

    //sweep_cut_min_conductance<Int>(neighbors, rep, r, reps);
	build_cluster<Int>(neighbors, rep, r, reps);
}

template<typename Int>
Int pr(FlatArray<Int>& neighbors, vector<Int>& reps) {
    
    const Int n = neighbors.size();
    int unclustered = 0;
    for (Int rep : reps)
        if (rep == -1)
            ++unclustered;
    if(unclustered == 0) {
        return -1;
	}
#ifdef LOGGING
	cerr << "Unclustered nodes: " << unclustered << endl;
#endif
    
    const double damping = 0.85;
    const double tol = 1e-12;
    const int maxIter = 1000;

    vector<double> pr(n, 0.0), next_pr(n, 0.0);
    vector<int> outdeg(n, 0);
    for (int i = 0; i < n; ++i) {
        if (reps[i] == -1) {
            pr[i] = 1.0 / unclustered;
            outdeg[i] = (int)neighbors.count(i);
        }
    }

    for (int iter = 0; iter < maxIter; ++iter) {
        fill(next_pr.begin(), next_pr.end(), 0.0);

        double dangling_sum = 0.0;
        for (int i = 0; i < n; ++i) if (outdeg[i] == 0) dangling_sum += pr[i];

        for (int u = 0; u < n; ++u) {
            if (outdeg[u] == 0) continue;
            double share = pr[u] / outdeg[u];
            const auto end = neighbors.cend(u);
            for (auto i = neighbors.cbegin(u); i < end; ++i)
                if (reps[*i] == -1)
                    next_pr[*i] += share;
        }

        double teleport = (1.0 - damping) / n;
        double dangling_share = damping * dangling_sum / n;

        double diff = 0.0;
        for (int i = 0; i < n; ++i) {
            if (reps[i] != -1)
                continue;
            double val = damping * next_pr[i] + dangling_share + teleport;
            diff += fabs(val - pr[i]);
            pr[i] = val;
        }

        if (diff < tol) break;
    }

    Int r = -1;
	double best = std::numeric_limits<double>::max();
    for (Int i = 0; i < n; ++i)
        if (reps[i] == -1 && pr[i] < best) {
            r = i;
            best = pr[i];
        }
        //std::cout << pr[i] << std::endl;
#ifdef LOGGING
    cout << r << ' ' << best << ' ' << neighbors.count(r) << endl;
#endif
    assert(r >= 0 && r < n);
    return r;
}

template<typename Int>
pair<FlatArray<Int>, Int> find_cc(FlatArray<Edge<Int>>& neighbors) {
    DSU dsu(neighbors.size());
	typename vector<Edge<Int>>::const_iterator end = neighbors.global_cend();
    for(auto it = neighbors.global_cbegin(); it < end; ++it)
		dsu.unite(it->node1, it->node2);
	Int cc = 0;
    FlatArray<Int> ccs;
	vector<pair<Int, Int>> cc_map;
    cc_map.reserve(neighbors.size());
    for (Int i = 0; i < neighbors.size(); ++i) {
        Int rep = dsu.find(i);
		cc_map.emplace_back(rep, i);
        if (rep == i)
            ++cc;
    }
    ips4o::parallel::sort(cc_map.begin(), cc_map.end());
    ccs.next();
    Int nontrivial = 0;
    for (Int i = 0; i < neighbors.size(); ++i) {
        if (i > 0 && cc_map[i].first != cc_map[i - 1].first) {
            ccs.next();
            if(ccs.count(ccs.size()-2) > 1)
				++nontrivial;
        }
        ccs.push_back(cc_map[i].second);
    }
    if (ccs.count(ccs.size() - 1) > 1)
        ++nontrivial;
    assert(cc = ccs.size());
    return make_pair(ccs, nontrivial);
}

template<typename Int>
static pair<FlatArray<Int>, vector<Int>> get_cc(Int cc, FlatArray<Edge<Int>>& neighbors, const FlatArray<Int>& ccs) {
    FlatArray<Int> r;
    vector<Int> id_map;
    unordered_map<Int, Int> rev_map;
    id_map.reserve(ccs.count(cc));
	rev_map.reserve(ccs.count(cc));
    for (auto i = ccs.cbegin(cc); i < ccs.cend(cc); ++i) {
        Int n1 = *i, n1_id = id_map.size();
        id_map.push_back(n1);
		rev_map[n1] = n1_id;
    }
	for (auto i = ccs.cbegin(cc); i < ccs.cend(cc); ++i) {
        r.next();
        const Int n1 = *i;
        assert(r.size() - 1 == rev_map.at(n1));
        for (auto j = neighbors.cbegin(n1); j < neighbors.cend(n1); ++j) {
			r.push_back(rev_map.at(j->node2));
        }
    }
    assert(r.size() == ccs.count(cc));
    return make_pair(r, id_map);
}

template<typename Int>
vector<Int> cluster_pr(FlatArray<Edge<Int>>& neighbors) {
    TaskTimer timer("Finding connected components");
    Int nontrivial;
    atomic<Int> clusters(0);
    FlatArray<Int> ccs;
    tie(ccs, nontrivial) = find_cc(neighbors);
    timer.finish();
    message_stream << "Connected components: " << nontrivial << endl;
    timer.go("Computing clustering");
    //vector<pair<Int, Int>> clustering;
    //clustering.reserve(neighbors.size());
	vector<Int> out(neighbors.size(), -1); // -1 means unclustered

    const int num_threads = std::thread::hardware_concurrency();
    std::atomic<Int> cc_index(0);
    std::vector<std::thread> threads;
    srand(1337);

    auto worker = [&](int tid) {
        while (true) {
            Int i = cc_index.fetch_add(1, std::memory_order_relaxed);
            if (i >= ccs.size()) break;
            #ifdef LOGGING
            cerr << "Processing component " << i << " of " << ccs.size() << " size=" << ccs.count(i) << " tid=" << tid << endl;
            #endif
            if (ccs.count(i) <= 1) {
				assert(ccs.count(i) == 1);
                //for (auto j = ccs.cbegin(i); j < ccs.cend(i); ++j)
                    //thread_clustering[tid].emplace_back(*j, *j);
                out[*ccs.cbegin(i)] = *ccs.cbegin(i);
                continue;
            }
            FlatArray<Int> cc_neighbors;
            vector<Int> id_map;
            tie(cc_neighbors, id_map) = get_cc(i, neighbors, ccs);
            vector<Int> reps(cc_neighbors.size(), -1);
            FlatArray<Int> clustering; // = mcl(cc_neighbors);
            clusters += clustering.size();
            for (Int i = 0; i < clustering.size(); ++i) {
				auto end = clustering.cend(i);
                const Int rep = id_map[clustering.cbegin(i)[rand() % clustering.count(i)]];
                for(auto it = clustering.cbegin(i); it < end; ++it) {
                    out[id_map[*it]] = rep;
				}
            }
            /*for (;;) {
                const Int rep = pr(cc_neighbors, reps);
                if (rep == -1) {
                    break;
                }
                personalized_pagerank<Int>(cc_neighbors, rep, reps);
                ++clusters;
            }
            for (Int j = 0; j < reps.size(); ++j)
                thread_clustering[tid].emplace_back(id_map[j], id_map[reps[j]]);*/
        }
        };

    for (int t = 0; t < num_threads; ++t) {
        threads.emplace_back(worker, t);
    }
    for (auto& th : threads) th.join();
    timer.finish();

    message_stream << "Clusters found: " << clusters << endl;
    /*ips4o::parallel::sort(clustering.begin(), clustering.end());
    vector<Int> out;
    out.reserve(neighbors.size());
    for (const auto& p : clustering)
        out.push_back(p.second);*/
    return out;
}

template vector<int32_t> cluster_pr(FlatArray<Edge<int32_t>>&);
template vector<int64_t> cluster_pr(FlatArray<Edge<int64_t>>&);

}}