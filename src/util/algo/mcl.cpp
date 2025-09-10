#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cmath>
#include <set>
#include "algo.h"
#include "basic/config.h"

using namespace std;
using Eigen::MatrixXd;

namespace Util { namespace Algo {

static void add_self_loops(MatrixXd& M, double loop_value) {
    int n = (int)M.rows();
    for (int i = 0; i < n; ++i) M(i, i) += loop_value;
}

static void normalize_columns(MatrixXd& M) {
    int n = (int)M.cols();
    for (int j = 0; j < n; ++j) {
        double s = M.col(j).sum();
        if (s > 0.0) M.col(j) /= s;
    }
}

static MatrixXd matrix_power(MatrixXd base, int power) {
    int n = (int)base.rows();
    MatrixXd result = MatrixXd::Identity(n, n);
    while (power > 0) {
        if (power & 1) result = result * base;
        base = base * base;
        power >>= 1;
    }
    return result;
}

static MatrixXd expand(const MatrixXd& M, int power) {
    if (power == 1) return M;
    return matrix_power(M, power);
}

static void inflate(MatrixXd& M, double inflation) {
    M = M.array().pow(inflation).matrix();
    normalize_columns(M);
}

static void prune(MatrixXd& M, double threshold) {
    for (int i = 0; i < M.rows(); ++i)
        for (int j = 0; j < M.cols(); ++j)
            if (M(i, j) < threshold) M(i, j) = 0.0;
}

static bool has_converged(const MatrixXd& A, const MatrixXd& B, double eps) {
    return ((A - B).cwiseAbs().maxCoeff() < eps);
}

template<typename Int>
static pair<FlatArray<Int>, Int> extract_clusters(const MatrixXd& M, double threshold = 1e-3) {
    int n = (int)M.rows();
#ifdef LOGGING
    cout << M << endl;
#endif
    DSU dsu(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j && M(i, j) >= threshold)
                dsu.unite(i, j);
    Int cc = 0;
    FlatArray<Int> ccs;
    vector<pair<Int, Int>> cc_map;
    cc_map.reserve(n);
    for (Int i = 0; i < n; ++i) {
        Int rep = dsu.find(i);
        cc_map.emplace_back(rep, i);
        if (rep == i)
            ++cc;
    }
    sort(cc_map.begin(), cc_map.end());
    ccs.next();
    Int nontrivial = 0;
    for (Int i = 0; i < n; ++i) {
        if (i > 0 && cc_map[i].first != cc_map[i - 1].first) {
            ccs.next();
            if (ccs.count(ccs.size() - 2) > 1)
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
static pair<FlatArray<Int>, Int> markov_clustering(MatrixXd M,
    double loop_value = 1.0,
    double prune_threshold = 1e-5,
    double conv_eps = 1e-5,
    int max_iter = 100) {
#ifdef LOGGING
    cout << M << endl;
#endif
    add_self_loops(M, loop_value);
    normalize_columns(M);
	
    MatrixXd prev = M;
    for (int iter = 0; iter < max_iter; ++iter) {
#ifdef LOGGING
		cerr << "Iteration " << (iter + 1) << "..." << endl;
#endif
        MatrixXd expanded = expand(M, config.cluster_mcl_expansion);
        inflate(expanded, config.cluster_mcl_inflation);
        //prune(expanded, prune_threshold);
        normalize_columns(expanded);

        if (has_converged(expanded, M, conv_eps)) {
            M = expanded;
#ifdef LOGGING
            cerr << "Converged at iteration " << (iter + 1) << endl;
#endif
                break;
        }
        M = expanded;
    }

    return extract_clusters<Int>(M, 1e-3);
}

template<typename Int>
static MatrixXd build_adjacency(const FlatArray<Int>& neighbors, double self_loop_weight = 0.0) {
	Int n = neighbors.size();
    MatrixXd A = MatrixXd::Zero(n, n);
    for (Int i = 0; i < n; ++i) {
        auto end = neighbors.cend(i);
        for (auto it = neighbors.cbegin(i); it < end; ++it) {
            int u = *it, v = i;
            if (u < 0 || u >= n || v < 0 || v >= n) throw runtime_error("");
            if (A(u, v) == 0.0) {
                A(u, v) += 1.0;
            }
            if (A(v, u) == 0.0) {
                A(v, u) += 1.0;
            }
        }
    }
    if (self_loop_weight != 0.0) for (int i = 0; i < n; ++i) A(i, i) += self_loop_weight;
    return A;
}

template<typename Int>
FlatArray<Int> mcl(
    FlatArray<Int>& neighbors) {

    MatrixXd A = build_adjacency(neighbors, 0.0);

    double loop_value = 1.0; // will be added inside markov_clustering
    double prune_threshold = 1e-5;
    double conv_eps = 1e-5;
    int max_iter = 200;

	FlatArray<Int> clusters;
    Int nontrivial;

    tie(clusters, nontrivial) = markov_clustering<Int>(A, loop_value, prune_threshold, conv_eps, max_iter);

#ifdef LOGGING
    cout << "Found " << clusters.size() << " clusters: ";
    for (size_t i = 0; i < clusters.size(); ++i) {
        cout << "Cluster " << i << ": ";
        for (auto it = clusters.cbegin(i); it < clusters.cend(i); ++it) cout << *it << " ";
        cout << ' ';
    }
#endif
    return clusters;
}

void mcl() {
    /*vector<mcxOption> options;
    mcxHash hashes;
    mclAlgParam param;
    mcxstatus status = mclAlgorithmInit(options.data(), &hashes, "", &param);*/
}

template FlatArray<int32_t> mcl(FlatArray<int32_t>& neighbors);
template FlatArray<int64_t> mcl(FlatArray<int64_t>& neighbors);

}}