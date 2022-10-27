#include <bits/stdc++.h>
#include <random>
#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <lemon/dijkstra.h>
#include <lemon/adaptors.h>
#include <lemon/concepts/path.h>
#include <lemon/concepts/graph.h>
#include <lemon/lp.h>

using namespace lemon;
using namespace std;


#ifndef ROBUST_PATH
#define ROBUST_PATH
class Graph{
	public:
		int n_;
		int edge_number_;
		double erdos_edge_possible_;
		ListDigraph g;
		
		Graph(int n, double erdos_p);
};


class Paths : public Graph {
    public:
		int people_n_;
		vector<pair<int,int>> paths_; //Starting and Edning vertices of the people travelling
		ListDigraph::ArcMap<double> arc_cost_q_{g};
		ListDigraph::ArcMap<double> arc_buy_p_{g};
		vector<vector<vector<double>>> utilitys_;
		vector<vector<int>> defining_polyhedra_q_;
		vector<vector<int>> defining_polyhedra_u_;
		double max_earn_;
		
		Paths(int people, int n, double erdos_p);
		
		pair<int,int> RandomPath();
		
		void RandomPaths();
		
		//ListDigraph:: ShortestPathCost();
		
		void Perturbation();
		
		void upper_add_line(double min_l, double max_l, map<int, Mip::Col> &f, Mip &mip);

		void lower_add_line(double min_l, map<int, Mip::Col> &f, Mip &mip);
		
		void subset_upper_value(double prob_in_subset, double upper_value, vector<vector<int>> &defining_polyhedra, map<int, Lp::Col> &f, Lp &lp, std::ostream &os = std::cerr);

		void CreateUpperLowerPolyhedra(map<int, Mip::Col> &f, Mip &mip);
		
		void CreateSubsetPolyhedra(int many_ineq_q, double prob_in_subset, double upper_value, vector<vector<int>> &defining_polyhedra, map<int, Lp::Col> &f, Lp &lp, std::ostream &os = std::cerr);
		
		void SettingQValue(int many_ineq_q, double prob_q, double max_q, std::ostream &os = std::cerr);
		
		void CreatingUPolyhedra(int many_ineq_u, double prob_u, double max_u, std::ostream &os = std::cerr);
		
		void IpSolve34(vector<double> &q_tariff, int big_M);
		
		void IPSolve35(int big_M);
		
		void FindingOptimalCost(int many_ineq_q, double prob_q, double max_q, int many_ineq_u, double prob_u, double max_u, std::ostream &os = std::cerr);
		
		void PrintData(std::ostream &os = std::cerr);
    
    
};

#endif ROBUST_PATH

#ifndef UTILITY_TOOLS
#define UTILITY_TOOLS
using ll = long long int;

const long long int INF = 10000000000;
const bool DEBUG = true;

#define all(x) begin(x), end(x)
#define FOR(i,n) for(int i = 0; i < (n); ++i)
#define FORO(i,n) for(int i = 1; i < (n); ++i)

template <class C>
void Print_vector(const C &Original, std::ostream &os = std::cerr) {
	for(const auto &v : Original) {
	    os << v << " ";
	}
	os << endl;
}

template <class C>
void Print_Matrix(const vector<vector<C>> &M, std::ostream &os = std::cerr) {
	for(auto &v : M) {
		for(auto &u : v) {
			os << u; os << " ";            
			}
	    os << endl;
	}
}

template<class T, class C>
void Print_pair(const pair<T,C> &M, std::ostream &os = std::cerr) {
    os << "(" << M.first << " , " << M.second << " ) ";
}

template <class C>
void Print_vector_pairs(const C &Original, std::ostream &os = std::cerr) {
	for(const auto &v : Original) {
	    Print_pair(v, os);
	}
	os << endl;
}
#endif //UTILITY_TOOLS
