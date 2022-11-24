#include <bits/stdc++.h>
#include <random>
#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <lemon/dijkstra.h>
#include <lemon/adaptors.h>
#include <lemon/concepts/path.h>
#include <lemon/concepts/graph.h>

#include <ilcplex/ilocplex.h>

using namespace lemon;
using namespace std;


#ifndef ROBUST_PATH
#define ROBUST_PATH

using NumVarMatrix = IloArray<IloNumVarArray>;
using NumMatrix = IloArray<IloNumArray>;
using NumMatrix3D = vector<NumMatrix>;

class Graph{
	protected:
		int n_;
		int edge_number_;
		double erdos_edge_possible_;
		ListDigraph g;
		
		Graph(int n, double erdos_p); //You read the input more specially D graph

		Graph(std::istream &is);
};


class Paths : public Graph {
    private:
        IloEnv env;
		int people_n_;
		vector<pair<int,int>> paths_; //Starting and Ending vertices of the people travelling
		ListDigraph::ArcMap<double> arc_cost_q_{g};
		ListDigraph::ArcMap<double> arc_buy_p_{g};

        IloModel polyhedra_q_{env}; IloNumVarArray q_{env};
        IloModel polyhedra_u_{env}; NumVarMatrix u_{env};
		NumMatrix3D set_of_utilities_;

		//vector<vector<int>> defining_polyhedra_q_;
		//vector<vector<int>> defining_polyhedra_u_;
		double leader_max_earn_;
		
		pair<int,int> RandomPath();
		
		void RandomPaths();
		
		
		void Perturbation();
		/*
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
		*/

		//void CreateSubsetPolyhedra(int many_ineq_q, double prob_in_subset, double upper_value, vector<vector<int>> &defining_polyhedra, std::ostream &os = std::cerr);

		//void CreatingUPolyhedra(int many_ineq_u, double prob_u, double max_u, std::ostream &os = std::cerr);

		void InitialQValue(std::ostream &os = std::cerr);

		double MinimizeLeadersEarning(const vector<double> &q_tariff, const int big_M, std::ostream &os = std::cerr);

		void FindingTariffWithFiniteUtilities(const int big_M, std::ostream &os = std::cerr);

	public:

		Paths(const int people, const int n, const double erdos_p);

		Paths(std::istream &is); //You read the input more specially D graph, paths of the people, arc_cost_buy_p, Q polyhedra, U polyhedra

        ~Paths();

		void GenerateProblem(const int seed);

		void FindingOptimalCost(std::ostream &os = std::cerr);

		void PrintData(std::ostream &os = std::cerr) const;
        
		void PrintDataRaw(std::ostream &os) const;

		void SaveGenerated(std::ostream &os);
    
    
};

#endif //ROBUST_PATH


#ifndef UTILITY_TOOLS
#define UTILITY_TOOLS
using ll = long long int;

const long long int INF = std::numeric_limits<long long int>::max();
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
		os << v.size() << " ";
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

template <class C>
void Print_vector_pairs_raw(const C &Original, std::ostream &os = std::cerr) {
	for(const auto &v : Original) {
	    os << v.first << " " << v.second << " ";
	}
	os << endl;
}
#endif //UTILITY_TOOLS
