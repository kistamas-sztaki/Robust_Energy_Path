#include <bits/stdc++.h>
#include <random>
#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <lemon/dijkstra.h>
#include <lemon/adaptors.h>
#include <lemon/concepts/path.h>
#include <lemon/concepts/graph.h>
#include <lemon/lp.h>

#include "robust_energy.h"

using namespace lemon;
using namespace std;

Graph::Graph(int n, double erdos_p) : n_{n}, erdos_edge_possible_{erdos_p}, edge_number_{0} {
		vector<ListDigraph::Node> nodes;
		FOR(i,n_) {
			nodes.push_back(g.addNode());
		}
		
		//Erdős-Renyi Graph Model
		std::random_device rd;  //Will be used to obtain a seed for the random number engine
		std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
		std::discrete_distribution<> distrib({1-erdos_p, erdos_p});
		FOR(i,n_) {
			FOR(j,n_) {
			//for(int j = i; j < n_; ++j) {
				if(distrib(gen) && i != j)
					{g.addArc(nodes[i], nodes[j]); ++edge_number_;}
			}
		}
		
}

Paths::Paths(int people, int n, double erdos_p): people_n_{people}, max_earn_{INF}, Graph{n, erdos_p} {
	RandomPaths();
}

pair<int,int> Paths::RandomPath() {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(1, n_-1);
    pair<int,int> r; r.first = distrib(gen);
    do{r.second = distrib(gen);}
    while(r.second == r.first);
    
    return r;
}

void Paths::RandomPaths() {
	FOR(i,people_n_) {
       paths_.push_back(RandomPath());
    }
}

void Paths::Perturbation() {
	vector<pair<int, ListDigraph::Arc>> pert;
	for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
		pert.push_back(make_pair(arc_cost_q_[e]-arc_buy_p_[e], e));
	}
}

void Paths::upper_add_line(double min_l, double max_l, map<int, Mip::Col> &f, Mip &mip) {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> distrib(min_l, max_l);
	for(auto v : f) {
		mip.colUpperBound(v.second, distrib(gen));
	}
}

void Paths::lower_add_line(double min_l, map<int, Mip::Col> &f, Mip &mip) {
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> distrib(0, min_l);
    for(auto v : f) {
		mip.colLowerBound(v.second, distrib(gen));
    }
}

void Paths::subset_upper_value(double prob_in_subset, double upper_value, vector<vector<int>> &defining_polyhedra, map<int, Lp::Col> &f, Lp &lp, std::ostream &os) {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::discrete_distribution<> distrib({1-prob_in_subset, prob_in_subset});
	
	//Creating the subset with at least one element in it
		vector<int> defining; //Storting the polyheadra for later use
		Lp::Expr expr;
		if(DEBUG) os << "	The followings are in the set: " << endl;
		std::uniform_int_distribution<> dist_alw(0, (int)f.size()-1);
		int in_the_expr = dist_alw(gen);
		defining.push_back(in_the_expr); //Storting the polyheadra for later use
		expr += f[in_the_expr];
		if(DEBUG) os << "		";
		for(auto v: f) {
			if(distrib(gen) && v.first != in_the_expr) {
				defining.push_back(v.first); //Storting the polyheadra for later use
				expr += v.second;
				if(DEBUG) os << v.first << " ";
			}
			if(v.first == in_the_expr) if(DEBUG) os << v.first << " ";
		}
		defining_polyhedra.push_back(defining); //Storting the polyheadra for later use
		if(DEBUG) os << endl;
	
	//Setting the upperlimit for the sum of subset
		std::uniform_real_distribution<> dist_m(0, upper_value);
		double upp_sub = dist_m(gen);
		if(DEBUG) os << "	The upper bound for the subset is: " << upp_sub << endl;
		lp.addRow(expr <= upp_sub);
}

void Paths::CreateUpperLowerPolyhedra(map<int, Mip::Col> &f, Mip &mip) {
	/*
    FOR(i,edge_number_) {
        upper_add_line(f);
        lower_add_line(f);
    }*/
}

void Paths::CreateSubsetPolyhedra(int many_ineq_q, double prob_in_subset, double upper_value, vector<vector<int>> &defining_polyhedra, map<int, Lp::Col> &f, Lp &lp, std::ostream &os) {
	FOR(i,many_ineq_q) {
		subset_upper_value(prob_in_subset, upper_value, defining_polyhedra, f, lp, os);
	}
}

void Paths::SettingQValue(int many_ineq_q, double prob_q, double max_q, std::ostream &os) {
	//Polyheadron creation
		Lp lp_q;
		map<int, Lp::Col> q_f;
		FOR(i,edge_number_) {q_f[i] = lp_q.addCol(); lp_q.colLowerBound(q_f[i], 0); lp_q.colUpperBound(q_f[i], max_q);}
		if(DEBUG) os << "Creating the polyheadra:\n";
		CreateSubsetPolyhedra(many_ineq_q, prob_q, max_q, defining_polyhedra_q_, q_f, lp_q, os);
		if(DEBUG) os << "Polyheadra for q is created\n";
	
	//Solving the Polyhedron creation
		lp_q.max();
		Lp::Expr expr; FOR(i,edge_number_) expr += q_f[i];
		lp_q.obj(expr);
		lp_q.solve();
		if (lp_q.primalType() == Lp::OPTIMAL || lp_q.primalType() == Lp::FEASIBLE)
			if(DEBUG) os << "lp_q is solved" << endl;
		else {if(DEBUG) os << "lp_q is unsolvable" << endl;}
		
	//Saving the solution to the class variable
		int i_e = 0;
		for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
			arc_cost_q_[e] = lp_q.primal(q_f[i_e++]);
		}
		
	//Printing out the q tariff
		if(DEBUG) {
			os << "A solution of q tariff" << endl;
			for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
				os << "(" << g.id(g.source(e)) << ", " << g.id(g.target(e)) << ") cost: " << arc_cost_q_[e] << endl;
			}
		}
}

void Paths::CreatingUPolyhedra(int many_ineq_u, double prob_u, double max_u, std::ostream &os) {
	Lp lp_u;
	map<int, Lp::Col> u_f; FOR(i,edge_number_) {u_f[i] = lp_u.addCol(); lp_u.addRow(u_f[i] >= 0);}
	if(DEBUG) os << "Creating the polyheadra:\n";
	CreateSubsetPolyhedra(many_ineq_u, prob_u, max_u, defining_polyhedra_u_, u_f, lp_u, os);
	if(DEBUG) os << "Polyheadra for u is created\n";
}

void Paths::IpSolve34(vector<double> &q_tariff, int big_M) {
	Mip mip;
	vector<map<int, Mip::Col>> u;
	vector<map<int, Mip::Col>> x; x.resize(people_n_);
	map<int, Mip::Col> a_plus;
	map<int, Mip::Col> a_minus;
	vector<map<ListDigraph::Node, Mip::Col>> beta;
	vector<map<int, Mip::Col>> gamma;
	
	FOR(j,people_n_) {
		a_plus[j] = mip.addCol(); mip.colLowerBound(a_plus, 0);
		a_minus[j] = mip.addCol(); mip.colLowerBound(a_minus, 0);
		FOR(i,edge_number_) {
			u[j][i] = mip.addCol();
			x[j][i] = mip.addCol();
			gamma[j][i] = mip.addCol(); mip.colLowerBound(gamma[j][i], 0);
			
			Mip::Col helper = mip.addCol(); mip.colType(helper, Mip::INTEGER); mip.colLowerBound(helper, 0); mip.colUpperBound(helper, 1);
			mip.addRow(gamma[j][i] <= helper*big_M); //\gamma_{j,e}
			Mip::Expr expr_t = 1 -x[j][i];
			mip.addRow(expr_t <= (1-helper)*big_M); //-x_j(e)+1 \geq 0
		}
	}
	
	//a_plus >= 0 and x_j(\rho(t_j)) >= 1
		FOR(j,people_n_) {
			Mip::Expr expr;
			for (ListDigraph::InArcIt a(g, g.nodeFromId(paths_[j].second)); a != INVALID; ++a) {
				expr += x[j][g.id(a)];
			}
			mip.addRow(expr-1 >= 0);
			//a_plus or  x_j(\rho(t_j)) -1  is zero
				Mip::Col helper = mip.addCol(); mip.colType(helper, Mip::INTEGER); mip.colLowerBound(helper, 0); mip.colUpperBound(helper, 1);
				mip.addRow(a_plus[j] <= helper*big_M); //a_plus
				
				mip.addRow(expr-1 <= (1-helper)*big_M); //x_j(\rho(t_j)) -1
		}
	
	//a_minus >= 0 and -x_j(\delta(t_j))+1 >= 0
		FOR(j,people_n_) {
			Mip::Expr expr;
			for (ListDigraph::OutArcIt a(g, g.nodeFromId(paths_[j].second)); a != INVALID; ++a) {
				expr -= x[j][g.id(a)];
			}
			mip.addRow(expr+1 >= 0);
			//a_minus or  -x_j(\delta(t_j))+1 >= 0  is zero
				Mip::Col helper = mip.addCol(); mip.colType(helper, Mip::INTEGER); mip.colLowerBound(helper, 0); mip.colUpperBound(helper, 1);
				mip.addRow(a_minus[j] <= helper*big_M); //a_minus
				
				mip.addRow(expr+1 <= (1-helper)*big_M); //-x_j(\delta(t_j))+1
		}
	
	//beta and flow conversing Segmentation fault from this
		FOR(j,people_n_) {
			Mip::Expr expr;
			for (ListDigraph::NodeIt u(g); u != INVALID; ++u) {
				if(g.nodeFromId(paths_[j].first) != u && g.nodeFromId(paths_[j].second) != u) {
					beta[j][u] = mip.addCol(); mip.colLowerBound(beta[j][u], 0);  // fault from this fault from this
					
					for(ListDigraph::InArcIt a(g, u); a != INVALID; ++a) {
						expr += x[j][g.id(a)];
					}
					for (ListDigraph::OutArcIt a(g, u); a != INVALID; ++a) {
						expr -= x[j][g.id(a)];
					}
					mip.addRow(expr >= 0);
					
					//Zero equality
					Mip::Col helper = mip.addCol(); mip.colType(helper, Mip::INTEGER); mip.colLowerBound(helper, 0); mip.colUpperBound(helper, 1);
					mip.addRow(beta[j][u] <= helper*big_M); //beta
				
					mip.addRow(expr <= (1-helper)*big_M); // x_j(\rho(v_j))-x_j(\delta(v_j)) \geq 0
					
				}
			}
		}
	
	//equality at the end
		FOR(j,people_n_) {
			FOR(i,edge_number_) {
				Mip::Expr expr;
				expr += a_plus[j] - a_minus[j] - gamma[j][i];
				ListDigraph::Arc e = g.arcFromId(i);
				if(g.target(e) != g.nodeFromId(paths_[j].first) && g.target(e) != g.nodeFromId(paths_[j].second)) {
					expr += beta[j][g.target(e)];
				}
				if(g.source(e) != g.nodeFromId(paths_[j].first) && g.source(e) != g.nodeFromId(paths_[j].second)) {
					expr -= beta[j][g.source(e)];
				}
				mip.addRow(expr-q_tariff[i] + u[j][i]  == 0);
			}
		}
	
	mip.min();
	Mip::Expr expr; FOR(j,people_n_) FOR(i,edge_number_) expr += q_tariff[i]*x[j][i];
	mip.obj(expr);
	mip.solve();
	
}

void Paths::IPSolve35(int big_M) {
	Mip mip;
	vector<map<int, Mip::Col>> u;
	vector<map<int, Mip::Col>> x; x.resize(people_n_);
	map<int, Mip::Col> a_plus;
	map<int, Mip::Col> a_minus;
	vector<map<ListDigraph::Node, Mip::Col>> beta;
	vector<map<int, Mip::Col>> gamma;
	
	FOR(j,people_n_) {
		a_plus[j] = mip.addCol(); mip.colLowerBound(a_plus, 0);
		a_minus[j] = mip.addCol(); mip.colLowerBound(a_minus, 0);
		FOR(i,edge_number_) {
			u[j][i] = mip.addCol();
			x[j][i] = mip.addCol();
			gamma[j][i] = mip.addCol(); mip.colLowerBound(gamma[j][i], 0);
			
			Mip::Col helper = mip.addCol(); mip.colType(helper, Mip::INTEGER); mip.colLowerBound(helper, 0); mip.colUpperBound(helper, 1);
			mip.addRow(gamma[j][i] <= helper*big_M); //\gamma_{j,e}
			Mip::Expr expr_t = 1 -x[j][i];
			mip.addRow(expr_t <= (1-helper)*big_M); //-x_j(e)+1 \geq 0
		}
	}
	
	//a_plus >= 0 and x_j(\rho(t_j)) >= 1
		FOR(j,people_n_) {
			Mip::Expr expr;
			for (ListDigraph::InArcIt a(g, g.nodeFromId(paths_[j].second)); a != INVALID; ++a) {
				expr += x[j][g.id(a)];
			}
			mip.addRow(expr-1 >= 0);
			//a_plus or  x_j(\rho(t_j)) -1  is zero
				Mip::Col helper = mip.addCol(); mip.colType(helper, Mip::INTEGER); mip.colLowerBound(helper, 0); mip.colUpperBound(helper, 1);
				mip.addRow(a_plus[j] <= helper*big_M); //a_plus
				
				mip.addRow(expr-1 <= (1-helper)*big_M); //x_j(\rho(t_j)) -1
		}
	
	//a_minus >= 0 and -x_j(\delta(t_j))+1 >= 0
		FOR(j,people_n_) {
			Mip::Expr expr;
			for (ListDigraph::OutArcIt a(g, g.nodeFromId(paths_[j].second)); a != INVALID; ++a) {
				expr -= x[j][g.id(a)];
			}
			mip.addRow(expr+1 >= 0);
			//a_minus or  -x_j(\delta(t_j))+1 >= 0  is zero
				Mip::Col helper = mip.addCol(); mip.colType(helper, Mip::INTEGER); mip.colLowerBound(helper, 0); mip.colUpperBound(helper, 1);
				mip.addRow(a_minus[j] <= helper*big_M); //a_minus
				
				mip.addRow(expr+1 <= (1-helper)*big_M); //-x_j(\delta(t_j))+1
		}
	
	//beta and flow conversing Segmentation fault from this
		FOR(j,people_n_) {
			Mip::Expr expr;
			for (ListDigraph::NodeIt u(g); u != INVALID; ++u) {
				if(g.nodeFromId(paths_[j].first) != u && g.nodeFromId(paths_[j].second) != u) {
					beta[j][u] = mip.addCol(); mip.colLowerBound(beta[j][u], 0);  // fault from this fault from this
					
					for(ListDigraph::InArcIt a(g, u); a != INVALID; ++a) {
						expr += x[j][g.id(a)];
					}
					for (ListDigraph::OutArcIt a(g, u); a != INVALID; ++a) {
						expr -= x[j][g.id(a)];
					}
					mip.addRow(expr >= 0);
					
					//Zero equality
					Mip::Col helper = mip.addCol(); mip.colType(helper, Mip::INTEGER); mip.colLowerBound(helper, 0); mip.colUpperBound(helper, 1);
					mip.addRow(beta[j][u] <= helper*big_M); //beta
				
					mip.addRow(expr <= (1-helper)*big_M); // x_j(\rho(v_j))-x_j(\delta(v_j)) \geq 0
					
				}
			}
		}
}

void Paths::FindingOptimalCost(int many_ineq_q, double prob_q, double max_q, int many_ineq_u, double prob_u, double max_u, std::ostream &os) {
	//Determining the polyhedron for q and getting a starting value
		SettingQValue(many_ineq_q, prob_q, max_q, os);
		
	//Creating the polyhedron for u
		CreatingUPolyhedra(many_ineq_u, prob_u, max_u, os);
		
	//Random p_buying value from original vendor
		double max_val_q = 10;
		std::random_device rd;  //Will be used to obtain a seed for the random number engine
		std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
		std::uniform_real_distribution<> distrib(0, max_val_q);
		for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
			std::uniform_real_distribution<> d_p(0, arc_cost_q_[e]);
			arc_buy_p_[e] = d_p(gen);
		}
	
	//Perturbation of q TODO
	//Section 3.4. Solver
		vector<double> q_tariff;
		for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
			q_tariff.push_back(arc_cost_q_[e]);
		}
		IpSolve34(q_tariff, 300); //hyperparam TODO to-tune
	//Section 3.5.
	//Section 3.6.
}

void Paths::PrintData(std::ostream &os) {
	os << "The number of vertices: " << n_ << endl;
	os << "The number of arcs: " << edge_number_ << endl;
	os << "The Arcs: " << endl;
	for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
		os << "	(" << g.id(g.source(e)) << ", " << g.id(g.target(e)) << ")" << " the vendor is buying it originally at price: " << arc_buy_p_[e] << endl;
	}
	os << "Number of people travelling: " << people_n_ << endl;
	os << "Their starting and end points: " << endl;
	Print_vector_pairs(paths_);
	
	os << "The defining polyhedra for q:" << endl;
	Print_Matrix(defining_polyhedra_q_);

	os << "The defining polyhedra for u:" << endl;
	Print_Matrix(defining_polyhedra_u_);
}

    
