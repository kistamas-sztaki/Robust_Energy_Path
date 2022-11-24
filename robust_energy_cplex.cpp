#include <bits/stdc++.h>
#include <random>
#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <lemon/dijkstra.h>
#include <lemon/adaptors.h>
#include <lemon/concepts/path.h>
#include <lemon/concepts/graph.h>

#include <ilcplex/ilocplex.h>

#include "robust_energy_cplex.h"

using namespace lemon;
using namespace std;



Graph::Graph(int n, double erdos_p) : n_{n},  edge_number_{0}, erdos_edge_possible_{erdos_p} {
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

Graph::Graph(std::istream &is) : edge_number_{0} {
	//Input format of the graph
	//First line is the number of vertcies n
	// Next is n lines in the ith one the outvertex number and the edgelist of i
	is >> n_;
	vector<ListDigraph::Node> nodes;
	FOR(i,n_) {
		nodes.push_back(g.addNode());
	}

	FOR(i,n_) {
		int outvertex; is >> outvertex;
		edge_number_ += outvertex;
		FOR(j,outvertex) {
			int out; is >> out;
			g.addArc(nodes[i], nodes[out]);
		}
	}
}

Paths::Paths(std::istream &is) : Graph(is), leader_max_earn_{0} {
	//The first input is the arc_buy_p_
	// then peoples destination with first how many people are there
	// The polyhedra of Q and U
	// A polyhedra is written in the following way
	// First line is the number of inequalities
	// the ith line starts with the number of variables in the inequalities, then
	// the indexes of the q variables in the inequality  // the numbers are in the following form "alpha x", where alpha is the coefficient of x 
	// the final number on the line is the upper limit of the inequality
    #if _DEBUG
    cerr << "Graph is in" << endl;
	PrintData();
    #endif

	#if _DEBUG
    cerr << "arc_buy_p" << endl;
    #endif

	//arc_buy_p_
		FOR(i,edge_number_) {
			double cost_p; is >> cost_p;
			arc_buy_p_[g.arcFromId(i)] = cost_p;
		}

	#if _DEBUG
    cerr << "-------------PEOPLE'S PATHS----------" << endl;
    #endif

	//Peoples paths
		is >> people_n_;
		FOR(i,people_n_) {
			int t1, t2; is >> t1 >> t2;
			paths_.push_back(make_pair(t1, t2));
		}

	#if _DEBUG
    cerr << "-------------Q POLYHEDRA READING IN-------------" << endl;
    #endif

	//Q
		int lines; is >> lines;
        //IloNumVarArray q(env); 
        FOR(_i,edge_number_) q_.add(IloNumVar(env, 0, +IloInfinity, ILOFLOAT));
		FOR(i,lines) {
			int variab; is >> variab;
            IloExpr expr(env);
			FOR(_j,variab-1) {
				int t; is >> t;
				expr += q_[t];
			}
            //ILOFLOAT;
            double maxi; is >> maxi;
            polyhedra_q_.add(expr <= maxi);
            expr.end();
		}

	#if _DEBUG
	IloCplex cplex(polyhedra_q_);
	cplex.exportModel("PROBLEM_Q.lp");
	cplex.end();
    cerr << "-------------U POLYHEDRA READING IN-------------" << endl;
	PrintData();
    #endif

	//U
		is >> lines; u_.setSize(people_n_);
		FOR(i,people_n_) {
			u_[i] = IloNumVarArray(env, edge_number_);
			FOR(j,edge_number_) {
				u_[i][j] = IloNumVar(env, 0., +IloInfinity, ILOFLOAT);
			}
		}
		
		FOR(i,lines) {
			int variab; is >> variab;
            IloExpr expr(env);
			FOR(_j,variab-1) {
				int t; is >> t;
				int remaind = t % edge_number_;
				int k = (t-remaind)/edge_number_;
				expr += u_[k][remaind];
			}
            double maxi; is >> maxi;
			polyhedra_u_.add(expr <= maxi);
            expr.end();
		}
		
	#if _DEBUG
	IloCplex cplex2(polyhedra_u_);
	cplex2.exportModel("PROBLEM_U.lp");
	cplex2.end();
    cerr << "-------------READ EVERYTHING IN-------------" << endl;
	PrintData();
    #endif
}

Paths::~Paths() {
	q_.end(); polyhedra_q_.end();
	u_.end(); polyhedra_u_.end();
    env.end();
}

void Paths::InitialQValue(std::ostream &os) {
	//Getting a q tariff that satisfies the constraints

	IloExpr expr(env); FOR(i,q_.getSize()) expr += q_[i];
	IloObjective obj(env, expr, IloObjective::Maximize);

	IloModel poly_q(env); poly_q.add(polyhedra_q_); poly_q.add(obj);
	IloCplex cplex_q(poly_q);
	#if _DEBUG
    //cplex_q.exportModel("PROBLEM_InitialQValue.lp");
    #endif
	cplex_q.solve();
	if(cplex_q.getStatus() == IloAlgorithm::Optimal) {
        IloNumArray qr(env); cplex_q.getValues(qr, q_);
		#if _DEBUG
		os << "The solution is for a sample q value:" << qr << endl;
		#endif
		FOR(i,qr.getSize()) {
			arc_cost_q_[g.arcFromId(i)] = qr[i];
		}
		qr.end();
    }
	else{os << "CAN NOT FIND SOLUTION FOR SettingQValue():" << endl;}

	expr.end();
	obj.end();
	poly_q.end();
	cplex_q.end();
}

double Paths::MinimizeLeadersEarning(const vector<double> &q_tariff, const int big_M, std::ostream &os) {
	//Given the tariff's on the roads, it gives the worst case for the leader.
	//Rerturn's the leader's minimal earning
	#if _DEBUG
    cerr << "\n\n-------MinimizeLeadersEarning-------" << endl;
    #endif

	IloModel model(env);
	// u \in U
		model.add(polyhedra_u_);
	#if _DEBUG
    cerr << "-------IMPORTED POLYHEDRA U-------" << endl;
    #endif

	NumVarMatrix x(env, people_n_);
	FOR(i,people_n_) {
		x[i] = IloNumVarArray(env, edge_number_);
		FOR(j,edge_number_) {
			x[i][j] = IloNumVar(env, 0., 1., ILOFLOAT);
		}
	}

	#if _DEBUG
    cerr << "-------DEFINING X-------" << endl;
	/*
	IloCplex cplex1(model);
	cplex1.exportModel("PROBLEM_1.lp");
	cplex1.end();
    cerr << "-------------PROBLEM1-------------" << endl;*/
    #endif

	IloNumVarArray alpha_plus(env, people_n_);
	IloNumVarArray alpha_negative(env, people_n_);
	FOR(i,alpha_plus.getSize()) {
		alpha_plus[i] = IloNumVar(env, 0., +INFINITY, ILOFLOAT);
		alpha_negative[i] = IloNumVar(env, 0., +INFINITY, ILOFLOAT);

		//a_plus >= 0 and x_j(\rho(t_j)) >= 1
			IloExpr expr(env);
			//x_j(\rho(t_j)) >= 1
			for (ListDigraph::InArcIt a(g, g.nodeFromId(paths_[i].second)); a != INVALID; ++a) {
				expr += x[i][g.id(a)];
			}
			model.add(expr >= 1);
			
			IloNumVar helper(env, 0, 1, ILOINT);
			model.add(alpha_plus[i] <= helper*big_M);
			model.add(expr-1 <= (1-helper)*big_M);
			//helper.end();
			expr.end();
			



		//a_minus >= 0 and -x_j(\delta(t_j))+1 >= 0
			IloExpr expr_2(env);
			//-x_j(\delta(t_j))+1 >= 0
			for (ListDigraph::OutArcIt a(g, g.nodeFromId(paths_[i].first)); a != INVALID; ++a) {
				expr_2 -= x[i][g.id(a)];
			}
			model.add(expr_2+1 >= 0);
			
			IloNumVar helper_2(env, 0, 1, ILOINT);
			model.add(alpha_negative[i] <= helper_2*big_M);
			model.add(expr_2+1 <= (1-helper_2)*big_M);
			expr_2.end();
			

	}

	#if _DEBUG
    cerr << "-------ALPHAS DEFINED-------" << endl;
	
	IloCplex cplex_2(model);
	cplex_2.exportModel("PROBLEM_2.lp");
	cplex_2.end();
    cerr << "-------------PROBLEM_2-------------" << endl;
    #endif

	NumVarMatrix gamma(env, people_n_);
	FOR(i,people_n_) {
		gamma[i] = IloNumVarArray(env, edge_number_);
		FOR(j,edge_number_) {
			gamma[i][j] = IloNumVar(env, 0., +IloInfinity, ILOFLOAT);
			IloNumVar helper(env, 0, 1, ILOINT);
			model.add(gamma[i][j] <= helper*big_M); //\gamma_{j,e}
			IloExpr expr_t = 1 - x[i][j];
			model.add(expr_t <= (1-helper)*big_M); //-x_j(e)+1 \geq 0
		}
	}


	#if _DEBUG
    cerr << "-------GAMMA DEFINED-------" << endl;
	
	IloCplex cplex_3(model);
	cplex_3.exportModel("PROBLEM_3.lp");
	cplex_3.end();
    cerr << "-------------PROBLEM_3-------------" << endl;
    #endif
	
	vector<map<int,IloNumVar>> beta(people_n_);
	FOR(j,people_n_) {
		FOR(v,n_) {
			ListDigraph::Node curr_v = g.nodeFromId(v);
			if(paths_[j].first != v && paths_[j].second != v) {
				beta[j][v] = IloNumVar(env, 0., +IloInfinity, ILOFLOAT);
				//beta[i].insert({j,IloNumVar(env, 0., +IloInfinity, ILOFLOAT)});
				IloExpr expr(env);
				int is_there{0};
				for(ListDigraph::InArcIt a(g, curr_v); a != INVALID; ++a) {
					expr += x[j][g.id(a)];
					++is_there;
				}
				for (ListDigraph::OutArcIt a(g, curr_v); a != INVALID; ++a) {
					expr -= x[j][g.id(a)];
					++is_there;
				}
				if(is_there) {
					model.add(expr >= 0);

					IloNumVar helper(env, 0, 1, ILOINT);
					model.add(beta[j][v] <= helper*big_M);
					model.add(expr <= (1-helper)*big_M);
					//helper.end();
				}
				
				expr.end();
			}
		}
	}
	

	#if _DEBUG
    cerr << "-------BETA DEFINED-------" << endl;
	IloCplex cplex_4(model);
	cplex_4.exportModel("PROBLEM_4.lp");
	cplex_4.end();
    cerr << "-------------PROBLEM_4-------------" << endl;
    #endif

	FOR(j,people_n_) {
		FOR(e,edge_number_) {
			IloExpr expr(env);
			int source = g.id(g.source(g.arcFromId(e))); int target = g.id(g.target(g.arcFromId(e)));
			if(target == paths_[j].second) expr += alpha_plus[j];
			if(source == paths_[j].first) expr -= alpha_negative[j];
			expr -= gamma[j][e];
			//+ beta[j][g.id(g.target(g.arcFromId(e)))] - beta[j][g.id(g.source(g.arcFromId(e)))] 
			if(paths_[j].first != target && paths_[j].second != target) {
				expr += beta[j][target];
			}
			if(paths_[j].first != source && paths_[j].second != source) {
				expr -= beta[j][source];
			}
			model.add(expr == q_tariff[e]- u_[j][e]);
			expr.end();
		}
	}

	#if _DEBUG
    cerr << "-------DUAL OPTIMALITY DEFINED-------" << endl;
	IloCplex cplex_5(model);
	cplex_5.exportModel("PROBLEM_5.lp");
	cplex_5.end();
    cerr << "-------------PROBLEM_5-------------" << endl;
    #endif

	
	IloExpr expr_obj(env);
	FOR(i,people_n_) {
		FOR(j,edge_number_) {
			expr_obj += q_tariff[j]*x[i][j];
		}
	}
	IloObjective obj(env, expr_obj, IloObjective::Minimize);
	model.add(obj);

	#if _DEBUG
    cerr << "-------OBJECTIVE DEFINED-------" << endl;
    #endif


	IloCplex cplex(model);
	#if _DEBUG
    cplex.exportModel("PROBLEM_MinimizeLeadersEarning.lp");
    #endif

	cplex.solve();
	double obj_value{0};
	/*
	if(cplex.getStatus() == IloAlgorithm::Optimal) {
		try{
			NumMatrix ur(env); ur.setSize(people_n_);
			FOR(i,people_n_) {
				ur[i] = IloNumArray(env, edge_number_);
				cplex.getValues(ur[i], u_[i]);
			}

			if(true) {//TODO
				set_of_utilities_.push_back(ur);
			}

			#if _DEBUG
			os << "The solution for u is:\n";
			FOR(i,people_n_)
				os << ur[i] << endl;
			#endif

			obj_value = cplex.getObjValue();
		}
		catch(IloAlgorithm::NotExtractedException) {
			cerr << "IPSOLVE34 VARIABLES ARE NOT RELATED TO THE OBJECTIVE\n"; 
		}
    }
	*/
	switch (cplex.getStatus())
	{

		case IloAlgorithm::Optimal:
			try{
			NumMatrix ur(env); ur.setSize(people_n_);
			FOR(i,people_n_) {
				ur[i] = IloNumArray(env, edge_number_);
				cplex.getValues(ur[i], u_[i]);
			}

			if(true) {//TODO
				set_of_utilities_.push_back(ur);
			}

			#if _DEBUG
			os << "The solution for u is:\n";
			FOR(i,people_n_)
				os << ur[i] << endl;
			#endif

			obj_value = cplex.getObjValue();
			}
			catch(IloAlgorithm::NotExtractedException) {
				cerr << "IPSOLVE34 VARIABLES ARE NOT RELATED TO THE OBJECTIVE\n"; 
			}
			break;
	
		case IloAlgorithm::Unbounded:
			os << "THE PROBLEM IS UNBOUNDED:\n";
			break;
		case IloAlgorithm::Infeasible:
			os << "THE PROBLEM IS Infeasible:\n";
			break;
		
		case IloAlgorithm::Error:
			os << "THE PROBLEM HAS Error:\n";
			break;
	
	}

	expr_obj.end();
	cplex.end();
	x.end();
	alpha_plus.end();
	alpha_negative.end();
	model.end();

	return obj_value;
}

void Paths::FindingTariffWithFiniteUtilities(const int big_M, std::ostream &os) {
	//Given a discrete set of utility vectors (set_of_utilities_) it determines the leader's best response

	IloModel model(env);

	// q \in Q
		model.add(polyhedra_q_);

	IloNumVar Z(env, -IloInfinity, +IloInfinity, ILOFLOAT);

	NumVarMatrix x(env, people_n_);
	FOR(i,people_n_) {
		x[i] = IloNumVarArray(env, edge_number_);
		FOR(j,edge_number_) {
			x[i][j] = IloNumVar(env, 0., 1., ILOFLOAT);
		}
	}

	#if _DEBUG
    cerr << "-------DEFINING X-------" << endl;
    #endif

	IloNumVarArray alpha_plus(env, people_n_);
	IloNumVarArray alpha_negative(env, people_n_);
	FOR(i,alpha_plus.getSize()) {
		alpha_plus[i] = IloNumVar(env, 0., +INFINITY, ILOFLOAT);
		alpha_negative[i] = IloNumVar(env, 0., +INFINITY, ILOFLOAT);

		//a_plus >= 0 and x_j(\rho(t_j)) >= 1
			IloExpr expr(env);
			//x_j(\rho(t_j)) >= 1
			for (ListDigraph::InArcIt a(g, g.nodeFromId(paths_[i].second)); a != INVALID; ++a) {
				expr += x[i][g.id(a)];
			}
			model.add(expr >= 1);
			
			IloNumVar helper(env, 0, 1, ILOINT);
			model.add(alpha_plus[i] <= helper*big_M);
			model.add(expr-1 <= (1-helper)*big_M);
			expr.end();
			



		//a_minus >= 0 and -x_j(\delta(t_j))+1 >= 0
			IloExpr expr_2(env);
			//-x_j(\delta(t_j))+1 >= 0
			for (ListDigraph::OutArcIt a(g, g.nodeFromId(paths_[i].second)); a != INVALID; ++a) {
				expr_2 -= x[i][g.id(a)];
			}
			model.add(expr_2+1 >= 0);
			
			IloNumVar helper_2(env, 0, 1, ILOINT);
			model.add(alpha_negative[i] <= helper_2*big_M);
			model.add(expr_2+1 <= (1-helper_2)*big_M);
			expr_2.end();
			

	}

	#if _DEBUG
    cerr << "-------ALPHAS DEFINED-------" << endl;
    #endif

	NumVarMatrix gamma(env, people_n_);
	FOR(i,people_n_) {
		gamma[i] = IloNumVarArray(env, edge_number_);
		FOR(j,edge_number_) {
			gamma[i][j] = IloNumVar(env, 0., +IloInfinity, ILOFLOAT);
			IloNumVar helper(env, 0, 1, ILOINT);
			model.add(gamma[i][j] <= helper*big_M); //\gamma_{j,e}
			IloExpr expr_t = 1 - x[i][j];
			model.add(expr_t <= (1-helper)*big_M); //-x_j(e)+1 \geq 0
		}
	}


	#if _DEBUG
    cerr << "-------GAMMA DEFINED-------" << endl;
    #endif
	/*
	IloArray<map<int,IloNumVar>> beta(env, people_n_);
	FOR(i,people_n_) {
		FOR(j,n_) {
			if(paths_[i].first != j && paths_[i].second != j) {
				//beta[i][j] = IloNumVar(env, 0., +IloInfinity, ILOFLOAT);
				beta[i].insert({j,IloNumVar(env, 0., +IloInfinity, ILOFLOAT)});
				IloNumVar helper(env, 0, 1, ILOINT);
				//model.add
				
			}
		}
	}
	*/

	#if _DEBUG
    cerr << "-------BETA DEFINED-------" << endl;
    #endif

	FOR(j,people_n_) {
		FOR(e,edge_number_) {
			IloExpr expr(env);
			expr += alpha_plus[j]- alpha_negative[j] - gamma[j][e];
			/*Beta*/ 
			model.add(expr == q_[e]- u_[j][e]);
			expr.end();
		}
	}

	#if _DEBUG
    cerr << "-------DUAL OPTIMALITY DEFINED-------" << endl;
    #endif

	#if _DEBUG
    cerr << "-------CONSTRAINTS ON Z DEFINED-------" << endl;
    #endif


	IloCplex cplex(model);

	#if _DEBUG
    cplex.exportModel("problem_3_4.lp");
    #endif

	cplex.solve();
	if(cplex.getStatus() == IloAlgorithm::Optimal) {
		
		try{
			NumMatrix ur(env); ur.setSize(people_n_);
			FOR(i,people_n_) {
				ur[i] = IloNumArray(env, edge_number_);
				cplex.getValues(ur[i], u_[i]);
			}

			if(true) {//TODO
				set_of_utilities_.push_back(ur);
			}

			#if _DEBUG
			os << "The solution for u is:\n";
			FOR(i,people_n_)
				os << ur[i] << endl;
			#endif
		}
		catch(IloAlgorithm::NotExtractedException) {
			cerr << "IPSOLVE34 VARIABLES ARE NOT RELATED TO THE OBJECTIVE\n"; 
		}
    }

	Z.end();
	//expr_obj.end();
	cplex.end();
	x.end();
	alpha_plus.end();
	alpha_negative.end();
	model.end();
}

void Paths::FindingOptimalCost(std::ostream &os) {
	//Getting an initial q value saving it in arc_cost_q_
		InitialQValue();
	
	//Section 3.4. Solver Finding worst case for leader
		vector<double> q_tariff;
		for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
			q_tariff.push_back(arc_cost_q_[e]);
		}
		int big_M = 300;
		double leader_earn = MinimizeLeadersEarning(q_tariff, big_M); //hyperparam TODO to-tune
	
	//Section 3.5.
		//FindingTariffWithFiniteUtilities(big_M);
}

void Paths::PrintData(std::ostream &os) const {
	os << "\n\n-----------------PRINTDATA----------------" << endl;
    os << "The number of vertices: " << n_ << endl;
	os << "The number of arcs: " << edge_number_ << endl;
	os << "The Arcs: " << endl;
	for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
		os << "	(" << g.id(g.source(e)) << ", " << g.id(g.target(e)) << ")" << " the vendor is buying it originally at price: " << arc_buy_p_[e] << endl;
	}
	os << "Number of people travelling: " << people_n_ << endl;
	os << "Their starting and end points: " << endl;
	Print_vector_pairs(paths_, os);
	
	os << "The defining polyhedra for q:" << endl;
	os << polyhedra_q_ << endl;

	os << "The defining polyhedra for u:" << endl;
	os << polyhedra_u_ << endl;

	os << "Current set of utilities: " << endl;
	Print_vector(set_of_utilities_, os);

	os << "Leader's current maximum profit: " << leader_max_earn_ << endl;

	os << "-----------------PRINTDATA----------------" << endl;
}

void Paths::PrintDataRaw(std::ostream &os) const {
	os << n_ << endl;
	FOR(i,n_) {
		int out = 0;
		for (ListDigraph::OutArcIt a(g, g.nodeFromId(i)); a != INVALID; ++a) {
			++out;
		}
		os << out << " ";
		for (ListDigraph::OutArcIt a(g, g.nodeFromId(i)); a != INVALID; ++a) {
			os << g.id(g.target(a)) << " ";
		}
		os << endl;
	}
	for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
		os << arc_buy_p_[e] << " ";
	}
	os << people_n_ << endl;
	Print_vector_pairs_raw(paths_, os);
	/*Print_Matrix(defining_polyhedra_q_, os);
	Print_Matrix(defining_polyhedra_u_, os);*/
}
