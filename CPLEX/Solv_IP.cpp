#include <bits/stdc++.h>
#include <random>

#include <ilcplex/ilocplex.h>

using namespace std;

int main(int argc, char** argv) {
    IloEnv env;
    assert(argc >= 2);
    #if _DEBUG
    cerr << "THE INPUT ARGV: " << argv[1] << endl;
    #endif
    IloModel model(env);
    IloCplex cplex(env);
    IloObjective   obj;
    IloNumVarArray var(env);
    IloRangeArray  rng(env);
    cplex.importModel(model, argv[1], obj, var, rng);
    cplex.extract(model);
    cplex.solve();

    switch (cplex.getStatus())
	{
		case IloAlgorithm::Optimal:
			cerr << "THE PROBLEM HAS OPTIMUM:\n";
            cerr << "THE OPTIMUM OBJECTIVE VALUE IS: " << cplex.getObjValue() << endl;
			break;
	
		case IloAlgorithm::Unbounded:
			cerr << "THE PROBLEM IS UNBOUNDED:\n";
			break;
		case IloAlgorithm::Infeasible:
			cerr << "THE PROBLEM IS Infeasible:\n";
			break;
		
		case IloAlgorithm::Error:
			cerr << "THE PROBLEM HAS Error:\n";
			break;
	
	}

    cplex.end();
    env.end();
}