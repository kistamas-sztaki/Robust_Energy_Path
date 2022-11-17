#include <bits/stdc++.h>
#include <random>
#include <ilcplex/ilocplex.h>

using namespace std;

#define all(x) begin(x), end(x)
#define FOR(i,n) for(int i = 0; i < (n); ++i)
#define FORO(i,n) for(int i = 1; i < (n); ++i)

using NumVarMatrix = IloArray<IloNumVarArray>;

int main() {
    IloEnv env;

    IloModel model_u(env); 
    int length = 4;
    int columns = 4;
    NumVarMatrix u_(env); u_.setSize(length);
    FOR(i,length) {
        u_[i] = IloNumVarArray(env, columns); //(env, 0, 0.3, ILOFLOAT);
        FOR(j,columns) {
            u_[i][j] = IloNumVar(env, 0, 0.3, ILOFLOAT);
        }
    }

    IloModel model(env); model.add(model_u);
    IloNumVarArray x(env);
    IloInt generators = 4;
    int max_l = 1;
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> distrib(0, max_l);
    for (IloInt j = 0; j < generators; ++j)
        x.add(IloSemiContVar(env, 0, distrib(gen)));
    
    IloExpr expr(env);
    std::uniform_int_distribution<> dist_i(1, 4);
    for (IloInt j = 0; j < generators; ++j) {
        expr += dist_i(gen)*x[j];
        cerr << "exprssion" << expr << endl;
    }
    FOR(i,length) 
        FOR(j,columns)
            expr += u_[i][j];

    IloObjective obj(env, expr, IloObjective::Maximize);
    //model.add(obj);
    model.add(x[0] + x[1] <= 0.5);
    cerr << "MODEL:\n" << model << endl;

    IloCplex cplex(model); cerr << "CPLEX:\n" << cplex.getModel() << endl;
    cplex.exportModel("problem.lp");
    model.add(obj);
    cplex.extract(model);
    cplex.exportModel("problem_2.lp");

    
    cplex.solve();
    
    if(cplex.getStatus() == IloAlgorithm::Optimal) {
        IloNumArray vr(env); cplex.getValues(vr, x);
        cout << "The solution for x:" << vr << endl;
        
        IloArray<IloNumArray> ur(env); ur.setSize(length); FOR(i,length) ur[i] = IloNumArray(env, columns);
        //cplex.getValues(ur[0], u_[0]);
        //cout << "The solution for u:" << ur[0] << endl;
        FOR(i,length) {
            cplex.getValues(ur[i], u_[i]);
            cout << "The solution for u:" << ur[i] << endl;
        }
    }

    obj.end();
    expr.end();
    x.end();
    model.end();
    env.end();
}
