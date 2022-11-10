#include <bits/stdc++.h>
#include <random>
#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <lemon/dijkstra.h>
#include <lemon/adaptors.h>
#include <lemon/concepts/path.h>
#include <lemon/concepts/graph.h>

#include "robust_energy.h"
#include "robust_energy.cpp"

using namespace std;
using namespace lemon;


int main() {
	/*
	ifstream fin("input_robust_path.txt");
	Paths Test(fin);
	fin.close();
	Test.PrintData();*/
	
	Paths Test(1, 10, 0.1);
	Test.PrintData();
	Test.FindingOptimalCost(4, 0.3, 5, 4, 0.2, 10, std::cerr);
	/*
	ofstream fout("out_robust_data.txt");
	Test.PrintData(fout);
	fout.close();
	*/
}
