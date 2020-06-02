#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <ilcplex\ilocplex.h>

class LinearProgram
{
protected:
	int p;
	int n;
	int m;
	std::vector<int> optDirection; // 1 if max, 0 if min

public:
	LinearProgram();
	LinearProgram(std::string file);

	int get_p();
	int get_n();
	int get_m();
	int get_objDir(int i);
};



class MathematicalFormulation : public LinearProgram
{
public:
	MathematicalFormulation (std::string file);

	void print_objective();
	void print_constraints();
	void formate_to_min();

	int get_objective(int i, int j);
	int get_constraint(int i, int j);
	int get_rhs(int j);
	int get_signCte(int j);

private:
	std::vector<std::vector<int>> objMatrix;
	std::vector<std::vector<int>> constrMatrix;
	std::vector<std::pair<int,int>> rhs; // first : sign cte; second : rhs
};



class CplexModel : public LinearProgram // put all cte in standard form before building ??
{
public:
	CplexModel();
	void build_dualBenson(MathematicalFormulation* LP); // creates template of the dual benson LP
	void build_bestValidPoint(MathematicalFormulation* LP); // creates template for finding furthest interior point
	void solveDualBenson(std::vector<double>& y);

private:
	IloEnv env;
	IloModel model;
	IloCplex cplex;
	IloNumVarArray x;
	IloObjective ptrObj;
};
