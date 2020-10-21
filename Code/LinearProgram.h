/**
* \author Nicolas Forget
*
* This class describe an extreme point of the lower bound set. Used in the linear relaxation (LP relax).
*/

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
	void formateToMin(); // automatically switch objectives in maximization to min -z_k(x)

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

	void buildDualBenson(MathematicalFormulation* LP); // creates template of the dual benson LP
	void solveDualBenson(std::vector<double>& y); // always feasible
	std::vector<double> getNormalVector();
	double getRhs(MathematicalFormulation& LP);

	void buildBestValidPoint(MathematicalFormulation* LP); // creates template for finding furthest interior point
	void solveBestValidPoint(std::vector<double>& s, std::vector<double>& phat); // always feasible
	std::vector<double> getBestValidPoint(std::vector<double>& s, std::vector<double>& phat);

	void buildFeasibility(MathematicalFormulation* LP); // creates template for checking feasibility of an extreme point
	bool solveFeasibility(std::vector<double>& s); // true is solved to optimality, false otherwise
	void retrieveSolutionFeasibility(std::vector<double>& s);

	void buildWeightedSumScalarization(MathematicalFormulation* LP); // creates template for solving a weighted sum scalarization
	void solveWeightedSumScalarization(MathematicalFormulation& LP, std::vector<double> lambda); // [toDo] true is solved to optimality, false otherwise
	double getWeightedSumValue(MathematicalFormulation& LP, std::vector<double> lambda);
	

	IloCplex* get_cplex();
	IloEnv* get_env();

private:
	IloEnv env;
	IloModel model;
	IloCplex cplex;
	IloNumVarArray x;
	IloObjective ptrObj;
	IloRangeArray ptrCtes;
};
