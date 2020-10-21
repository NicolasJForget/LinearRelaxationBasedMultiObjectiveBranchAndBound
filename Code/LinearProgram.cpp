#include "LinearProgram.h"

/* ==========================================================
		Constructors
 ========================================================= */

LinearProgram::LinearProgram() : p(0), n(0), m(0), optDirection(0) {}

LinearProgram::LinearProgram(std::string file) {

	n = 0;
	m = 0;
	p = 0;

	std::ifstream f(file);

	if (f) {

		std::string str;
		std::string line;

		// line 1

		getline(f, line);
		std::istringstream iss(line);
		iss >> n;
		iss >> m;
		iss >> p;

		// line 2 : optimization direction

		getline(f, line);
		getline(f, line);
		iss.clear();
		iss.str(line);
		std::vector<int> TempOptDirection(p);
		
		for (int i = 0; i < p; i++) {
			iss >> str;
			if (str == "maxsum") {
				TempOptDirection[i] = 1;
			}
			else if (str == "minsum") {
				TempOptDirection[i] = -1;
			}
		}
		optDirection = TempOptDirection;
	}
	else {
		std::cout << "Problem with file opening" << std::endl;
	}

	f.close();
};

MathematicalFormulation::MathematicalFormulation(std::string file) : LinearProgram(file), objMatrix(0), constrMatrix(0), rhs(0) {

	std::ifstream f(file);

	if (f) {

		std::string line;
		int i;
		std::string coef;
		getline(f, line);
		getline(f, line);
		getline(f, line);
		getline(f, line);
		std::istringstream iss(line);

		// objective coefficients matrix

		for (int k = 0; k < p; k++) {
			objMatrix.push_back(std::vector<int>(n));
			getline(f, line);
			iss.clear();
			iss.str(line);
			int i = 0;
			while (i < n) {
				iss >> coef;
				objMatrix[k][i] = stoi(coef);
				++i;
			}
		}

		// constraint coefficients matrix

		getline(f, line);
		for (int j = 0; j < m; j++) {
			constrMatrix.push_back(std::vector<int>(n));
			getline(f, line);
			iss.clear();
			iss.str(line);
			i = 0;
			while (i < n) {
				iss >> coef;
				constrMatrix[j][i] = stoi(coef);
				++i;
			}
		}

		// right-hand side of the constraints

		getline(f, line);
		for (int j = 0; j < m; j++) {
			rhs.push_back(std::pair<int,int>());
			getline(f, line);
			iss.clear();
			iss.str(line);
			iss >> coef;
			rhs[j].first = stoi(coef);
			iss >> coef;
			rhs[j].second = stoi(coef);
		}
	}
	else {
		std::cout << "Problem with file opening 2" << std::endl;
	}
};

CplexModel::CplexModel() : LinearProgram() {
	model = IloModel(env);
	cplex = IloCplex(model);
	x = IloNumVarArray(env);
	ptrObj = IloObjective(env);
	ptrCtes = IloRangeArray(env);
	//ptrCtes = std::vector<IloConstraint>(p);

	cplex.setOut(env.getNullStream());
}

/* ==========================================================
		Regular Methods
 ========================================================= */

void MathematicalFormulation::print_objective() {

	std::cout << "\n\n";
	for (int k = 0; k < p; k++) {
		std::cout << "objective " << k + 1 << " : ";
		for (int i = 0; i < n; i++) {
			std::cout << objMatrix[k][i] << " ";
		}
		std::cout << "\n";
	}
}

void MathematicalFormulation::print_constraints() {

	std::cout << "\n\n";
	for (int k = 0; k < m; k++) {
		std::cout << "constraint " << k + 1 << " : ";
		for (int i = 0; i < n; i++) {
			std::cout << constrMatrix[k][i] << " ";
		}
		if (rhs[k].first == 1) {
			std::cout << "  <=   ";
		}
		else if (rhs[k].first == 2) {
			std::cout << "  =   ";
		}
		else {
			std::cout << "  >=   ";
		}
		std::cout << rhs[k].second << "\n";
	}
}

void MathematicalFormulation::formateToMin() {

	for (int k = 0; k < p; k++) {
		if (optDirection[k] == 1) {
			for (int i = 0; i < n; i++) {
				objMatrix[k][i] *= -1;
			}
		}
	}

}


// -------------- Dual Benson -----------------

void CplexModel::buildDualBenson(MathematicalFormulation* LP) {

	try{
		// building LP caracteristics

		p = LP->get_p();
		n = LP->get_m() + LP->get_n() + LP->get_p(); //LP->get_m() + LP->get_p()
		m = LP->get_n() + 1;
		optDirection.push_back(1);

		// building variables

		for (int i = 0; i < LP->get_m(); i++) { // u variables

			// take care of sign of constraints in primal for sign of variables in dual

			if (LP->get_signCte(i) == 0) x.add(IloNumVar(env, 0, IloInfinity, ILOFLOAT)); // >= cte
			else if(LP->get_signCte(i) == 1) x.add(IloNumVar(env, -IloInfinity, 0, ILOFLOAT)); // <= cte
			else x.add(IloNumVar(env, -IloInfinity, IloInfinity, ILOFLOAT)); // = cte
		
		}
		for (int i = 0; i < LP->get_n(); i++) { // v variables (for UB ctes on variables)
			x.add(IloNumVar(env, 0, IloInfinity, ILOFLOAT));
		}
		for (int i = 0; i < LP->get_p(); i++) { // w variables 
			x.add(IloNumVar(env, 0, IloInfinity, ILOFLOAT));
		}

		// building objective function
	
		IloExpr obj(env);
		for (int i = 0; i < LP->get_m(); i++) { // u variables
			obj += LP->get_rhs(i) * x[i];
		}
		for (int i = 0; i < LP->get_n(); i++) { // v variables (UB on variables)
			obj -= 1.0 * x[LP->get_m() + i];
		}
		for (int i = 0; i < LP->get_p(); i++) { // w variables
			obj -= 1.0 * x[LP->get_m() + LP->get_n() + i]; // replace 1 with value of y when using template
		}

		ptrObj = IloObjective(env, obj, IloObjective::Maximize);
		model.add(ptrObj);
		obj.end();

		// building constraints

		std::vector<IloExpr> lhs(0);
		for (int i = 0; i < LP->get_n(); i++) { // 1st set of ctes
			lhs.push_back(IloExpr(env));
			for (int j = 0; j < LP->get_m(); j++) { // A^T.u
			lhs[i] += LP->get_constraint(j, i) * x[j];
			}
			lhs[i] -= 1.0 * x[LP->get_m() + i]; // - v
			for (int j = 0; j < LP->get_p(); j++) { // - C^T.w
				lhs[i] -= (double) LP->get_objective(j, i) * x[LP->get_m() + LP->get_n() + j];
			}
			model.add(lhs[i] <= 0);
		}

		//std::cout << LP->get_p() << "we are here\n";

		for (int i = 0; i < LP->get_p(); i++) {
			//std::cout << i << std::endl;
			lhs.push_back(IloExpr(env));
			lhs[m - 1] += x[LP->get_m() + LP->get_n() + i];
		}
		//std::cout << "what about here?\n";
		model.add(lhs[m - 1] == 1);
		//std::cout << "mmmmmyeah?\n";
	}
	catch (IloException& e)
	{
		std::cout << e << std::endl;
		e.end();
	}
}

void CplexModel::solveDualBenson(std::vector<double>& y) {

	try
	{
		// Modify objective coefficients
		for (int k = 0; k < p; k++) {
			ptrObj.setLinearCoef(x[n - p + k], -y[k]);
		}

		//cplex.exportModel("debug.lp");
		// solve
		cplex.solve();
		//std::cout << cplex.getStatus() << std::endl;

		// Retrieve solution
		//IloNumArray rlt(env);
		//cplex.getValues(rlt, x);

		// print solution
		//std::cout << "\n\n\n";
		//for (int k = 0; k < n; k++) {
		//	std::cout << "x[" << k << "] = " << rlt[k] << std::endl;
		//}
		//std::cout << "\n\n\n";
	}
	catch (IloException& e)
	{
		std::cout << e << std::endl;
		e.end();
	}

	//cplex.exportModel("check2.lp");

}

std::vector<double> CplexModel::getNormalVector() {

	std::vector<double> theVector(p);

	IloNumArray rlt(env);
	cplex.getValues(rlt, x);

	for (int k = 0; k < p; k++) {
		theVector[k] = rlt[n - p + k];
	}

	return theVector;
}

double CplexModel::getRhs(MathematicalFormulation& LP) {

	double theRightHandSide = 0;

	IloNumArray rlt(env);
	cplex.getValues(rlt, x);

	for (int i = 0; i < LP.get_m(); i++) {
		theRightHandSide += LP.get_rhs(i) * rlt[i];
	}
	for (int i = 0; i < LP.get_n(); i++) {
		theRightHandSide -= rlt[LP.get_m() + i];
	}

	return theRightHandSide;
}


// -------------- Best Valid Point --------------

void CplexModel::buildBestValidPoint(MathematicalFormulation* LP) {

	try{
		// building LP caracteristics

		p = LP->get_p();
		n = LP->get_n() + 1;
		m = LP->get_m() + LP->get_p();
		optDirection.push_back(1);

		// building variables

		for (int i = 0; i < LP->get_n(); i++) {
			x.add(IloNumVar(env, 0, 1, ILOFLOAT));
		}
		x.add(IloNumVar(env, 0, 1, ILOFLOAT));

		// builind objective function

		model.add(IloMaximize(env, x[n - 1]));

		// building constraints

		std::vector<IloExpr> lhs(0);
		for (int i = 0; i < LP->get_m(); i++) {
			lhs.push_back(IloExpr(env));
			for (int j = 0; j < LP->get_n(); j++) {
				lhs[i] += LP->get_constraint(i, j) * x[j];
			}
			//model.add(lhs[i] >= LP->get_rhs(i));
			if (LP->get_signCte(i) == 0) model.add(lhs[i] >= LP->get_rhs(i)); // >= cte
			else if (LP->get_signCte(i) == 1) model.add(lhs[i] <= LP->get_rhs(i)); // <= cte
			else model.add(lhs[i] == LP->get_rhs(i)); // = cte
		}

		for (int i = 0; i < LP->get_p(); i++) {
			ptrCtes.add(IloRange(env, -IloInfinity, 0)); // change rhs (p)
			for (int j = 0; j < LP->get_n(); j++) {
				ptrCtes[i].setLinearCoef(x[j], LP->get_objective(i, j));
			}
			ptrCtes[i].setLinearCoef(x[n - 1], 1); // change coef of last var (lambda)
			model.add(ptrCtes[i]);
		}
	}
	catch (IloException& e)
	{
		std::cout << e << std::endl;
		e.end();
	}

}

void CplexModel::solveBestValidPoint(std::vector<double>& s, std::vector<double>& phat) {

	try{
		// adjust constraints
		for (int i = 0; i < p; i++) {
			ptrCtes[i].setUB(phat[i]);
			ptrCtes[i].setLinearCoef(x[n - 1], phat[i] - s[i]);
		}

		// solve
		cplex.solve();

		// retreive results
		//IloNumArray rlt(env);
		//cplex.getValues(rlt, x);

		// print solution
		//for (int k = 0; k < n; k++) {
		//	std::cout << "x[" << k << "] = " << rlt[k] << std::endl;
		//}
	}
	catch (IloException& e)
	{
		std::cout << e << std::endl;
		e.end();
	}

	//cplex.exportModel("check3.lp");
}

std::vector<double> CplexModel::getBestValidPoint(std::vector<double>& s, std::vector<double>& phat) {

	solveBestValidPoint(s, phat);

	//IloNumArray rlt(env);
	//cplex.getValues(rlt, x);

	double lambda = cplex.getObjValue(); //rlt[n - 1];
	std::vector<double> thePoint(p);
	for (int k = 0; k < p; k++) {
		thePoint[k] = lambda * s[k] + (1 - lambda) * phat[k];
	}

	return thePoint;
}


// -------------- Feasibility Checks --------------

void CplexModel::buildFeasibility(MathematicalFormulation* LP) {

	// building LP caracteristics

	p = LP->get_p();
	n = LP->get_n();
	m = LP->get_m() + LP->get_p();

	// building variables

	for (int i = 0; i < LP->get_n(); i++) {
		x.add(IloNumVar(env, 0, 1, ILOFLOAT));
	}

	// building objective function

	model.add(IloMaximize(env, 0));

	// building constraints

	std::vector<IloExpr> lhs(0);
	for (int i = 0; i < LP->get_m(); i++) {
		lhs.push_back(IloExpr(env));
		for (int j = 0; j < LP->get_n(); j++) {
			lhs[i] += LP->get_constraint(i, j) * x[j];
		}

		if (LP->get_signCte(i) == 0) model.add(lhs[i] >= LP->get_rhs(i)); // >= cte
		else if (LP->get_signCte(i) == 1) model.add(lhs[i] <= LP->get_rhs(i)); // <= cte
		else model.add(lhs[i] == LP->get_rhs(i)); // = cte
	}

	for (int i = 0; i < LP->get_p(); i++) {
		ptrCtes.add(IloRange(env, -IloInfinity, 0)); // change rhs (p)
		for (int j = 0; j < LP->get_n(); j++) {
			ptrCtes[i].setLinearCoef(x[j], LP->get_objective(i, j));
		}
		model.add(ptrCtes[i]);
	}

}

bool CplexModel::solveFeasibility(std::vector<double>& s) {

	try {
		// adjust constraints
		for (int i = 0; i < p; i++) {
			ptrCtes[i].setUB(s[i]); // not a valid epsilon
		}

		// solve
		bool solved;
		solved = cplex.solve();

		/*
		if (solved) {
			IloNumArray rlt(env);
			cplex.getValues(rlt, x);

			// print solution
			//for (int k = 0; k < n; k++) {
			//	std::cout << "x[" << k << "] = " << rlt[k] << std::endl;
			//}
		}
		else
		{
			std::cout << "Not feasible !" << std::endl;
		}*/
		return solved;
		
	}
	catch (IloException& e)
	{
		std::cout << e << std::endl;
		e.end();
	}

	//cplex.exportModel("check4.lp");
}

void CplexModel::retrieveSolutionFeasibility(std::vector<double>& s) {

	IloNumArray rlt(env);
	cplex.getValues(rlt, x);
	for (int k = 0; k < n; k++) {
		s[k] = rlt[k];
	}
}


// -------------- Weighted Sum Scalarizations --------------

void CplexModel::buildWeightedSumScalarization(MathematicalFormulation* LP) {

	// building LP caracteristics

	p = LP->get_p();
	n = LP->get_n();
	m = LP->get_m();

	// building variables

	for (int i = 0; i < LP->get_n(); i++) {
		x.add(IloNumVar(env, 0, 1, ILOFLOAT));
	}

	// building objective function

	IloExpr obj(env);
	for (int i = 0; i < LP->get_n(); i++) {
		obj += 1.0 * x[i]; // to be modified before solving
	}

	ptrObj = IloObjective(env, obj, IloObjective::Minimize);
	model.add(ptrObj);
	obj.end();

	// building constraints

	std::vector<IloExpr> lhs(0);
	for (int i = 0; i < LP->get_m(); i++) {
		lhs.push_back(IloExpr(env));
		for (int j = 0; j < LP->get_n(); j++) {
			lhs[i] += LP->get_constraint(i, j) * x[j];
		}

		if (LP->get_signCte(i) == 0) model.add(lhs[i] >= LP->get_rhs(i)); // >= cte
		else if (LP->get_signCte(i) == 1) model.add(lhs[i] <= LP->get_rhs(i)); // <= cte
		else model.add(lhs[i] == LP->get_rhs(i)); // = cte
	}
}

void CplexModel::solveWeightedSumScalarization(MathematicalFormulation& LP, std::vector<double> lambda) {

	try
	{
		// Modify objective coefficients
		double coef;
		for (int i = 0; i < n; i++) {
			coef = 0;
			for (int k = 0; k < p; k++) {
				coef += lambda[k] * LP.get_objective(k, i);
			}
			ptrObj.setLinearCoef(x[i], coef);
		}

		// solve
		cplex.solve();

		// Retrieve solution
		// IloNumArray rlt(env);
		// cplex.getValues(rlt, x);

		// print solution
		//for (int k = 0; k < n; k++) {
		//	std::cout << "x[" << k << "] = " << rlt[k] << std::endl;
		//}
	}
	catch (IloException& e)
	{
		std::cout << e << std::endl;
		e.end();
	}

	//cplex.exportModel("check5.lp");
}

double CplexModel::getWeightedSumValue(MathematicalFormulation& LP, std::vector<double> lambda) {

	solveWeightedSumScalarization(LP,lambda);

	return cplex.getObjValue();
}



 /* ==========================================================
		 Getters
  ========================================================= */

int LinearProgram::get_p() {
	return p;
}

int LinearProgram::get_n() {
	return n;
}

int LinearProgram::get_m() {
	return m;
}

int LinearProgram::get_objDir(int i) {
	return optDirection[i];
}

int MathematicalFormulation::get_objective(int i, int j) {
	return objMatrix[i][j];
}

int MathematicalFormulation::get_constraint(int i, int j) {
	return constrMatrix[i][j];
}

int MathematicalFormulation::get_rhs(int j) {
	return rhs[j].second;
}

int MathematicalFormulation::get_signCte(int j) {
	return rhs[j].first;
}

IloCplex* CplexModel::get_cplex() {
	return &cplex;
}

IloEnv* CplexModel::get_env() {
	return &env;
}

/* ==========================================================
		Setters
 ========================================================= */









// objective coefficient matrix

//getline(f, line);
//std::vector<std::vector<int>> objMatrix(m);

//for (int i = 0; i < m; i++) {
//	objMatrix[i]
//}