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

void MathematicalFormulation::formate_to_min() {

	for (int k = 0; k < p; k++) {
		if (optDirection[k] == 1) {
			for (int i = 0; i < n; i++) {
				objMatrix[k][i] *= -1;
			}
		}
	}

}

void CplexModel::build_dualBenson(MathematicalFormulation* LP) {

	// building LP caracteristics

	p = LP->get_p();
	n = LP->get_m() + LP->get_p();
	m = LP->get_n() + 1;
	optDirection.push_back(1);

	// building variables

	for (int i = 0; i < LP->get_m(); i++) {

		// take care of sign of constraints in primal for sign of variables in dual

		if (LP->get_signCte(i) == 0) x.add(IloNumVar(env, 0, IloInfinity, ILOFLOAT));
		else if(LP->get_signCte(i) == 1) x.add(IloNumVar(env, -IloInfinity, 0, ILOFLOAT));
		else x.add(IloNumVar(env, -IloInfinity, IloInfinity, ILOFLOAT));
		
	}
	for (int i = 0; i < LP->get_p(); i++) {
		x.add(IloNumVar(env, 0, IloInfinity, ILOFLOAT));
	}

	// building objective function

	
	IloExpr obj(env);
	for (int i = 0; i < LP->get_m(); i++) {
		obj += LP->get_rhs(i) * x[i];
	}
	for (int i = 0; i < LP->get_p(); i++) {
		obj -= 1.0 * x[LP->get_m() + i]; // replace 1 with value of y when using template
	}
	//IloObjective* objective;
	ptrObj = IloObjective(env, obj, IloObjective::Maximize); //IloMaximize IloObjective objective
	//ptrObj = &IloObjective(env, obj, IloObjective::Maximize);
	model.add(ptrObj);
	//ptrObj = objective;
	obj.end();

	// building constraints

	std::vector<IloExpr> lhs(0);
	for (int i = 0; i < LP->get_n(); i++) {
		lhs.push_back(IloExpr(env));
		for (int j = 0; j < LP->get_m(); j++) {
			lhs[i] += LP->get_constraint(j, i) * x[j];
		}
		for (int j = 0; j < LP->get_p(); j++) {
			lhs[i] -= (double) LP->get_objective(j, i) * x[LP->get_m() + j];
		}
		model.add(lhs[i] <= 0);

	}

	for (int i = 0; i < LP->get_p(); i++) {
		lhs.push_back(IloExpr(env));
		lhs[m - 1] += x[LP->get_m() + i];
	}
	model.add(lhs[m - 1] == 1);

	
	cplex.exportModel("check.lp");
}

void CplexModel::build_bestValidPoint(MathematicalFormulation* LP) {

	// building LP caracteristics

	p = LP->get_p();
	n = LP->get_n() + 1;
	m = LP->get_m() + LP->get_p();
	optDirection.push_back(1);

	// building variables

	for (int i = 0; i < LP->get_n(); i++) {
		x.add(IloNumVar(env, 0, IloInfinity, ILOFLOAT));
	}
	x.add(IloNumVar(env, 0, 1, ILOFLOAT));

	// builind objective function

	model.add(IloMaximize(env, x[n-1]));

	// building constraints

	std::vector<IloExpr> lhs(0);
	for (int i = 0; i < LP->get_m(); i++) {
		lhs.push_back(IloExpr(env));
		for (int j = 0; j < LP->get_n(); j++) {
			lhs[i] += LP->get_constraint(i, j) * x[j];
		}
		model.add(lhs[i] >= LP->get_rhs(i));
	}

	for (int i = 0; i < LP->get_p(); i++) {
		lhs.push_back(IloExpr(env));
		for (int j = 0; j < LP->get_n(); j++) {
			lhs[LP->get_m() + i] += LP->get_objective(i, j) * x[j];
		}
		lhs[LP->get_m() + i] += 1 * x[n - 1]; // change coef of last var (lambda)
		model.add(lhs[LP->get_m() + i] <= 0); // change rhs (p)
	}

	//cplex.exportModel("check.lp");
}

void CplexModel::solveDualBenson(std::vector<double>& y) {

	std::cout << "Final obj : " << cplex.getObjective() << std::endl;
	std::cout << "solve 1" << std::endl;
	cplex.solve();


	//model.getObject();
	std::cout << "Modifying the obj \n";
	for (int k = 0; k < p; k++) {
		std::cout << "Modifying variable " << n - p + k << "\n";
		ptrObj.setLinearCoef(x[n - p + k], -y[k]);
		//cplex.getObjective().setLinearCoef(x[n - p + k], -y[k]);
		//model.setLinearCoef(x[n - p + k], -y[k]);
		//std::cout << y[k] << " => " << cplex.getObjective() << std::endl;
	}
	std::cout << "Final obj : " << cplex.getObjective() << std::endl;

	//std::cout << "Obj : " << model.getObject() << std::endl;




	std::cout << "solve 2" << std::endl;

	cplex.solve();
	IloNumArray rlt(env);
	cplex.getValues(rlt, x);
	for (int k = 0; k < n; k++) {
		//rlt = cplex.getValue(k);
		std::cout << "x[" << k << "] = " << rlt[k] << std::endl;
	}
	std::cout << "z(x) = " << cplex.getObjValue() << std::endl;

	cplex.exportModel("check2.lp");
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

/* ==========================================================
		Setters
 ========================================================= */









// objective coefficient matrix

//getline(f, line);
//std::vector<std::vector<int>> objMatrix(m);

//for (int i = 0; i < m; i++) {
//	objMatrix[i]
//}