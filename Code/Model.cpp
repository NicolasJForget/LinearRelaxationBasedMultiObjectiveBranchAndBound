#include "Model.h"

/* ==========================================================
		Constructors
 ========================================================= */

/*! \brief Default constructor of the model.
 */
Model::Model() : p(0), n(0), m(0), optDirection(0) {}

/* ==========================================================
		 Getters
  ========================================================= */

/*! \brief Returns the number of objectives
 *
 * \return the number of objectives, as an int.
 */
int Model::get_p() {
	return p;
}

/*! \brief Returns the number of variables
 *
 * \return the number of variables, as an int.
 */
int Model::get_n() {
	return n;
}

/*! \brief Returns the number of constraints
 *
 * \return the number of constraints, as an int.
 */
int Model::get_m() {
	return m;
}

/*! \brief Returns the optimization direction of a given objective
 *
 * \param i integer. The index of the objective to look at.
 * \return the optimization direction objective i, as an int.
 */
int Model::get_objDir(int i) {
	return optDirection[i];
}




// ===============================================================================================================================
//							MathematicalModel
// ===============================================================================================================================

/* ==========================================================
		Constructors
 ========================================================= */

 /*! \brief Default constructor of a mathematical model.
  */
MathematicalModel::MathematicalModel() : Model(), C(0), A(0), constrSign(0), b(0) {};

/*! \brief Constructor of a mathematical model for a given instance.
 *
 * \param file string. The path and name of the instance file.
 */
MathematicalModel::MathematicalModel(std::string file) : MathematicalModel() {
	fill(file);
	formateToMin();
}

/* ==========================================================
		Regular Methods
 ========================================================= */

/*! \brief Fill the MathematicalModel with the instance given in an instance file.
 *
 * \param file string. The path and name of the instance file.
 */
void MathematicalModel::fill(std::string file) {

	std::ifstream f(file);

	if (f) {

		std::string str;
		std::string line;
		std::string coef;

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

		optDirection.resize(p);
		for (int i = 0; i < p; i++) {
			iss >> str;
			if (str == "maxsum") {
				optDirection[i] = 1;
			}
			else if (str == "minsum") {
				optDirection[i] = -1;
			}
		}

		// objective coefficients matrix

		getline(f, line);
		for (int k = 0; k < p; k++) {
			C.push_back(std::vector<int>(n));
			getline(f, line);
			iss.clear();
			iss.str(line);
			int i = 0;
			while (i < n) {
				iss >> coef;
				C[k][i] = stoi(coef);
				++i;
			}
		}

		// constraint coefficients matrix

		getline(f, line);
		for (int j = 0; j < m; j++) {
			A.push_back(std::vector<int>(n));
			getline(f, line);
			iss.clear();
			iss.str(line);
			int i = 0;
			while (i < n) {
				iss >> coef;
				A[j][i] = stoi(coef);
				++i;
			}
		}

		// right-hand side of the constraints

		getline(f, line);
		constrSign.resize(m);
		b.resize(m);
		for (int j = 0; j < m; j++) {
			getline(f, line);
			iss.clear();
			iss.str(line);
			iss >> coef;
			constrSign[j] = stoi(coef);
			iss >> coef;
			b[j] = stoi(coef);
		}

		// close file
		f.close();
	}
	else {
		std::cout << "Problem with file opening" << std::endl;
	}
}

/*! \brief Print the objective matrix.
*/
void MathematicalModel::printObjective() {

	std::cout << "\n\n";
	for (int k = 0; k < p; k++) {
		std::cout << "objective " << k + 1 << " : ";
		for (int i = 0; i < n; i++) {
			std::cout << C[k][i] << " ";
		}
		std::cout << "\n";
	}
}

/*! \brief Print the constraint matrix.
 */
void MathematicalModel::printConstraints() {

	std::cout << "\n\n";
	for (int k = 0; k < m; k++) {
		std::cout << "constraint " << k + 1 << " : ";
		for (int i = 0; i < n; i++) {
			std::cout << A[k][i] << " ";
		}
		if (constrSign[k] == 1) {
			std::cout << "  <=   ";
		}
		else if (constrSign[k] == 2) {
			std::cout << "  =   ";
		}
		else {
			std::cout << "  >=   ";
		}
		std::cout << b[k] << "\n";
	}
}

/*! \brief Set all objective into minimization form.
 */
void MathematicalModel::formateToMin() {
	for (int k = 0; k < p; k++) {
		if (optDirection[k] == 1) {
			for (int i = 0; i < n; i++) {
				C[k][i] *= -1;
			}
		}
	}
}

/* ==========================================================
		 Getters
  ========================================================= */

/*! \brief Returns the value of the coefficient of variable $x_j$ in objective $i$.
 *
 * \param i integer. The index of the objective.
 * \param j integer. The index of the variable.
 * \return the value of the coefficient, as an int.
 */
int MathematicalModel::get_objective(int i, int j) {
	return C[i][j];
}

/*! \brief Returns the value of the coefficient of variable $x_j$ in constraint $i$.
 *
 * \param i integer. The index of the constraint.
 * \param j integer. The index of the variable.
 * \return the value of the coefficient, as an int.
 */
int MathematicalModel::get_constraint(int i, int j) {
	return A[i][j];
}

/*! \brief Returns the value of the right-hand side of constraint $j$.
 *
 * \param j integer. The index of the constraint.
 * \return the value of the right-hand side, as an int.
 */
int MathematicalModel::get_rhs(int j) {
	return b[j];
}

/*! \brief Returns the sign of constraint $j$.
 *
 * \param j integer. The index of the constraint.
 * \return the sign, as an int.
 */
int MathematicalModel::get_signCte(int j) {
	return constrSign[j];
}



// ===============================================================================================================================
//							CplexModel
// ===============================================================================================================================

/* ==========================================================
		Constructors
 ========================================================= */

/*! \brief Default constructor of a cplex model.
 */
CplexModel::CplexModel() : Model() {
	model = IloModel(env);
	cplex = IloCplex(model);
	cplex.setOut(env.getNullStream());
}

/* ==========================================================
		 Getters
  ========================================================= */

/*! \brief Returns a pointer to the IloCplex object.
 *
 * \return a pointer to the IloCplex object.
 */
IloCplex* CplexModel::get_cplex() {
	return &cplex;
}

/*! \brief Returns a pointer to the IloEnv object.
 *
 * \return a pointer to the IloEnv object.
 */
IloEnv* CplexModel::get_env() {
	return &env;
}



// ===============================================================================================================================
//							DualBensonModel
// ===============================================================================================================================

/* ==========================================================
		Constructors
 ========================================================= */

/*! \brief Default constructor of a DualBenson's model.
 */
DualBensonModel::DualBensonModel() : CplexModel() {
	u = IloNumVarArray(env);
	vUB = IloNumVarArray(env);
	vLB = IloNumVarArray(env);
	w = IloNumVarArray(env);
	ptrObj = IloObjective(env);
}

/* ==========================================================
		Regular Methods
 ========================================================= */

/*! \brief Build the Dual Benson model.
 *
 * \param LP MathematicalModel. This is the model of the initial problem, from which the dual is built.
 */
void DualBensonModel::build(MathematicalModel& LP) {

	try {
		// building LP caracteristics

		p = LP.get_p();
		n = LP.get_m() + LP.get_n() + LP.get_p(); //LP->get_m() + LP->get_p()
		m = LP.get_n() + 1;
		optDirection.push_back(1);

		// building variables

		for (int i = 0; i < LP.get_m(); i++) { // u variables

			// take care of sign of constraints in primal for sign of variables in dual

			if (LP.get_signCte(i) == 0) u.add(IloNumVar(env, 0, IloInfinity, ILOFLOAT)); // >= cte
			else if (LP.get_signCte(i) == 1) u.add(IloNumVar(env, -IloInfinity, 0, ILOFLOAT)); // <= cte
			else u.add(IloNumVar(env, -IloInfinity, IloInfinity, ILOFLOAT)); // = cte

		}
		for (int i = 0; i < LP.get_n(); i++) { // v variables (for UB ctes on variables)
			vUB.add(IloNumVar(env, 0, IloInfinity, ILOFLOAT));
		}
		for (int i = 0; i < LP.get_n(); i++) {
			vLB.add(IloNumVar(env, 0, IloInfinity, ILOFLOAT));
		}
		for (int i = 0; i < LP.get_p(); i++) { // w variables 
			w.add(IloNumVar(env, 0, IloInfinity, ILOFLOAT));
		}

		// building objective function

		IloExpr obj(env);
		for (int i = 0; i < LP.get_m(); i++) { // u variables
			obj += LP.get_rhs(i) * u[i];
		}
		for (int i = 0; i < LP.get_n(); i++) { // vUB variables (UB on variables)
			obj -= 1.0 * vUB[i];
		}
		for (int i = 0; i < LP.get_n(); i++) { // vLB variables (LB on variables)
			obj += 0.0 * vLB[i];
		}
		for (int i = 0; i < LP.get_p(); i++) { // w variables
			obj -= 1.0 * w[i]; // replace 1 with value of y when using template
		}

		ptrObj = IloObjective(env, obj, IloObjective::Maximize);
		model.add(ptrObj);
		obj.end();

		// building core constraints

		std::vector<IloExpr> lhs(0);
		for (int i = 0; i < LP.get_n(); i++) {
			lhs.push_back(IloExpr(env));
			for (int j = 0; j < LP.get_m(); j++) { // A^T.u
				lhs[i] += LP.get_constraint(j, i) * u[j];
			}
			lhs[i] -= 1.0 * vUB[i]; // - vUB
			lhs[i] += 1.0 * vLB[i]; // + vLB
			for (int j = 0; j < LP.get_p(); j++) { // - C^T.w
				lhs[i] -= (double)LP.get_objective(j, i) * w[j];
			}
			model.add(lhs[i] <= 0);
		}

		// constraint on the normal vector

		lhs.push_back(IloExpr(env));
		for (int i = 0; i < LP.get_p(); i++) {
			lhs[m - 1] += w[i];
		}
		model.add(lhs[m - 1] == 1);
	}
	catch (IloException& e)
	{
		std::cout << e << std::endl;
		e.end();
	}
}

/*! \brief Solve the Dual Benson model.
 *
 * \param y vector of doubles. It is the point located on the facet that is being computed.
 */
void DualBensonModel::solve(std::vector<double>& y) {

	try
	{
		// Modify objective coefficients

		for (int k = 0; k < p; k++) {
			ptrObj.setLinearCoef(w[k], -y[k]);
		}

		// solve

		cplex.solve();
	}
	catch (IloException& e)
	{
		std::cout << e << std::endl;
		e.end();
	}
}

/*! \brief Extracts the normal vector of the facet computed.
 *
 * \return the normal vector, as a vector of double.
 */
std::vector<double> DualBensonModel::extractNormalVector() {

	std::vector<double> theVector(p);
	//std::cout << "status: " << cplex.getStatus() << std::endl;

	//cplex.exportModel("pksabugsvp.lp");
	IloNumArray rlt(env);
	cplex.getValues(rlt, w);

	for (int k = 0; k < p; k++) {
		theVector[k] = rlt[k];
	}

	return theVector;
}

/*! \brief Extracts the constant of the equation of the facet.
 *
 * \return the constant of the equation of the facet.
 */
double DualBensonModel::extractConstant(MathematicalModel& LP, BranchingDecisions* branchDec) {

	double theRightHandSide = 0;

	IloNumArray rlt(env);
	cplex.getValues(rlt, u);
	IloNumArray rlt2(env);
	cplex.getValues(rlt2, vUB);
	IloNumArray rlt3(env);
	cplex.getValues(rlt3, vLB);

	//IloExpr coef = ptrObj.getExpr();
	//IloExpr::LinearIterator it = coef.getLinearIterator();

	for (int i = 0; i < LP.get_m(); i++) {
		theRightHandSide += LP.get_rhs(i) * rlt[i];
		//theRightHandSide += it.getCoef() * rlt[i];
		//it.operator++();
	}
	for (int i = 0; i < LP.get_n(); i++) {
		theRightHandSide -= branchDec->ub[i] * rlt2[i];
		//theRightHandSide += it.getCoef() * rlt2[i]; //LP.get_m() + 
		//it.operator++();
	}
	for (int i = 0; i < LP.get_n(); i++) {
		theRightHandSide += branchDec->lb[i] *rlt3[i]; //LP.get_m() + 
		//it.operator++();
	}

	return theRightHandSide;
}

/*! \brief Adjust the bounds of the variables & constraints given some branching decisions
 *
 * \param bd BranchingDecisions. The data structure that describes the branching decisions to apply.
 */
void DualBensonModel::adjustBounds(BranchingDecisions& bd) {
	for (int i = 0; i < m - 1; i++) { // m - 1 = n of initial pb
		ptrObj.setLinearCoef(vUB[i], -bd.ub[i]);
		ptrObj.setLinearCoef(vLB[i], bd.lb[i]);
	}
}



// ===============================================================================================================================
//							FurthestFeasiblePointModel
// ===============================================================================================================================

/* ==========================================================
		Constructors
 ========================================================= */

 /*! \brief Default constructor of Furthest Feasible Point.
	  *
	  * This function creates an empty model.
	  */
FurthestFeasiblePointModel::FurthestFeasiblePointModel() : CplexModel() {

	x = IloNumVarArray(env);
	lambda = IloNumVar(env);
	ptrCtes = IloRangeArray(env);
}

/* ==========================================================
		Regular Methods
 ========================================================= */

/*! \brief Build the Furthest Feasible Point model.
 *
 * \param LP MathematicalModel. This is the model of the initial problem, from which this linear program is built.
 */
void FurthestFeasiblePointModel::build(MathematicalModel& LP) {

	try {
		// building LP caracteristics

		p = LP.get_p();
		n = LP.get_n() + 1;
		m = LP.get_m() + LP.get_p();
		optDirection.push_back(1);

		// building variables

		for (int i = 0; i < LP.get_n(); i++) {
			x.add(IloNumVar(env, 0, 1, ILOFLOAT));
		}
		lambda = IloNumVar(env, 0, 1, ILOFLOAT);

		// builind objective function

		model.add(IloMaximize(env, lambda));

		// building the constraints extracted from the initial problem (stored in LP)

		std::vector<IloExpr> lhs(0);
		for (int i = 0; i < LP.get_m(); i++) { // for each constraint
			lhs.push_back(IloExpr(env));
			for (int j = 0; j < LP.get_n(); j++) {
				lhs[i] += LP.get_constraint(i, j) * x[j];
			}
			if (LP.get_signCte(i) == 0) model.add(lhs[i] >= LP.get_rhs(i)); // >= cte
			else if (LP.get_signCte(i) == 1) model.add(lhs[i] <= LP.get_rhs(i)); // <= cte
			else model.add(lhs[i] == LP.get_rhs(i)); // = cte
		}

		// building the constraints on the objectives

		for (int i = 0; i < LP.get_p(); i++) { // for each objective
			ptrCtes.add(IloRange(env, -IloInfinity, 0)); // change rhs (phat)
			for (int j = 0; j < LP.get_n(); j++) {
				ptrCtes[i].setLinearCoef(x[j], LP.get_objective(i, j));
			}
			ptrCtes[i].setLinearCoef(lambda, 1); // change coef of lambda
			model.add(ptrCtes[i]);
		}
	}
	catch (IloException& e)
	{
		std::cout << e << std::endl;
		e.end();
	}
}

/*! \brief Solve the Dual Benson model.
 *
 * \param s vector of doubles. It is the point of the edge that is not located in $\mathcal{L} + \mathbb{R}^p$
 * \param phat vector of doubles. It is the point of the edge that is located in $\mathcal{L} + \mathbb{R}^p$
 */
void FurthestFeasiblePointModel::solve(std::vector<double>& s, std::vector<double>& phat) {

	try {
		// adjust constraints
		for (int i = 0; i < p; i++) {
			ptrCtes[i].setUB(phat[i]);
			ptrCtes[i].setLinearCoef(lambda, phat[i] - s[i]);
		}

		// solve
		cplex.solve();

	}
	catch (IloException& e)
	{
		std::cout << e << std::endl;
		e.end();
	}
}

/*! \brief Extracts the point searched from Cplex.
 *
 * \param s vector of double. It is the point of the edge that is not located in $\mathcal{L} + \mathbb{R}^p$
 * \param phat vector of doubles. It is the point of the edge that is located in $\mathcal{L} + \mathbb{R}^p$
 * \return the objective vector of the computed point, as a vector of double.
 */
std::vector<double> FurthestFeasiblePointModel::extractPoint(std::vector<double>& s, std::vector<double>& phat) {

	solve(s, phat);

	double lambda = cplex.getObjValue();
	std::vector<double> thePoint(p);
	for (int k = 0; k < p; k++) {
		thePoint[k] = lambda * s[k] + (1 - lambda) * phat[k];
	}

	return thePoint;
}

/*! \brief Adjust the bounds of the variables & constraints given some branching decisions
 *
 * \param bd BranchingDecisions. The data structure that describes the branching decisions to apply.
 */
void FurthestFeasiblePointModel::adjustBounds(BranchingDecisions& bd) {
	for (int i = 0; i < n - 1; i++) {
		x[i].setBounds(bd.lb[i], bd.ub[i]);
	}
}



// ===============================================================================================================================
//							FeasibleCheckModel
// ===============================================================================================================================

/* ==========================================================
		Constructors
 ========================================================= */

 /*! \brief Default constructor of Feasibility Check.
  *
  * This function creates an empty model.
  */
FeasibilityCheckModel::FeasibilityCheckModel() : CplexModel() {
	x = IloNumVarArray(env);
	ptrCtes = IloRangeArray(env);
}

/* ==========================================================
		Regular Methods
 ========================================================= */

 /*! \brief Build the Feasibility Check model.
  *
  * \param LP MathematicalModel. This is the model of the initial problem, from which this linear program is built.
  */
void FeasibilityCheckModel::build(MathematicalModel& LP) {

	try {

		// building LP caracteristics

		p = LP.get_p();
		n = LP.get_n();
		m = LP.get_m() + LP.get_p();

		// building variables
	
		for (int i = 0; i < LP.get_n(); i++) {
			x.add(IloNumVar(env, 0, 1, ILOFLOAT));
		}
	
		// building objective function
	
		model.add(IloMaximize(env, 0));
	
		// building constraints from the original problem
	
		std::vector<IloExpr> lhs(0);
		for (int i = 0; i < LP.get_m(); i++) {
			lhs.push_back(IloExpr(env));
			for (int j = 0; j < LP.get_n(); j++) {
				lhs[i] += LP.get_constraint(i, j) * x[j];
			}

			if (LP.get_signCte(i) == 0) model.add(lhs[i] >= LP.get_rhs(i)); // >= cte
			else if (LP.get_signCte(i) == 1) model.add(lhs[i] <= LP.get_rhs(i)); // <= cte
			else model.add(lhs[i] == LP.get_rhs(i)); // = cte
		}

		// building constraints on the objective functions
	
		for (int i = 0; i < LP.get_p(); i++) {
			ptrCtes.add(IloRange(env, -IloInfinity, 0)); // change rhs (p)
			for (int j = 0; j < LP.get_n(); j++) {
				ptrCtes[i].setLinearCoef(x[j], LP.get_objective(i, j));
			}
			model.add(ptrCtes[i]);
		}
	}
	catch (IloException& e)
	{
		std::cout << e << std::endl;
		e.end();
	}
}

/*! \brief Solve the Feasibility Check model.
 *
 * \param s vector of doubles. The used reference point of the objective space.
 * \return true if it is feasible and solved to optimality, false otherwise
 */
bool FeasibilityCheckModel::solve(std::vector<double>& s) {

	try {
		// update constraints
		for (int i = 0; i < p; i++) {
			ptrCtes[i].setUB(s[i]);
		}

		// solve
		bool solved;
		solved = cplex.solve();

		return solved;

	}
	catch (IloException& e)
	{
		std::cout << e << std::endl;
		e.end();
	}
}

/*! \brief retreive the pre-image computed.
 *
 * \param s vector of doubles. The vector in which we store the computed pre-image of the reference point.
 */
void FeasibilityCheckModel::retrieveSolutionFeasibility(std::vector<double>& s) {

	IloNumArray rlt(env);
	cplex.getValues(rlt, x);
	for (int k = 0; k < n; k++) {
		s[k] = rlt[k];
	}
}

/*! \brief Adjust the bounds of the variables & constraints given some branching decisions
 *
 * \param bd BranchingDecisions. The data structure that describes the branching decisions to apply.
 */
void FeasibilityCheckModel::adjustBounds(BranchingDecisions& bd) {
	for (int i = 0; i < n; i++) {
		x[i].setBounds(bd.lb[i], bd.ub[i]);
	}
}



// ===============================================================================================================================
//							WeightedSumModel
// ===============================================================================================================================

/* ==========================================================
		Constructors
 ========================================================= */

 /*! \brief Default constructor of Weighted Sum.
  */
WeightedSumModel::WeightedSumModel() {
	x = IloNumVarArray(env);
	ptrObj = IloObjective(env);
}

/* ==========================================================
		Regular Methods
 ========================================================= */

/*! \brief Build the Feasibility Check model.
 *
 * \param LP MathematicalModel. This is the model of the initial problem, from which this linear program is built.
 */
void WeightedSumModel::build(MathematicalModel& LP) {

	// building LP caracteristics

	p = LP.get_p();
	n = LP.get_n();
	m = LP.get_m();

	// building variables

	for (int i = 0; i < LP.get_n(); i++) {
		x.add(IloNumVar(env, 0, 1, ILOFLOAT));
	}

	// building objective function

	IloExpr obj(env);
	for (int i = 0; i < LP.get_n(); i++) {
		obj += 1.0 * x[i]; // to be modified before solving
	}

	ptrObj = IloObjective(env, obj, IloObjective::Minimize);
	model.add(ptrObj);
	obj.end();

	// building constraints

	std::vector<IloExpr> lhs(0);
	for (int i = 0; i < LP.get_m(); i++) {
		lhs.push_back(IloExpr(env));
		for (int j = 0; j < LP.get_n(); j++) {
			lhs[i] += LP.get_constraint(i, j) * x[j];
		}

		if (LP.get_signCte(i) == 0) model.add(lhs[i] >= LP.get_rhs(i)); // >= cte
		else if (LP.get_signCte(i) == 1) model.add(lhs[i] <= LP.get_rhs(i)); // <= cte
		else model.add(lhs[i] == LP.get_rhs(i)); // = cte
	}
}

/*! \brief Solve the Weighted Sum model.
 *
 * \param LP MathematicalModel. This is the model of the initial problem, from which this linear program is built.
 * \param lambda vector of double. This is the weight vector.
 * \return true if it is feasible and solved to optimality, false otherwise [ToDo]
 */
void WeightedSumModel::solve(MathematicalModel& LP, std::vector<double> lambda) {

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
	}
	catch (IloException& e)
	{
		std::cout << e << std::endl;
		e.end();
	}
}

/*! \brief retrieve the objective value of the weighted sum.
 *
 * \param LP MathematicalModel. This is the model of the initial problem, from which this linear program is built.
 * \param lambda vector of double. This is the weight vector.
 * \return the optimal objective value, as a double.
 */
double WeightedSumModel::retrieveObjectiveValue(MathematicalModel& LP, std::vector<double> lambda) {
	solve(LP, lambda);
	return cplex.getObjValue();
}

/*! \brief Adjust the bounds of the variables & constraints given some branching decisions
 *
 * \param bd BranchingDecisions. The data structure that describes the branching decisions to apply.
 */
void WeightedSumModel::adjustBounds(BranchingDecisions& bd) {
	for (int i = 0; i < n; i++) {
		x[i].setBounds(bd.lb[i], bd.ub[i]);
	}
}