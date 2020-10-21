#pragma once
/**
* \author Nicolas Forget
*
* The class Model gather all kind of model to be solved. Two kind of models can be identified:
* 1) a mathematical model, stored as vectors and matrices of coefficients.
* 2) a cplex model, stored using ilocplex's interface. One subclass is created for each new model.
*/

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <ilcplex\ilocplex.h>

class Model {
protected:
	int p; //!< the number of objectives
	int n; //!< the number of variables
	int m; //!< the number of constraints
	std::vector<int> optDirection; //!< the optimization direction for each objective: 1 for maximization, 0 for minimization.

public:

	/*! \brief Default constructor of the model.
	 *
	 * This function creates an empty model.
	 */
	Model();

	/*! \brief Returns the number of objectives
	 *
	 * This function returns the right-hand side of the hyperplane.
	 * \return the number of objectives, as an int.
	 */
	int get_p();

	/*! \brief Returns the number of variables
	 *
	 * This function returns the number of variables.
	 * \return the number of variables, as an int.
	 */
	int get_n();

	/*! \brief Returns the number of constraints
	 *
	 * This function returns the number of constraints.
	 * \return the number of constraints, as an int.
	 */
	int get_m();

	/*! \brief Returns the optimization direction of a given objective
	 *
	 * This function returns the optimization direction of objective i.
	 * \param i integer. The index of the objective to look at.
	 * \return the optimization direction objective i, as an int.
	 */
	int get_objDir(int i);
};



// ===============================================================================================================================
//							MathematicalModel
// ===============================================================================================================================

class MathematicalModel : public Model {
private:
	std::vector<std::vector<int>> C; //!< Matrix of objective coefficients, as a vector of vector
	std::vector<std::vector<int>> A; //!< Matrix of constraint coefficients, as a vector of vector
	std::vector<int> constrSign; //!< Vector for the constraint signs.
	std::vector<int> b; //!< Vector for the right-hand side of the constraints.

public:
	
	/*! \brief Default constructor of a mathematical model.
	 *
	 * This function creates an empty MathematicalModel.
	 */
	MathematicalModel();

	/*! \brief Fill the MathematicalModel with the instance given in an instance file.
	 *
	 * This function fills the MathematicalModel object that corresponds to the instance found in the file.
	 * \param file string. The path and name of the instance file.
	 */
	void fill(std::string file);

	/*! \brief Print the objective matrix.
	 *
	 * This function prints the coefficients of the objective matrix in a user-friendly manner in the command line.
	 */
	void printObjective();

	/*! \brief Print the constraint matrix.
	 *
	 * This function prints the coefficients of the constraint matrix in a user-friendly manner in the command line.
	 */
	void printConstraints();

	/*! \brief Set all objective into minimization form.
	 *
	 * This function turns the objective in maximization to objective in minimization, i.e. instead of maximizing $z(x)$,
	 * it minimizes $-z(x)$. To do so, it just multiplies each coefficient of each objective k = 1,...,p by -optDirection[k]
	 */
	void formateToMin();

	/*! \brief Returns the value of the coefficient of variable $x_j$ in objective $i$.
	 *
	 * This function returns the value of the coefficient of variable $x_j$ in objective $i$ by reading the matrix C.
	 * \param i integer. The index of the objective.
	 * \param j integer. The index of the variable.
	 * \return the value of the coefficient, as an int.
	 */
	int get_objective(int i, int j);

	/*! \brief Returns the value of the coefficient of variable $x_j$ in constraint $i$.
	 *
	 * This function returns the value of the coefficient of variable $x_j$ in constraint $i$ by reading the matrix A.
	 * \param i integer. The index of the constraint.
	 * \param j integer. The index of the variable.
	 * \return the value of the coefficient, as an int.
	 */
	int get_constraint(int i, int j);

	/*! \brief Returns the value of the right-hand side of constraint $j$.
	 *
	 * This function returns the value of the right-hand side of constraint $j$ by reading the vector b.
	 * \param j integer. The index of the constraint.
	 * \return the value of the right-hand side, as an int.
	 */
	int get_rhs(int j);

	/*! \brief Returns the sign of constraint $j$.
	 *
	 * This function returns the sign of constraint $j$ by reading the vector constrSign.
	 * constrSign = 0 corresponds to >= constraint
	 * constrSign = 1 corresponds to <= constraint
	 * constrSign = 2 corresponds to = constraint
	 * \param j integer. The index of the constraint.
	 * \return the sign, as an int.
	 */
	int get_signCte(int j);
};



// ===============================================================================================================================
//							CplexModel
// ===============================================================================================================================

class CplexModel : public Model {
protected:
	IloEnv env; //!< cpelx object
	IloModel model; //!< cpelx object
	IloCplex cplex; //!< cpelx object

public:
	/*! \brief Default constructor of a cplex model.
	 *
	 * This function creates an empty CplexModel.
	 */
	CplexModel();

	/*! \brief Returns a pointer to the IloCplex object.
	 *
	 * This function returns a pointer to the IloCplex object.
	 * \return a pointer to the IloCplex object.
	 */
	IloCplex* get_cplex();

	/*! \brief Returns a pointer to the IloEnv object.
	 *
	 * This function returns a pointer to the IloEnv object.
	 * \return a pointer to the IloEnv object.
	 */
	IloEnv* get_env();
};



// ===============================================================================================================================
//							DualBensonModel
// ===============================================================================================================================

/**
 * This model is the dual model defined in Benson's method. It is used to determine the hyperplanes of the linear relaxation.
 * Given a point y on the boundary of $\mathcal{L} + \mathbb{R}^p$, it returns the equation of the facet y is located on.
 */

class DualBensonModel : public CplexModel {
private:
	IloNumVarArray u; //!< first set of variables, dual of the initial problem's constraints (A >= b)
	IloNumVarArray v; //!< second set of variables, dual of the upper bounds on the constraints
	IloNumVarArray w; //!< third set of variables, dual of the constraints on the objectives (Cx - et <= y)
	IloObjective ptrObj; //!< the objective function object, explicitely stored here to be modified prior to solving

public:
	/*! \brief Default constructor of Dual Benson Model.
	 *
	 * This function creates an empty model.
	 */
	DualBensonModel();

	/*! \brief Build the Dual Benson model.
	 *
	 * This function set up the variables, constraints and objective function of the linear program.
	 * \param LP MathematicalModel. This is the model of the initial problem, from which the dual is built.
	 */
	void build(MathematicalModel& LP);

	/*! \brief Solve the Dual Benson model.
	 *
	 * This function updates the objective coefficients that need to be changed and solve the model.
	 * \param y vector of doubles. It is the point located on the facet that is being computed.
	 */
	void solve(std::vector<double>& y);

	/*! \brief Extracts the normal vector of the facet computed.
	 *
	 * This function extract the normal vector of the facet computed with this model. It is given by the w vector.
	 * \return the normal vector, as a vector of double.
	 */
	std::vector<double> extractNormalVector();

	/*! \brief Extracts the constant of the equation of the facet.
	 *
	 * This function extract the constant of the equation of the facet computed. It given by $b^Tu - e^Tv$.
	 * \return the constant of the equation of the facet.
	 */
	double extractConstant(MathematicalModel& LP);
};



// ===============================================================================================================================
//							FurthestFeasiblePointModel
// ===============================================================================================================================

/**
 * Given two point $s,y \in \mathbb{R}^p$ such that $y \in \mathcal{L} + \mathbb{R}^p$ and $s \notin \mathcal{L} + \mathbb{R}^p$,
 * this model search for the feasible point of the edge defined by s and y that is the closest to s. In other words, it omputes
 * the feasible point of this edge the furthest away from y.
 */

class FurthestFeasiblePointModel : public CplexModel {
private:
	IloNumVarArray x; //!< variables from the initial problem
	IloNumVar lambda; //!< the value that we are looking for to compute the furthest feasible point
	IloRangeArray ptrCtes; //!< a pointer to the constraints on the objectives, to update them as necessary later

public:
	/*! \brief Default constructor of Furthest Feasible Point.
	 *
	 * This function creates an empty model.
	 */
	FurthestFeasiblePointModel();

	/*! \brief Build the Furthest Feasible Point model.
	 *
	 * This function set up the variables, constraints and objective function of the linear program.
	 * \param LP MathematicalModel. This is the model of the initial problem, from which this linear program is built.
	 */
	void build(MathematicalModel& LP);

	/*! \brief Solve the Furthest Feasible Point model.
	 *
	 * This function updates the constraints that need to be changed and solves the model.
	 * \param s vector of doubles. It is the point of the edge that is not located in $\mathcal{L} + \mathbb{R}^p$
	 * \param phat vector of doubles. It is the point of the edge that is located in $\mathcal{L} + \mathbb{R}^p$
	 */
	void solve(std::vector<double>& s, std::vector<double>& phat);

	/*! \brief Extracts the point searched from Cplex.
	 *
	 * This function extract the point we are looking for by solving this model. It is retrieved by extracting the value of
	 * lambda from the solved model, and by computing the convex combination of s and phat using lambda as the weight.
	 * Note that it calls the function solve.
	 * \param s vector of double. It is the point of the edge that is not located in $\mathcal{L} + \mathbb{R}^p$
	 * \param phat vector of doubles. It is the point of the edge that is located in $\mathcal{L} + \mathbb{R}^p$
	 * \return the objective vector of the computed point, as a vector of double.
	 */
	std::vector<double> extractPoint(std::vector<double>& s, std::vector<double>& phat);
};



// ===============================================================================================================================
//							FeasibilityCheckModel
// ===============================================================================================================================

/**
 * This model checks whether there exists a point feasible for the initial problem that weakly dominates a given point $s$ of
 * the objective space.
 * Note: currently there is no objective function. It may be changed to min a sum of slack variables, to make sure it reaches
 * an already existing extreme point, that may speed-up the warmstarting procedure in the branch-and-bound algorithm.
 */

class FeasibilityCheckModel : public CplexModel {
private:
	IloNumVarArray x; //!< variables from the initial problem
	IloRangeArray ptrCtes; //!< a pointer to the constraints on the objectives, to update them as necessary later

public:
	/*! \brief Default constructor of Feasibility Check.
	 *
	 * This function creates an empty model.
	 */
	FeasibilityCheckModel();

	/*! \brief Build the Feasibility Check model.
	 *
	 * This function set up the variables, constraints and objective function of this linear program.
	 * \param LP MathematicalModel. This is the model of the initial problem, from which this linear program is built.
	 */
	void build(MathematicalModel& LP);

	/*! \brief Solve the Feasibility Check model.
	 *
	 * This function updates the constraints that need to be changed and solves the model.
	 * \param s vector of doubles. The used reference point of the objective space.
	 * \return true if it is feasible and solved to optimality, false otherwise
	 */
	bool solve(std::vector<double>& s);

	/*! \brief retreive the pre-image computed.
	 *
	 * This function retreive from cplex the pre-image of the reference point used to solve the model by reading the x vector
	 * and storing it in s.
	 * \param s vector of doubles. The vector in which we store the computed pre-image of the reference point.
	 */
	void retrieveSolutionFeasibility(std::vector<double>& s);
};



// ===============================================================================================================================
//							WeightedSumScalarizationModel
// ===============================================================================================================================

/**
 * This model computes the weighted sum scalarization of the initial problem.
 */

class WeightedSumModel : public CplexModel {
private:
	IloNumVarArray x; //!< variables from the initial problem
	IloObjective ptrObj; //!< the objective function object, explicitely stored here to be modified prior to solving

public:
	/*! \brief Default constructor of Weighted Sum.
	 *
	 * This function creates an empty model.
	 */
	WeightedSumModel();

	/*! \brief Build the Feasibility Check model.
	 *
	 * This function set up the variables, constraints and objective function of this linear program.
	 * \param LP MathematicalModel. This is the model of the initial problem, from which this linear program is built.
	 */
	void build(MathematicalModel& LP);

	/*! \brief Solve the Weighted Sum model.
	 *
	 * This function updates the objective coefficient by making a weighted sum of the objective functions, given the weight
	 * vector lambda.
	 * \param LP MathematicalModel. This is the model of the initial problem, from which this linear program is built.
	 * \param lambda vector of double. This is the weight vector.
	 * \return true if it is feasible and solved to optimality, false otherwise [ToDo]
	 */
	void solve(MathematicalModel& LP, std::vector<double> lambda); // [toDo] true is solved to optimality, false otherwise

	/*! \brief retrieve the objective value of the weighted sum.
	 *
	 * This function retreive from cplex the objective value of the weighted sum. It in particular gives the right-hand side
	 * of the equation of the hyperplane defined by the normal vector lambda and that contains at least one feasible point.
	 * \param LP MathematicalModel. This is the model of the initial problem, from which this linear program is built.
	 * \param lambda vector of double. This is the weight vector.
	 * \return the optimal objective value, as a double.
	 */
	double retrieveObjectiveValue(MathematicalModel& LP, std::vector<double> lambda);
};
