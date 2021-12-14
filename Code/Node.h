/**
* \author Nicolas Forget
*
* This class describe a node of the branch-and-bound tree.
*/

#pragma once
#include <stdlib.h>
#include <list>
#include <vector>
#include "Model.h"
#include "LB2.h"
#include "UB.h"
#include "GlobalConstants.h"
#include "Stat.h"
#include "timer.hpp"
#include "SLUB.h"
#include "Lub.h"
#include "solution.h"
#include "Tree.h"

class TreeManager;

class Node
{
private:
	MathematicalModel* P; //!< a pointer to the mathematical formulation to the initial problem.
	ProbingModel* prob; //!< a pointer to the probing (cplex) model.
	VariableFixingModel* varFix; //!< a pointer to the (cplex) model used to fix variables when probing.
	WeightedSumModel* ws; //!< a pointer to the (cplex) model that is used to compute weighted sums,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
	Parameters* param; //!< a pointer to the parameters of the algorithm
	LowerBoundSet* LB; //!< Lower bound set of this node. As a pointer for inheritance and virtual functions to work properly.
	UpperBoundSet* UB; //!< A pointer to the global upper bound set.
	double score; //!< a score, used for node selection
	BranchingDecisions branchingDecision; //!< a data structure that save the branching decisions made in this node. See GlobalConstants.h for more details.
	int splittingIndex; //!< the index of the last variable splitted, i.e. the branching decision that led to the creation of this node.
	std::list<int> ndLub; //!< list of non-dominated local upper bounds.
	std::list<LocalUpperBound*> lubDomi; //!< ptrs to the dominated lubs of the node
	Statistics* stat; //!< a pointer to the data structure managing statistics
	//StatProbing* sp;
	int depth;
	Timer timeLBComputation; // depreciated
	Timer timeDominanceTest; // depreciated
	Timer timeUpdateUB; // 
	std::vector<std::vector<int>> coverCuts; // description of cover cuts. First column is rhs, reminder are indices of variables in the cut. One row per cut.
	//int nbCuts; // number of cuts generated at this node
	int iteration; //!< the current node is the iteration-th node explored

	// used for variable selection
	double sLbFatherNode; //!< the score value of the lower bound set from the father node
	double sLbCurrentNode; //!< the score value of the lower bound set of the current node
	double sXiFatherNode; //!< the score of the branching variable in the father node
	std::vector<double> sXiCurrentNode; //!< the score of each variable in the current node

public:
	/*! \brief Default constructor of a node.
	 *
	 * This function builds an empty node.
	 */
	Node();

	/* \brief Destructor of the node.
	 *
	 * This function destroy the node. It in particular desallocate the lower bound set, which is specifically defined for
	 * this node, unlike the other pointer-type members.
	 */
	~Node();

	/*! \brief Default constructor of the root node.
	 *
	 * This function builds the root node.
	 * \param lp MathematicalModel*. A pointer to the initial problem.
	 * \param par Parameters*. The parameters of the algorithm.
	 * \param stat Statistics*. A pointer to the statistics of the BB
	 */
	Node(MathematicalModel* lp, Parameters* par, Statistics* stat, UpperBoundSet* U);

	/*! \brief Creates a node given a splitting index and a bound.
	 *
	 * This function creates a new node that is a sub-problem of node nd, splitted on variable x[i].
	 * \param nd Node*. A pointer to the node it is created from.
	 * \param index int. The index of the variable to split.
	 * \param bound int. States the sign of the new constraint: are we adding an upper or a lower bound value on a variable?
	 * \param val int. The value of the bound added on the variable, i.e. the right-hand side of the new branching constraint.
	 * \param slub SLUB. The super local upper bound that defines the sub-problem in the objective space.
	 */
	Node(Node& nd, int index, int bound, int val, SLUB& slub);

	/*! \brief Creates a node given a splitting index and a bound.
	 *
	 * This function creates a new node that is a sub-problem of node nd, splitted on variable x[i].
	 * \param nd Node*. A pointer to the node it is created from.
	 * \param newbd. New branching decisions computed beforehand. Used when probing.
	 * \param index int. The index of the variable to split.
	 * \param bound int. States the sign of the new constraint: are we adding an upper or a lower bound value on a variable?
	 * \param val int. The value of the bound added on the variable, i.e. the right-hand side of the new branching constraint.
	 * \param slub SLUB. The super local upper bound that defines the sub-problem in the objective space.
	 */
	Node(Node& nd, BranchingDecisions* newbd, int index, SLUB& slub);

	/*! \brief Process the node
	 *
	 * This function process the node by computing a lower bound set, updating the lower bound set, and detecting 
	 * if the node can be fathomed.
	 * \param U UpperBoundSet. The upper bound set of the B&B.
	 * \param iteration int. The number of the iteration in the B&B.
	 */
	void process(UpperBoundSet& U, int iteration);

	/*! \brief Check whether the node is fathomed.
	 *
	 * This function checks whether this node is fathomed.
	 * \return true is the node is fathomed, false otherwise.
	 */
	bool isFathomed();

	/*! \brief Split the problem in the objective space
	 *
	 * This function splits the problem contained in this node in the objective space. For each sub-problem created, it calls
	 * the split function in the variable space. Store the new nodes created in the list Q.
	 * \param Q pointer to a list of pointers of Nodes. The list in which the new nodes created are stored.
	 * \param U pointer to the upper bound set. Used for computing objective branching.
	 */
	//void splitOS(std::list<Node*>* Q, UpperBoundSet* U, int iteration);
	void splitOS(TreeManager* T, UpperBoundSet* U, int iteration);

	/*! \brief Split the problem in the variable space
	 *
	 * This function splits the problem contained in this node in the variable space. Store the new nodes created in the list Q.
	 * \param Q pointer to a list of pointers of Nodes. The list in which the new nodes created are stored.
	 * \param slub SLUB. Defines the sub-problem in the objective space split in the decision space.
	 */
	//void splitVS(std::list<Node*>* Q, SLUB& slub, int iteration, UpperBoundSet* U);
	void splitVS(TreeManager* T, SLUB& slub, int iteration, UpperBoundSet* U);

	/*! \brief Prints the node
	 *
	 * This function prints the nodem by showing its branching decisions, and those of its father nodes.
	 */
	void print();

	/*! \brief Prints the status of the node.
	 *
	 * This function prints the status of the node, by reading the status of its lower bound set.
	 */
	void showStatus();

	/*! \brief Write statistics of this node
	 *
	 * This function writes the statistics of this node.
	 * \param stat Statistics*. Pointer to the struct where statistics are recorded.
	 */
	//void getStatistics(Statistics* stat);

	/*! \brief Rreturn the status of the node
	 *
	 * Returns the status of the node (infeasible, optimal, dominated...). See GlobalConstants.h for status.
	 * \return status of the node, as an integer
	 */
	int getStatus();

	/*! \brief Rreturn the split index.
	 *
	 * Returns the last split index in this node.
	 * \return the index, as an integer
	 */
	int getSplittingIndex();

	/*! \brief Rreturn the depth of the node.
	 *
	 * Returns the depth fo the node in the BB tree
	 * \return the depth, as an integer
	 */
	int getDepth();

	/*! \brief Print the LB of the node.
	 */
	void showLB();

	bool isOurCulprit();

	/*! \brief Try to fix integer bounds to each variable by solving LPs and looking at their feasibility
	 *
	 * \param BranchingDecisions newBd. Describes a potential futur subproblem.
	 * \return true if this node is detected as infeasible after checking the LPs.
	 */
	bool fixBounds(BranchingDecisions* newBd0, BranchingDecisions* newBd1, UpperBoundSet* U, int iteration);

	bool checkSubProblem(BranchingDecisions* newBd, int i, int v);

	/* Branch on the variable that has an ideal point with a dimension the closest to the slub
	*/
	int selectBranchingIndex1(BranchingDecisions* newBd0, BranchingDecisions* newBd1, std::vector<BranchingDecisions>& dec0, std::vector<BranchingDecisions>& dec1);

	/* Branch on the variable with the smallest volume of the cube defined by (yI,slub)
	*/
	//int selectBranchingIndex2(BranchingDecisions* newBd0, BranchingDecisions* newBd1);

	/* Branch on the variable that has an ideal point with a dimension the closest to the slub
	*/
	int selectBranchingMinScore(BranchingDecisions* newBd0, BranchingDecisions* newBd1, std::vector<BranchingDecisions>& dec0, std::vector<BranchingDecisions>& dec1);

	/* \brief Generate cuts to the node.
	 *
	 * Generate cover inequalities for minimization problems for OB constraints. Add them to the probing model and to the linear relaxation.
	 */
	void generateCuts();

	/* \brief Compute the cover inequalities
	*
	* Compute cover inequalities for minimization problems, and add them to coverCuts, their description. The cover cut is computed by considering variables with the highest to the smallest
	* objective coefficient.
	*/
	void computeCoverCuts();

	/* \brief Generate cuts to the node.
	 *
	 * Generate cover inequalities for minimization problems for OB constraints. Add them to the probing model and to the linear relaxation.
	 */
	void generateCuts(BranchingDecisions& bd);

	/* \brief Compute the cover inequalities
	*
	* Compute cover inequalities for minimization problems, and add them to coverCuts, their description.
	*/
	void computeCoverCuts(BranchingDecisions& bd);

	/* \brief Compute the cover inequalities
	*
	* Compute cover inequalities for minimization problems, and add them to coverCuts, their description. The cover cut is computed by considering variables from the most to the least often fractional
	* in the extreme points of the lower bound set.
	*/
	void computeMofCoverCuts(BranchingDecisions& bd);

	/*! \brief Try to fix integer bounds to each variable by solving LPs and looking at their feasibility
	 *
	 * \param BranchingDecisions newBd0. Describes a potential futur subproblem where a variable is set to 0.
	 * \param BranchingDecisions newBd1. Describes a potential futur subproblem where a variable is set to 1.
	 * \return true if this node is detected as infeasible after checking the LPs.
	 */
	bool fixAndSelectBranchingVariable(BranchingDecisions* newBd0, BranchingDecisions* newBd1, UpperBoundSet* U);

	/*! \brief Try to fix integer bounds to each variable by solving LPs and looking at their feasibility
	 *
	 * \param BranchingDecisions newBd0. Describes a potential futur subproblem where a variable is set to 0.
	 * \param BranchingDecisions newBd1. Describes a potential futur subproblem where a variable is set to 1.
	 * \return true if this node is detected as infeasible after checking the LPs.
	 */
	bool fixAndSelectBranchingVariable1(BranchingDecisions* newBd0, BranchingDecisions* newBd1, UpperBoundSet* U);

	/*! \brief Try to fix integer bounds to each variable by solving LPs and looking at their feasibility
	 *
	 * \param BranchingDecisions newBd0. Describes a potential futur subproblem where a variable is set to 0.
	 * \param BranchingDecisions newBd1. Describes a potential futur subproblem where a variable is set to 1.
	 * \return true if this node is detected as infeasible after checking the LPs.
	 */
	bool fixAndSelectBranchingVariable2(BranchingDecisions* newBd0, BranchingDecisions* newBd1, UpperBoundSet* U);

	/* \brief Tries to fix variables manually by looking at the constraints
	 *
	 * \param BranchingDecisions newBd0. Describes a potential futur subproblem where a variable is set to 0.
	 * \param BranchingDecisions newBd1. Describes a potential futur subproblem where a variable is set to 1.
     */
	bool manualInspection(BranchingDecisions* newBd0, BranchingDecisions* newBd1, std::vector<std::vector<bool>>* toTest);

	/* \brief Tries to fix variable i manually by looking at the constraints
	 *
	 * \param BranchingDecisions newBd0. Describes a potential futur subproblem where a variable is set to 0.
	 * \param BranchingDecisions newBd1. Describes a potential futur subproblem where a variable is set to 1.
	 * \param int i. Index of the variable to fix.
	 */
	bool manualInspection(BranchingDecisions* newBd0, BranchingDecisions* newBd1, std::vector<std::vector<bool>>* toTest, int i, std::vector<bool>* inList, std::list<int>* L);

	/* \brief Tries to fix variable i manually by looking at the constraints
	 *
	 * \param BranchingDecisions newBd0. Describes a potential futur subproblem where a variable is set to 0.
	 * \param BranchingDecisions newBd1. Describes a potential futur subproblem where a variable is set to 1.
	 * \param int i. Index of the variable to fix.
	 * \param bool modified. Set to true if a variable is fixed through preprocessing.
	 */
	bool manualInspection(BranchingDecisions* newBd0, BranchingDecisions* newBd1, int i, bool* modified);

	/* \brief Update the list of variable to potentially fix (L) after variable x_i was fixed in the model.
	 *
     * \param vector of vector of bool, toTest. Used to know which variables should be tested.
	 * \param vector of bool, inList. Vector that states whether each variable is already in the list of variables to explore or not. Used to avoid redundancies.
	 * \param list of int. List of variables to explore for fixing.
	 * \param int i. Index of the variable fixed.
	 */
	void updateListVarToFix(std::vector<std::vector<bool>>* toTest, std::vector<bool>* inList , std::list<int>* L, int i);

	/* \brief Record at which iteration this node is explored.
	 *
	 * \param int it. The number of the iteration.
     */
	void setIteration(int it);

	/* \brief Compute and returns the pseudo-cost of the last branching variable in the current node.
	 */
	double getBranchingVariableScore();

	/* \brief Returns the agregated value of variable i in the current node.
	 */
	double getVariableScoreValue(int i);

	/* \brief returns the last branching variable.
	*/
	int getLastBranchingVariable();

	/* \brief returns the last branching value.
	 * Note: valid for binary problems only
	 */
	int getLastBranchingValue();

	/* Get the percentage of integer extreme point in the lower bound set computed at this node.
	 *
	 * \return the percentage as a double between 0 and 1.
	 */
	double getPercentageIntegralityLB();

	/* Get the percentage of extreme point in the lower bound set computed at this node that satisfy the branching decisions of this node.
	 *
	 * \return the percentage as a double between 0 and 1.
	 */
	double getPercentageFeasibilityLB();

	/*! \brief Compute the rhs of the weighted sum scalarization when the node is created.
	 */
	void computeWsScore();

	/*! \brief Compute the the minimum difference between the rhs of an hyperplane and the w.s. of an lub using the normal vector of the hyperplane,
	 * for each Lub. Then, return the maximal value obtained.
	 */
	void computeMaxMinGapScore();

	/* Compute idx with a full strong branching approach on a w.s. in function of the branching decisions bd.
	 *
	 * \param branchingDecisions* bd. The branching decisions made so far.
	 */
	int getFullStrongBranchingWSIndex(BranchingDecisions* bd);

	int getLargestFractionalReducedCost(BranchingDecisions* bd, std::vector<double>& rc, std::vector<double>& sol0, std::vector<double>& sol1);

	void setUpScoreWS();

	double getScore();

	bool isMorePromising(Node* nd);

	void mergeSlubs(std::list<SLUB*>* S);

	void printCteFixed(std::vector<bool>& free, int j);

	bool operator>(const Node& nd) const;
};