/**
* \author Nicolas Forget
*
* This class describes the branch-and-bound algorithm.
*/

#pragma once
#include <stdlib.h>
#include <list>
#include <vector>
#include <string>
#include "Model.h"
#include "Node.h"
#include "UB.h"
#include "GlobalConstants.h"
#include "timer.hpp"

class BranchAndBound
{
private:
	MathematicalModel lp; //!< a mathematical model, that describes the problem solved
	std::list<Node*> queue; //!< list of non-explored node of the branch-and-bound tree
	UpperBoundSet U; //!< upper bound set
	Parameters P; //!< parameters of the branch-and-bound algorithm
	Statistics stat; //!< statistics about the branch-and-bound

public:
	/*! \brief Default constructor of a branch-and-bound
	 *
	 * This function creates a branch-and-bound for the problem described in file.
	 * \param instance string. The path to the file that contains the instance solved by this branch and bound
	 */
	BranchAndBound(std::string instance);

	/*! \brief Runs the branch-and-bound algorithm with the given parameters.
	 *
	 * This function runs the branch-and-bound algorithm with the given parameters, which are:
	 * 1) the lower bound set used
	 * 2) the node selection strategy
	 * 3) the variable selection strategy
	 * 4) the objective branching strategy
	 * See GlobalConstants.h for the available identifiers.
	 * \param lb int. The identifier of the lower bound set.
	 * \param nodeSel int. The identifier of the node selection strategy.
	 * \param varSel int. The identifier of the variable selection strategy.
	 * \param ob int. The identifier of the objective branching strategy.
	 */
	void run(int lb, int nodeSel, int varSel, int ob);

	/*! \brief Select the next node to be explored in the tree
	 *
	 * This function return a pointer to the next node to be exlored in the tree, by searching for the appropriate node
	 * in the queue given the node selection rule used.
	 * \return a pointer to a Node.
	 */
	Node* selectNode();

	/*! \brief Print out YN
	 *
	 * This function prints the set of non-dominated point, by printing the upper bound set U.
	 */
	void printYN();

	/*! \brief Print out statistics
	 *
	 * This function prints various statistics about this branch-and-bound run.
	 */
	void printStatistics();

	/*! \brief Reset the branch-and-bound
	 *
	 * This function reset the branch-and-bound by deleting all its elements except the linear program.
	 */
	void reset();
};

