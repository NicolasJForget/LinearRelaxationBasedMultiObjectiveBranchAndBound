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
#include <ctime>
#include <cstdlib>
#include "Model.h"
#include "Node.h"
#include "UB.h"
#include "GlobalConstants.h"
#include "Stat.h"
#include "timer.hpp"
#include "Tree.h"

class BranchAndBound
{
private:
	MathematicalModel lp; //!< a mathematical model, that describes the problem solved
	std::list<Node*> queue; //!< list of non-explored node of the branch-and-bound tree
	TreeManager T; //!< tree, manages the non-explored nodes.
	UpperBoundSet U; //!< upper bound set
	Parameters P; //!< parameters of the branch-and-bound algorithm
	Statistics stat; //!< statistics about the branch-and-bound
	//StatProbing sp; //!< statistics for probing
	std::string inst; //!< the instance solved using this BB
	std::string pbClass; //!< the problem class of the instance solved
	std::vector<std::vector<double>> scoreVariable; //!< the score of each variable, used for pseudo-cost branching (psy in Archetberg)
	std::vector<std::vector<double>> nbFeasibleBranching; //!< the number of node where x_i = v led to a feasible node for each variable i=1,...,n and value v = 0,1 , used for pseudo-cost branching (eta in Archetberg)

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
	 * \param limitSubPb int. Maximum number of subproblem allowed for objective branching. Relevant for LIMITED_OBJEECTIVE_BRANCHING only.
	 */
	void run(int lb, int nodeSel, int varSel, int ob, int valBranch, int timeout, int versionProbing, int limitSubPb, int domiVarFix, int adjustBBWS, int cuts);

	/*! \brief Runs the branch-and-bound algorithm with the given parameters.
	 *
	 * This function runs the branch-and-bound algorithm with the given parameters, which are:
	 * 1) the lower bound set used
	 * 2) the node selection strategy
	 * 3) the branching value selection strategy
	 * 4) the time limit
	 * See GlobalConstants.h for the available identifiers.
	 * \param lb int. The identifier of the lower bound set.
	 * \param nodeSel int. The identifier of the node selection strategy.
	 * \param ob int. The identifier of the objective branching strategy.
	 * \param timeout int. The time limit.
	 */
	void run(int lb, int nodeSel, int valBranch, int timeout);

	/*! \brief Runs the branch-and-bound algorithm with the given parameters.
	 *
	 * This function runs the branch-and-bound algorithm with the given parameters, which are:
	 * 1) the lower bound set used
	 * 2) the node selection strategy
	 * 4) the time limit
	 * See GlobalConstants.h for the available identifiers.
	 * \param lb int. The identifier of the lower bound set.
	 * \param nodeSel int. The identifier of the node selection strategy.
	 * \param ob int. The identifier of the objective branching strategy.
	 * \param timeout int. The time limit.
	 */
	void run(int lb, int nodeSel, int ob, int valBranch, int limitSubPb, int timeout);

	/*! \brief Select the next node to be explored in the tree
	 *
	 * This function return a pointer to the next node to be exlored in the tree, by searching for the appropriate node
	 * in the queue given the node selection rule used.
	 * \return a pointer to a Node.
	 */
	Node* selectNode();

	/*! \brief Update the scores of the variables
	 */
	void updateScores(Node* nd);

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

	/*! \brief Write the statistic of this run to the stat file (results.txt / .csv)
	 *
	 * This function write the statistics of this BB run to the stat file.
	 */
	void writeStatistics();

	/*! \brief Write the statistic of this run to the stat file (results.txt / .csv)
	 *
	 * This function write the statistics of this BB run to the stat file.
	 */
	void writeDepthStatistics();

	/*! \brief Write the statistic of this run to the stat file (results.txt / .csv)
	 *
	 * This function write the statistics of this BB run to the stat file.
	 */
	void writeUB();

	/*! \brief Write the statistic of this run to the stat file (results.txt / .csv)
	 *
	 * This function write the statistics of this BB run to the stat file.
	 */
	void writeStatProbing();

	/*! \brief Reset the branch-and-bound
	 *
	 * This function reset the branch-and-bound by deleting all its elements except the linear program.
	 */
	void reset();

	/* Compute the weight vector used for BEST_BOUND strategy
	 */
	void computeBestBoundDirection();

	void getYnFromFile();

	UpperBoundSet* getUb();
};

