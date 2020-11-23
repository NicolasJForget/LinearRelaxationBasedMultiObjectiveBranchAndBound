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
#include "LB.h"
#include "UB.h"
#include "GlobalConstants.h"
#include "timer.hpp"

class Node
{
private:
	MathematicalModel* P; //!< a pointer to the mathematical formulation to the initial problem.
	Parameters* param; //!< a pointer to the parameters of the algorithm
	LowerBoundSet* LB; //!< Lower bound set of this node. As a pointer for inheritance and virtual functions to work properly.
	double score; //!< a score, used for node selection
	BranchingDecisions branchingDecision; //!< a data structure that save the branching decisions made in this node. See GlobalConstants.h for more details.
	int splittingIndex; //!< the index of the last variable splitted, i.e. the branching decision that led to the creation of this node.
	//int status; //!< defines the status of the node
	Timer timeLBComputation;
	Timer timeDominanceTest;
	Timer timeUpdateUB;

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
	 * \param par Parameters. The parameters of the algorithm.
	 */
	Node(MathematicalModel* lp, Parameters* par);

	/*! \brief Creates a node given a splitting index and a bound.
	 *
	 * This function creates a new node that is a sub-problem of node nd, splitted on variable x[i].
	 * \param nd Node*. A pointer to the node it is created from.
	 * \param index int. The index of the variable to split.
	 * \param bound int. States the sign of the new constraint: are we adding an upper or a lower bound value on a variable?
	 * \param val int. The value of the bound added on the variable, i.e. the right-hand side of the new branching constraint.
	 */
	Node(Node& nd, int index, int bound, int val);

	/*! \brief Process the node
	 *
	 * This function process the node by computing a lower bound set, updating the lower bound set, and detecting 
	 * if the node can be fathomed.
	 * \param U UpperBoundSet. The upper bound set of the B&B.
	 * \param stat Statistics*. Pointer to statistics.
	 */
	void process(UpperBoundSet& U);

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
	 */
	void splitOS(std::list<Node*>* Q);

	/*! \brief Split the problem in the variable space
	 *
	 * This function splits the problem contained in this node in the variable space. Store the new nodes created in the list Q.
	 * \param Q pointer to a list of pointers of Nodes. The list in which the new nodes created are stored.
	 */
	void splitVS(std::list<Node*>* Q);

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
	void getStatistics(Statistics* stat);

	/*! \brief Rreturn the status of the node
	 *
	 * Returns the status of the node (infeasible, optimal, dominated...). See GlobalConstants.h for status.
	 * \return status of the node, as an integer
	 */
	int getStatus();
};

