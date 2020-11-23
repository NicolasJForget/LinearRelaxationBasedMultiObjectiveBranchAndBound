/**
* \author Nicolas Forget
*
* This class describe an upper bound set: the incumbent set.
*/

#pragma once
#include <stdlib.h>
#include <list>
#include <vector>
#include "solution.h"
#include "Model.h"
#include "Lub.h"




// ===============================================================================================================================
//							UpperBoundSet
// ===============================================================================================================================

class UpperBoundSet {
private:
	std::list<Solution*> incumbentSet; //!< list of feasible solutions, that represents the imcumbent set.
	std::list<LocalUpperBound> lubSet; //!< set of local upper bounds, that corresponds to the current imcumbet set.
	int lastGeneratedId; //!< last id generated for a local upper bound
	MathematicalModel* lp; //!< a pointer to the problem corresponding to this upper bound set

public:
	/*! \brief Default constructor of an upper bound set.
	 *
	 * This function creates an enpty upper bound set.
	 * \param lp MathematicalModel. Used to get the number of objectives.
	 */
	UpperBoundSet(MathematicalModel& lp);

	/*! \brief Generate a new unique id for a local upper bound.
	 *
	 * This function generates a new unique id for a local upper bound. It simply returns the last id generated + 1.
	 */
	int newId();

	/*! \brief Update the upper bound set with a new point.
	 *
	 * This function updates the upper bound set with a new point P. First, it creates a new Solution corresponding to P, then
	 * add it to the upper bound set if it is not dominated by an existing point and also discard the points that are dominated
	 * by this new solution.
	 * \param P Point. The new point added to the upper bound set. Assumed to be integer.
	 */
	void updateUB(Point& P);

	/*! \brief Prints the upper bound set.
	 *
	 * This function prints the objective vector of each solution of the upper bound set.
	 */
	void print();

	/*! \brief Prints the set of local upper bounds.
	 *
	 * This function prints the set of local upper bounds, by showing their coordindates in the objective space
	 */
	void printLub();

	/*! \brief Update the set of local upper bounds
	 *
	 * This function updates the set of local upper bounds, given the new Solution s added to the upper bound set.
	 * \param s Solution. The new solution added to the upper bound set.
	 */
	void updateLUB(Solution* s);

	/*! \brief Generate a new unique id for a new local upper bound
	 *
	 * This function generates a new id for a new local upper bound.
	 * \return the id as an int.
	 */
	int genNewId();

	/*! \brief Provides a pointer to the list of local upper bounds.
	 *
	 * This function returns a pointer to the list of local upper bounds of this upper bound set.
	 * \return a pointer to the list of local upper bounds.
	 */
	std::list<LocalUpperBound>* getLubs();

	/*! \brief Returns the number of non-dominated points
	 *
	 * This function returns the number of non-dominated points, by returning the size of incumbentSet.
	 * \return number of non-dominated points, as an int.
	 */
	int getNbNonDominatedPoints();

	/*! \brief Reset the UpperBoundSet.
	 *
	 * This function reset the upper bound set by deleting all its components.
	 */
	void reset();
};
