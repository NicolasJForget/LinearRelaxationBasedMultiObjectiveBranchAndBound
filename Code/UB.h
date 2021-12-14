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
#include "timer.hpp"
#include "SLUB.h"

//class SLUB;


// ===============================================================================================================================
//							UpperBoundSet
// ===============================================================================================================================

class UpperBoundSet {
private:
	std::list<Solution*> incumbentSet; //!< list of feasible solutions, that represents the imcumbent set.
	std::list<LocalUpperBound> lubSet; //!< set of local upper bounds, that corresponds to the current imcumbet set.
	int lastGeneratedId; //!< last id generated for a local upper bound
	MathematicalModel* lp; //!< a pointer to the problem corresponding to this upper bound set
	Timer cpuSol; //!< a timer to compute the time of creation of each solution
	FeasibilityCheckModel* fcm; //!< a pointer to the feasibility check model
	ProbingModel* pm; //!< a pointer to the probing model

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

	/*! \brief Update the upper bound set with a new point.
	 *
	 * This function updates the upper bound set with a new point P. First, it creates a new Solution corresponding to P, then
	 * add it to the upper bound set if it is not dominated by an existing point and also discard the points that are dominated
	 * by this new solution.
	 * \param P Point. The new point added to the upper bound set. Assumed to be integer.
	 */
	void updateUB(BranchingDecisions* bd);

	/*! \brief Update the upper bound set with a new candidate objective vector. Also check integer feasibility.
	 *
	 * This function updates the upper bound set with a new point P. First, it creates a new Solution corresponding to P, then
	 * add it to the upper bound set if it is not dominated by an existing point and also discard the points that are dominated
	 * by this new solution.
	 * \param P Point. The new point added to the upper bound set. Assumed to be integer.
	 */
	void updateUB(std::vector<double>& y);

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

	/*! \brief Provides a pointer to the list of local upper bounds.
	 *
	 * This function returns a pointer to the list of local upper bounds of this upper bound set.
	 * \return a pointer to the list of local upper bounds.
	 */
	std::list<Solution*>* getYN();

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

	/*! \brief Check the largest difference between the ws value w, and the weighted sum of a lub in cone s with weights l.
	 *
	 * \param double w. The weighted sum value.
	 * \param vector of double l. The weight vector of the objective functions.
	 * \param branchingDecisions bd. Include the slub s.
	 */
	double getLargestGap(double w, std::vector<double>& l, BranchingDecisions* bd);

	/* Test whether the weighted sum with weights l of the lubs included in the OB cone defined in bd is lower or greater
	 * than value ws. Return true if at least one of them has a greater weighted sum value, false otherwise.
	 *
	 * \param BranchingDecisions* bd. The branching decisions that define the cone we are located in in the objective space.
	 * \param vector of double l. The weight vector of the weighted sums considered.
	 * \param double ws. The value used for comparison of the weighted sums.
	 */
	bool testWeightedSumValue(BranchingDecisions* bd, std::vector<double>& l, double w);

	/* Get the largest value for a local upper bound with weight l and in the slub defined in bd.
	 *
 	 * \param BranchingDecisions* bd. The branching decisions that define the cone we are located in in the objective space.
	 * \param vector of double l. The weight vector of the weighted sums considered.
     */
	double getLargestWeightedSumValue(BranchingDecisions* bd, std::vector<double>& l);

	std::vector<double> getLubWeigth(BranchingDecisions* bd); // NOT FINISHED !!

	void setCplexModels(FeasibilityCheckModel* ff, ProbingModel* pp);

	void startTimerUB();

	void stopTimerUB();

	void addSol(Solution* s);

	std::list<Solution*>* getIncumbentSet();
};
