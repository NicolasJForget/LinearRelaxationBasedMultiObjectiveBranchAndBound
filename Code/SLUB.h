#pragma once
#include <stdlib.h>
#include <list>
#include <vector>
#include "GlobalConstants.h"
#include "Lub.h"
#include "Point.h"

class SLUB
{
private:
	std::vector<int> coord; //!< coordinates of the super local upper bound
	int dim; //!< number of dimensions, i.e. of coordinates

public:
	/*! \brief Constructor of slub given branching decisions.
	 *
	 * This function creates a slub, by assigning the slub field of the BranchingDecision to coord.
	 * \param branchDec BranchingDecision. The branching decision from which we read the slub field.
	 */
	SLUB(BranchingDecisions& branchDec, Parameters* param);

	/*! \brief Constructor of slub given a local upper bound.
	 *
	 * This function creates a slub, by assigning the coordinates of the local upper bound to the coord field.
	 * \param lub LocalUpperBound*. A pointer to the local upper bound from which we read the coordinates.
	 */
	SLUB(LocalUpperBound& lub);

	SLUB(std::vector<int>& y);

	/*! \brief Constructor for an empty slub in dimension p.
	 *
	 * This function creates a slub in dimension p with the value INT_MIN for each component.
	 * \param p int. The dimension of the slub.
	 */
	SLUB(int p);

	/*! \brief Merges two slub.
	 *
	 * This function merges this slub with the slub given in parameter. It replaces coord of this slub by the nadir point of
	 * this slub and the other.
	 * \param slub SLUB. The slub that has to be merged with this one.
	 */
	void merge(SLUB& slub);

	/*! \brief Merges a slub and a lub.
	 *
	 * This function merges this slub with the lub given in parameter. It replaces coord of this slub by the nadir point of
	 * this slub and the lub.
	 * \param lub LocalUpperBound. The lub that has to be merged with this slub.
	 */
	void merge(LocalUpperBound& lub);

	/*! \brief Return the value of a specific coordinate.
	*
	* This function returns the value of a specific coordinate.
	* \param k int. The index of the coordinate.
	* \return the value of the coordinate
	*/
	int get_coordinate(int k);

	/*! \brief Check whether this SLUB is dominated the point pts.
	 *
	 * This function check whether this SLUB is dominated the point pts, by doing component-wise comparisons
	 * \param pts Point*. A pointer to the point used for comparison.
	 * \return true if this slub is dominated by the point pts, i.e. pts is located in the cone defined by this SLUB.
	 */
	bool dominated(Point* pts);
	bool dominated(LocalUpperBound& u);
	bool strictlyDominated(Point* pts);
	bool dominatedFix(Point* pts);

	/* \brief Print the slub
	*/
	void print();

	/*! \brief Compute the distance between this SLUB and s.
	 *
	 * This function computes the euclidian distance between the two SLUBs.
	 * \param s SLUB*. The SLUB used for comparison.
	 * \return the distance, as a double.
	 */
	double distance(SLUB* s);
};

