/**
* \author Nicolas Forget
*
* This class describe a feasible solution of the initial problem.
*/

#pragma once
#include <stdlib.h>
#include <list>
#include <vector>
#include "Point.h"

class Solution {
private:
	std::vector<int> x; //!< solution vector of the solution. Always integer.
	std::vector<int> y; //!< objective vector of the solution. Always integer, as the objective coefficients are integer.

public:
	/*! \brief Default constructor of a solution.
	 *
	 * This function creates an empty solution.
	 */
	Solution();

	/*! \brief Construct a solution from an extreme point of the LinearRelaxation
	 *
	 * This function creates a solution from an extreme point of the LinearRelaxation, by extracting its pre-image and
	 * objective vector.
	 * \param pts Point. The point the solution is created from.
	 */
	Solution(Point& pts);
};