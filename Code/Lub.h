#pragma once
#include <stdlib.h>
#include <list>
#include <vector>
#include "solution.h"
#include "Model.h"
#include "Hyperplane.h"

// ===============================================================================================================================
//							LocalUpperBound
// ===============================================================================================================================

class LocalUpperBound {
private:
	std::vector<int> coordinates; //!< coordinates of the local upper bound in the objective space.
	std::vector<std::list<Solution*>> definingPoints; //!< list of defining points for each objective.
	int id; //!< a unique id for each local upper bound. Used to recognize the non-dominated ones prior to the dominance tests.
	bool discarded;

public:
	/*! \brief Default constructor of a local upper bound.
	 *
	 * This function creates a local upper bounds with empty coordinates and defining points.
	 */
	LocalUpperBound();

	/*! \brief Constructor of a local upper bound that is used in the initilalization of UB set.
	 *
	 * This function creates a local upper bounds with coordinates (M,...,M) and no defining point. Note that M is fixed
	 * to INT_MAX, the maximal possible value of an int.
	 * \param lp MathematicalModel. Used to get the number of objectives.
	 */
	LocalUpperBound(MathematicalModel& lp);

	/*! \brief Constructor of a local upper bound defined by an old one u and a new point y on objective j.
	 *
	 * This function creates a local upper bounds with the same coordinates and defining points as u, except on objective
	 * j, where z is used as defining point and defines the jth coordinates too.
	 * \param y Solution. The new solution that defines the local upper bound.
	 * \param u LocalUpperBound. The old local upper boud used to defines the other p - 1 coordinates.
	 * \param j int. The objective defined by the new solution y.
	 * \param id int. Id of the new local upper bound.
	 */
	LocalUpperBound(Solution* y, LocalUpperBound& u, int j, int id);

	/*! \brief Prints the local upper bound.
	 *
	 * This function prints the coordinates of this local upper bound
	 */
	void print();

	/*! \brief Check whether this Local Upper Bound is dominated by a Solution y.
	 *
	 * This function checks whether the Solution y dominates this Local Upper Bound, by doing component-wise comparisons.
	 * \param y Solution. The Solution used for comparison.
	 * \return true if this local upper bound is strictly dominated by y
	 */
	bool isStrictlyDominated(Solution* y);

	/*! \brief Check whether this Local Upper Bound is defined by a Solution y on a given component k.
	 *
	 * This function checks whether the Solution y defines coordinates k of this Local Upper Bound, by checking
	 * whether this local upper bound is equal to the solution on coordinate k, and if all other coordinates
	 * are pairwise strictly lower than those of this local upper bound.
	 * \param y Solution. The Solution used for comparison.
	 * \param k int. The index of the coordinate.
	 * \return true if this local upper bound is defined by y
	 */
	bool isDefinedBy(Solution* y, int k);

	/*! \brief Add a new defining component of this local upper bound on coordinate k
	 *
	 * This function adds the solution y as new defining component of this local upper bound on objective k.
	 * \param y Solution. The Solution used for comparison.
	 * \param k int. The index of the objective.
	 */
	void addDefiningComponent(Solution* y,int k);

	/*! \brief Filter the discarded defining solutions on coordinate k
	 *
	 * This function remove from the list of defining solutions of coordinate k the solution that are discarded, i.e. dominated.
	 * \param k int. The index of the objective.
	 */
	void filterDefiningComponent(int k);

	/*! \brief Compute the critical value z^{max}_j in the computation of the local upper bounds
	 *
	 * This function computes the critical value z^{max}_j from the method of Klamroth, Lacour and Vanderpooten from 2015,
	 * used to check whether a new local upper bound should be computed in direction of objective j.
	 * \param j int. The index of the coordinate we are looking at.
	 * \return the value of the critical value.
	 */
	int computeCriticalValue(int j);

	/*! \brief Set this local upper bound as discarded.
	*
	* This function notes this local upper bound as discarded, by setting discarded to true.
	*/
	void becomesDiscarded();

	/*! \brief Check whether the local upper bound is already discarded.
	*
	* This function check whether the local upper bound is already discarded or not, by looking at discarded.
	* \return true if the local upper bound is already discarded.
	*/
	bool isDiscarded();

	/*! \brief Return the number of objectives.
	*
	* This function returns the number of coordinates this local upper bound have.
	*/
	int get_dim();

	/*! \brief Return the value of a specific coordinate.
	*
	* This function returns the value of a specific coordinate.
	* \param k int. The index of the coordinate.
	* \return the value of the coordinate
	*/
	int get_coordinate(int k);

	/*! \brief Return the id of the local upper bound.
	*
	* This function returns the id of this local upper bounx
	* \return the value of the id, as an int
	*/
	int get_id();

	/*! \brief Return the defining points of a specific coordinate.
	*
	* This function returns a list of pointers to the solutions that defines coordinate k.
	* \param k int. The index of the coordinate.
	* \return a pointer to the list of defining Solutions.
	*/
	std::list<Solution*>* get_definingPoints(int k);

	/*! \brief Check whether a local upper bound is located above the hyperplane.
	*
	* This function checks whether this local upper bound is located above (or located on) the hyperplane H. It substracts
	* the greatest common divisor for each objective, as no new integer solution can be found in these areas.
	* \param H Hyperplane. A pointer to the hyperplane used for comparison.
	* \param P Parameters. A pointer to the parameters of the algorithm, to access the GCD of each objective.
	* \return a pointer to the list of defining Solutions.
	*/
	bool above(Hyperplane* H, Parameters* P);
	bool above(Hyperplane* H);

	/*! \brief Compute the weigthed sum of the components in the objective space.
	 *
	 * \param vector of double l. The weight vector.
	 */
	double getWeightedSum(std::vector<double> l);

	bool isRedundant(std::list<LocalUpperBound>& NU);
};

