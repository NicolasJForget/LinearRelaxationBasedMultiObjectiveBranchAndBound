/**
* \author Nicolas Forget
*
* This class describe an extreme point of the lower bound set. Used in the linear relaxation (LP relax).
*/

#pragma once
#include <stdlib.h>
#include "Model.h"
#include "Hyperplane.h"
#include <list>
#include <vector>
#include "Point.h"


class LowerBoundSet {
protected:
	MathematicalModel* lp; //!< a pointer to the initial problem

public:
	/*! \brief Default constructor of a lower bound set.
	 *
	 * This function creates a lower bound set with a reference to the corresponding problem (lp)
	 * \param lp MathematicalModel*. A pointer to the problem this lower bound set refers to.
	 */
	LowerBoundSet(MathematicalModel* lp);

	/*! \brief A virtual function, to call the correct computation.
	 *
	 * This function will be implemented for each lower bound set in the subclasses. If the program enters in this compute(),
	 * it throws a message.
	 */
	virtual void compute();

	/*! \brief A virtual function that prints the lower bound set.
	 *
	 * This function will be implemented for each lower bound set in the subclasses. If the program enters in this print(),
	 * it throws a message.
	 */
	virtual void print();
};

// ===============================================================================================================================
//							LinearRelaxation
// ===============================================================================================================================

/**
* \author Nicolas Forget
*
* This class describe the linear relaxation. Note that it actually store $\mathcal{L} + \mathbb{R}^p$. It is computed
* using benson's outer approximation algorithm.
*/

class LinearRelaxation : public LowerBoundSet {
private:
	// Cplex models
	WeightedSumModel* weightedSum; //!< pointer to the cplex model that computes weithed sums
	FeasibilityCheckModel* feasibilityCheck; //!< pointer to the cplex model that computes checks for the feasibility of a point
	DualBensonModel* dualBenson; //!< pointer to the cplex model that computes the dual from benson's method
	FurthestFeasiblePointModel* furthestFeasiblePoint; //!< pointer to the cplex model that computes the points on the boundary of the LB set in benson's method

	// Needed tools
	std::vector <double> antiIdealPoint; //!< anti-ideal point (it is actually a point weakly dominated by the anti-ideal)
	std::vector <double> interiorPoint; //!< interior point from benson's method, such that it dominates antiIdealPoint and is weakly dominated by at least one feasible point.
	std::vector <Hyperplane*> boundingBox; //!< list of p hyperplanes that defines the bounding box of the linear relaxation. It corresponds to the dominated cone defined by antiIdeal
	bool warmstarted; //!< true if the linear relaxation is warmstarted, false otherwise.

	// Representations
	std::list <Hyperplane*> facets; //!< list of hyperplanes, that defines the facet of the linear relaxation.
	std::list <Point*> extrPoints; //!< list of points, that contains the extreme points of the linear relaxation.

	// Statistics
	int nbIterations;

public:
	/*! \brief Default constructor of a linear relaxation.
	 *
	 * This function initialize the linear relaxation to a simplex as defined in benson's outer approximation algorithm.
	 * Note that it calls the initialize function, where the simplex is built. This constructor is called if the linear
	 * relaxation is not warmstarted.
	 * \param lp MathematicalModel*. A pointer to the problem this lower bound set refers to.
	 * \param ws WeightedSumModel*. A pointer to the model that solves weighted sums.
	 * \param feas FeasibilityCheckModel*. A pointer to the cplex model that computes checks for the feasibility of a point.
	 * \param db DualBensonModel*. A pointer to the cplex model that computes the dual from benson's method.
	 * \param bvpFurthestFeasiblePointModel*. A pointer to the cplex model that computes the points on the boundary of the LB set in benson's method,
	 */
	LinearRelaxation(MathematicalModel* lp, WeightedSumModel* ws, FeasibilityCheckModel* feas, DualBensonModel* db, FurthestFeasiblePointModel* bvp);

	/*! \brief The copy constructor for a linear relaxation
	*
	 * This function initialize a new linear relaxation, such that the polyhedron it defines is the same than LPrelax.
	 * It in particular take care that new memory is allocated for all the pointers and that these pointers does not link to
	 * objects, or to objects used in another LinearRelaxation. Also, as this LinearRelaxation is not initialized with a
	 * simplex as in Benson's outer approximation algorithm, it enables warm-starting.
	 * \param LPrelax LinearRelaxation. The LinearRelaxation this object is a copy of.
	 */
	LinearRelaxation(LinearRelaxation& LPrelax);


	/*! \brief Initialise the linear relaxation by building a simplex that contains it.
	 *
	 * This function initialize the linear relaxation by creating a simplex that contains it.
	 * \param linprog MathematicalModel*. A pointer to the problem this lower bound set refers to.
	 */
	void initialize(MathematicalModel& lp);

	/*! \brief This function computes the linear relaxation
	 *
	 * This function computes the linear relaxation by performing steps from Benson's outer approximation algorithm.
	 * It assumes that it has been initialized by a polyhedron that is an outer approximation of the lienar relaxation.
	 */
	virtual void compute();

	/*! \brief Update the representations of the linear relaxation by adding a new hyperplane.
	 *
	 * This function computes the new polyhedron obtained by intersecting by the half-space defined by H. It updates
	 * both representations, i.e. extrPoints and facets. It is done by using the on-line vertex enumeration algorithm
	 * from Chen, Hansen and Jaumard (1993).
	 * \param H Hyperplane*. A pointer to the new hyperplane added to the representation.
	 */
	void updatePolyhedron(Hyperplane* H);

	/*! \brief Identify the actual extreme point of the lower bound set
	 *
	 * This function identifies the actual extreme points of the lower bound set, i.e. not those that are located on the
	 * bounding box. To do so, it just check for each extreme point if it share one component with antiIdealPoint.
	 */
	void filterExtremePoints();

	/*! \brief Prints the lower bound set.
	 *
	 * This function prints the linear relaxation, both as a list of hyperplanes and a list of extreme points.
	 */
	virtual void print();
};