/**
* \author Nicolas Forget
*
* This class describe an extreme point of the lower bound set. Used in the linear relaxation (LP relax).
*/

#pragma once
#include <stdlib.h>
#include <list>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <vector>
#include "Model.h"
#include "Hyperplane.h"
#include "Point.h"
#include "UB.h"
#include "GlobalConstants.h"
#include "Stat.h"
#include "SLUB.h"
#include "timer.hpp"


class LowerBoundSet {
protected:
	MathematicalModel* lp; //!< a pointer to the initial problem
	int status; //!< status of the lower bound set
	BranchingDecisions* branchDec; //!< a pointer to the branching decision of the node it is contained in
	int iteration; //!< number of the iteration of the B&B tree

public:
	/*! \brief Default constructor of a lower bound set.
	 *
	 * This function creates a lower bound set with a reference to the corresponding problem (lp)
	 * \param lp MathematicalModel*. A pointer to the problem this lower bound set refers to.
	 * \param branchingDec BranchingDecisions. The data structure that describes the branching decisions to apply to this LB set.
	 */
	LowerBoundSet(MathematicalModel* lp, BranchingDecisions* branch);

	virtual ~LowerBoundSet();

	/*! \brief Return the current status of the lower bound set.
	 *
	 * This function returns the current status of the lower bound set. See GlobalConstants.h for the possible status.
	 * \return the identifier of the status, as an int.
	 */
	int getStatus();

	/*! \brief A virtual function, to call the correct computation.
	 *
	 * This function will be implemented for each lower bound set in the subclasses. If the program enters in this compute(),
	 * it throws a message.
	 */
	virtual void compute(UpperBoundSet* U, std::list<LocalUpperBound*> lubDomi);

	/*! \brief A virtual function, to call the correct computation.
	 *
	 * This function will be implemented for each lower bound set in the subclasses. If the program enters in this compute(),
	 * it throws a message.
	 */
	virtual void computeFull();

	/*! \brief A virtual function that prints the lower bound set.
	 *
	 * This function will be implemented for each lower bound set in the subclasses. If the program enters in this print(),
	 * it throws a message.
	 */
	virtual void print();

	/*! \brief A virtual function that updates the upper bound set with integer solutions from the lower bound set.
	 *
	 * This function search for new integer solutions in the lower bound set, and update the upper bound set U if they
	 * are non-dominated with respect to the points in U. If the program enters this gatherIntegerSolutions(), it throws
	 * a message.
	 * \param U UpperBoundSet. The upper bound set updated.
	 */
	virtual void gatherIntegerSolutions(UpperBoundSet& U);

	/*! \brief A virtual function that adjust the bounds of the variable given the bounds in the node nd.
	 *
	 * This function adjust the bounds of the variables in the Cplex models according to the branching decisions made
	 * in node nd. It also adjusts according to the objective branching constraints. Throw an error if the algorithm
	 * enters in this version of the method.
	 */
	virtual void applyBranchingDecisions(); //BranchingDecisions& branchingDec

	/*! \brief A virtual function that checks whether this lower bound set is dominated by the upper bound set U.
	 *
	 * This function checks whether the upper bound set U dominates this lower bound set. Throw an error if the algorithm
	 * enters in this version of the method. Depending on the parameters of the algorithm, all the local upper bounds may
	 * have to be tested for dominance, even though it is already determined that the node cannot be fathomed by dominance.
	 * \param U UpperBoundSet. The upper bound set used to do the dominance test. Note that the local upper bounds are shifted
	 * of -1 for each component when tested for dominance.
	 * \param param Parameters. Used to check whether all the local upper bounds have to be tested.
	 * \param ndLub list of LocalUpperBound. List of ids of the local upper bounds known as non-dominated by this LB set.
	 */
	virtual void applyDominanceTest(UpperBoundSet& U, Parameters* param, std::list<int>& ndLub, std::list<LocalUpperBound*>* domiLub);


	/*! \brief A virtual function that checks whether this lower bound set the super local upper bound a point y.
	 *
	 * This function checks whether the point y is dominated by the lower bound set. Throw an error if the algorithm
	 * enters in this version of the method. Returns true if this is the case.
	 * \param U UpperBoundSet. The upper bound set used to do the dominance test.
	 * \return true if y is dominated by the lower bound set.
	 */
	virtual bool dominates(std::vector<int>& y);

	/*! \brief Set the number of the iteration.
	 *
	 * \param it int. The number of the iteration.
	 * \return true if y is dominated by the lower bound set.
	 */
	void setIteration(int it);
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
	bool firstGeneration; //!< true if this linear relaxation object is solved for the first time. False if it has been solved at least once (i.e. compute has been called at least once).
	bool checkPointDestroyed; //!< true if the check point has been destroyed, which means that a new one should be extracted
	bool attemptedDeleteNewHpp;

	// Representations
	std::list <Hyperplane*> facets; //!< list of hyperplanes, that defines the facet of the linear relaxation.
	std::list <Point*> extrPoints; //!< list of points, that contains the extreme points of the linear relaxation.
	std::list <Point*> extrRays; //!< list of points, that contains the extreme rays of the linear relaxation.

	// Statistics
	int nbIterations;
	//LPRelaxStat stat;
	int call; //!< for debug purpose - number of calls to lb set solved so far
	//int iteration;
	Statistics* S;

public:

	/* \brief Destructor of a LinearRelaxation
	 *
	 * This functions destroy the linear relaxation. It in particular desallocates all the points and hyperplanes, which are
	 * specific to this lower bound set, unlike the Cplex models.
	 */
	virtual ~LinearRelaxation();

	/*! \brief Default constructor of a linear relaxation, with models built externally
	 *
	 * This function initialize the linear relaxation to a simplex as defined in benson's outer approximation algorithm.
	 * Note that it calls the initialize function, where the simplex is built. This constructor is called if the linear
	 * relaxation is not warmstarted.
	 * \param lp MathematicalModel*. A pointer to the problem this lower bound set refers to.
	 * \param ws WeightedSumModel*. A pointer to the model that solves weighted sums.
	 * \param feas FeasibilityCheckModel*. A pointer to the cplex model that computes checks for the feasibility of a point.
	 * \param db DualBensonModel*. A pointer to the cplex model that computes the dual from benson's method.
	 * \param bvpFurthestFeasiblePointModel*. A pointer to the cplex model that computes the points on the boundary of the LB set in benson's method,
	 * \param branchingDec BranchingDecisions. The data structure that describes the branching decisions to apply to this LB set.
	 * \param S Statistics*. A pointer to the statistics of the BB
	 */
	LinearRelaxation(MathematicalModel* lp, WeightedSumModel* ws, FeasibilityCheckModel* feas, DualBensonModel* db, FurthestFeasiblePointModel* bvp, BranchingDecisions* branch, Statistics* S, Parameters* param);

	/*! \brief Default constructor of a linear relaxation, with models built internally.
	 *
	 * This function initialize the linear relaxation to a simplex as defined in benson's outer approximation algorithm.
	 * Note that it calls the initialize function, where the simplex is built. This constructor is called if the linear
	 * relaxation is not warmstarted.
	 * \param lp MathematicalModel*. A pointer to the problem this lower bound set refers to.
	 * \param branchingDec BranchingDecisions. The data structure that describes the branching decisions to apply to this LB set.
	 * \param S Statistics*. A pointer to the statistics of the BB
	 */
	LinearRelaxation(MathematicalModel* lp, BranchingDecisions* branch, Statistics* S, Parameters* param);

	/*! \brief The copy constructor for a linear relaxation
	*
	 * This function initialize a new linear relaxation, such that the polyhedron it defines is the same than LPrelax.
	 * It in particular take care that new memory is allocated for all the pointers and that these pointers does not link to
	 * objects, or to objects used in another LinearRelaxation. Also, as this LinearRelaxation is not initialized with a
	 * simplex as in Benson's outer approximation algorithm, it enables warm-starting.
	 * \param LPrelax LinearRelaxation. The LinearRelaxation this object is a copy of.
	 * \param branchingDec BranchingDecisions. The data structure that describes the branching decisions to apply to this LB set.
	 */
	LinearRelaxation(LinearRelaxation* LPrelax, BranchingDecisions* branch, Parameters* param);


	/*! \brief Initialise the linear relaxation by building a simplex that contains it.
	 *
	 * This function initialize the linear relaxation by creating a simplex that contains it.
	 * \param linprog MathematicalModel*. A pointer to the problem this lower bound set refers to.
	 */
	void initialize(); //MathematicalModel& lp

	/*! \brief Initialise the linear relaxation by building a simplex that contains it.
	 *
	 * This function initialize the linear relaxation by creating a simplex that contains it.
	 * \param linprog MathematicalModel*. A pointer to the problem this lower bound set refers to.
	 */
	void initializeFull(); //MathematicalModel& lp

	/*! \brief This function computes the linear relaxation
	 *
	 * This function computes the linear relaxation by performing steps from Benson's outer approximation algorithm.
	 * It assumes that it has been initialized by a polyhedron that is an outer approximation of the lienar relaxation.
	 */
	virtual void compute(UpperBoundSet* U, std::list<LocalUpperBound*> lubDomi);

	/*! \brief This function computes the linear relaxation
	 *
	 * This function computes the linear relaxation by performing steps from Benson's outer approximation algorithm.
	 * It assumes that it has been initialized by a polyhedron that is an outer approximation of the lienar relaxation.
	 */
	virtual void computeFull();

	/*! \brief Update the representations of the linear relaxation by adding a new hyperplane.
	 *
	 * This function computes the new polyhedron obtained by intersecting by the half-space defined by H. It updates
	 * both representations, i.e. extrPoints and facets. It is done by using the on-line vertex enumeration algorithm
	 * from Chen, Hansen and Jaumard (1993).
	 * \param H Hyperplane*. A pointer to the new hyperplane added to the representation.
	 * \param ptsDiscardedByH Point*. A pointer to the point that is cut out of the lower bound set by the new hyperplane H.
	 */
	void updatePolyhedron(Hyperplane* H, Point* ptsDiscardedByH);

	/*! \brief Update the representations of the linear relaxation by adding a new hyperplane.
	 *
	 * This function computes the new polyhedron obtained by intersecting by the half-space defined by H. It updates
	 * both representations, i.e. extrPoints and facets. It is done by using the on-line vertex enumeration algorithm
	 * from Chen, Hansen and Jaumard (1993).
	 * \param H Hyperplane*. A pointer to the new hyperplane added to the representation.
	 * \param ptsDiscardedByH Point*. A pointer to the point that is cut out of the lower bound set by the new hyperplane H.
	 */
	void updatePolyhedron2(Hyperplane* H, Point* ptsDiscardedByH);

	/*! \brief Update the representations of the linear relaxation by adding a new hyperplane.
	 *
	 * This function computes the new polyhedron obtained by intersecting by the half-space defined by H. It updates
	 * both representations, i.e. extrPoints and facets. It is done by using the on-line vertex enumeration algorithm
	 * from Chen, Hansen and Jaumard (1993).
	 * \param H Hyperplane*. A pointer to the new hyperplane added to the representation.
	 * \param ptsDiscardedByH Point*. A pointer to the point that is cut out of the lower bound set by the new hyperplane H.
	 */
	void updatePolyhedron3(Hyperplane* H, Point* ptsDiscardedByH);

	/*! \brief Update the representations of the linear relaxation by adding a new hyperplane.
	 *
	 * This function computes the new polyhedron obtained by intersecting by the half-space defined by H. It updates
	 * both representations, i.e. extrPoints and facets. It is done by using the on-line vertex enumeration algorithm
	 * from Chen, Hansen and Jaumard (1993).
	 * \param H Hyperplane*. A pointer to the new hyperplane added to the representation.
	 * \param ptsDiscardedByH Point*. A pointer to the point that is cut out of the lower bound set by the new hyperplane H.
	 */
	void updatePolyhedronFull(Hyperplane* H, Point* ptsDiscardedByH);

	/*! \brief Correct the set of facet-defining hyperplanes by removing only face-defining ones.
	 *
	 * CAN BE MORE EFFICIENT BY FILTERING EACH TIME A NEW REDUNDANT HPP IS FOUND.
	 * This function checks whether there are any hyperplanes remaning that defines a face (and not a facet) of the LB set.
	 * If one is found, it is removed. Once all hyperplanes have been checked, the adjacency list of the extreme points are
	 * consequently corrected.
	 */
	void correctHyperplanes();

	/*! \brief Identify the actual extreme point of the lower bound set
	 *
	 * This function identifies the actual extreme points of the lower bound set, i.e. not those that are located on the
	 * bounding box. To do so, it just check for each extreme point if it share one component with antiIdealPoint.
	 */
	void filterExtremePoints();

	/*! \brief Check whether the problem solved is feasible.
	 *
	 * This function checks whether the linear relaxation is feasible. It calls a FeasibilityCheckModel on the
	 * anti-ideal point.
	 * \return true is it is feasible, false otherwise.
	 */
	bool isFeasible();

	/*! \brief Prints the lower bound set.
	 *
	 * This function prints the linear relaxation, both as a list of hyperplanes and a list of extreme points.
	 */
	virtual void print();

	/*! \brief Updates the upper bound set with integer solutions found in the linear relaxation
	 *
	 * This function check whether each extreme point that is not located on the bounding box is integer. If yes, it is
	 * added to the upper bound set U if it is not dominated by any point from U. Also, the points from U dominated by
	 * the new feasible points found are discarded.
	 * \param U UpperBoundSet. The upper bound set updated.
	 */
	virtual void gatherIntegerSolutions(UpperBoundSet& U);

	/*! \brief A virtual function that adjust the bounds of the variable given the bounds in the node nd.
	 *
	 * This function adjust the bounds of the variables in the Cplex models according to the branching decisions made
	 * in node nd. It also adjusts according to the objective branching constraints.
	 */
	virtual void applyBranchingDecisions(); // BranchingDecisions& branchingDec

	/*! \brief A virtual function that checks whether this lower bound set is dominated by the upper bound set U.
	 *
	 * This function checks whether the upper bound set U dominates this lower bound set. This is done by checking whether there
	 * exists at least one local upper bound located above each hyperplane of the linear relaxation. If this is not the case,
	 * the node can be fathomed by dominance. Depending on the parameters of the algorithm, all the local upper bounds may
	 * have to be tested for dominance, even though it is already determined that the node cannot be fathomed by dominance.
	 * \param U UpperBoundSet. The upper bound set used to do the dominance test.
	 * \param param Parameters. Used to check whether all the local upper bounds have to be tested.
	 */
	virtual void applyDominanceTest(UpperBoundSet& U, Parameters* param, std::list<int>& ndLub, std::list<LocalUpperBound*>* domiLub);

	/*! \brief A virtual function that checks whether this lower bound set the super local upper bound a point y.
	 *
	 * This function checks whether the point y is dominated by the lower bound set. Throw an error if the algorithm
	 * enters in this version of the method. Returns true if this is the case. Specific the the linear relaxation: 
	 * it is checked on which side of each hyperplane that constitute the LB set the point y is located. If it is not
	 * above at least one of these, the point is not dominated and the function returns false.
	 * \param U UpperBoundSet. The upper bound set used to do the dominance test.
	 * \return true if y is dominated by the lower bound set.
	 */
	virtual bool dominates(std::vector<int>& y);

	/*! \brief Computes the index of the most often fractional variable among the extreme points.
	 *
	 * This functions returns the index of the most often fractional variable in the lower bound set located in the area
	 * defined by the super local upper bound slub. This is done by checking each extreme point that dominates slub.
	 * If there is no fractional variable, it returns the index of the most fractional in average (i.e. avg value closest
	 * to 0.5). If there is no point in the given area, it returns the index of the first free variable it finds.
	 * \param slub SLUB. The slub that defines the part of the objective space to search in.
	 * \return the index of the most often fractional variable, as an int.
	 */
	int computeMostOftenFractionalIndex(SLUB& slub);

	/*! \brief Computes the index of the most often fractional variable among the extreme points.
	 *
	 * This functions returns the index of the most often fractional variable in the lower bound set located in the area
	 * defined by the super local upper bound slub. This is done by checking each extreme point that dominates slub.
	 * If there is no fractional variable, it returns the index of the most fractional in average (i.e. avg value closest
	 * to 0.5). If there is no point in the given area, it returns the index of the first free variable it finds.
	 * \param slub SLUB. The slub that defines the part of the objective space to search in.
	 * \param branching decision. Used to make sure the choice of index is not redundant after probing was performed.
	 * \return the index of the most often fractional variable, as an int.
	 */
	int computeMostOftenFractionalIndex(SLUB& slub, BranchingDecisions* newbd);

	/*! \brief Compute the index that violates the most often the branching decisions among the extreme points of the lower bound set.
	 *
	 * \param BranchingDecisions* newbd. A pointer to the branching decisions used for violation test.
	 * \return the index of the variable chosen, as an int.
	 */
	int computeMostViolatingIndex(BranchingDecisions* newbd);

	/*! \brief Computes the median value of variable $x_i$ among the extreme points of the linear relaxation.
	 *
	 * This functions returns the median value of the variable $x_i$ among the extreme points of the LP relax.
	 * \param slub SLUB. The slub that defines the part of the objective space to search in.
	 * \param i int. Index of the splitting variable.
	 * \return the index of the most often fractional variable, as an int.
	 */
	int computeMedianSplittingValue(SLUB& slub, int i, int maxUb);

	/*! \brief Computes the median value of variable $x_i$ among the extreme points of the linear relaxation.
	 *
	 * This functions returns the median value of the variable $x_i$ among the extreme points of the LP relax. The integer values
	 * are ignored in the first place. If there is no decimal values, the integer ones are then considered.
	 * \param slub SLUB. The slub that defines the part of the objective space to search in.
	 * \param i int. Index of the splitting variable.
	 * \return the index of the most often fractional variable, as an int.
	 */
	int computeMedianSplittingValue2(SLUB& slub, int i);

	/*! \brief Computes most often fractional value of variable $x_i$ among the extreme points of the linear relaxation.
	 *
	 * This functions returns the most often fractional value of the variable $x_i$ among the extreme points of the LP relax. If there is no fractional
	 * value, the median splitting rule is called.
	 * \param slub SLUB. The slub that defines the part of the objective space to search in.
	 * \param i int. Index of the splitting variable.
	 * \return the index of the most often fractional variable, as an int.
	 */
	int computeMostOftenFractionalSplittingValue(SLUB& slub, int i, int maxUb);

	/*! \brief Computes most often fractional value of variable $x_i$ among the extreme points of the linear relaxation.
	 *
	 * This functions returns the most often fractional value of the variable $x_i$ among the extreme points of the LP relax. If there is no fractional
	 * value, the median splitting rule is called.
	 * \param slub SLUB. The slub that defines the part of the objective space to search in.
	 * \param i int. Index of the splitting variable.
	 * \return the index of the most often fractional variable, as an int.
	 */
	int computeRandomSplittingValue(SLUB& slub, int i, int maxUb);

	/*! \brief Check whether a point y is weakly dominated by an existing extreme point
	 *
	 * \param y Point. A pointer to the point to look at.
	 * \return true if there exists an extreme point that weakly dominates y.
	 */
	bool isWeeklyDominated(Point* y);

	/*! \brief Compute the extreme rays derived from hyperplane H and attached to vertex y.
	*
	* \param vtx Point*. A pointer to the vertex.
	* \param H Hyperplane*. A pointer to the hyperplane potentially generating extreme rays.
	* \param allNewRays. A pointer to the vector that will contains all the new rays.
	*/
	void computeExtremeRays(Point* vtx, Hyperplane* H, std::vector<Point*>* allNewRays);

	/*! \brief Check whether this extreme ray is valid.
	*
	* Search among the other new rays whether there is a ray such that the set of active hyperplane of this point is included in
	* the set of active hyperplanes of the other point.
	* \param vtx Point*. A pointer to the vertex.
	* \param allNewRays. A pointer to the vector that will contains all the new rays.
	* \return true if it is not discarded.
	*/
	bool filterExtremeRays(Point* vtx, std::vector<Point*>* allNewRays);

	/* Generate the 2^p points of the box defined by antiIdealPoint and yI.
	 *
	 * Done in a recursive way.
	 * \param yI. The ideal point, used for the definition of the box.
	 * \param coord. Coordinates of the the point being built.
	 * \param int index. Used for the recursion.
	 * \param H. Vector of pointers to the hyperplanes that defines the ideal point yI.
	 */
	void generateBox(std::vector<double> yI, std::vector<double> coord, int index, std::vector<Hyperplane*> H);

	/*! \brief Write into stat the statistics about this LP-relax
	 *
	 * This function update the lpStat field of stat.
	 * \param statBB Statistics*. A pointer to the data structure that store statistics for the whole Branch and Bound
	 */
	//void getStatistics(Statistics* statBB);

	/*! \brief Return the pointer to the weighted sum model.
	 *
	 * This function returns a pointer to the weighted sum model.
	 * \return a pointer to a WeightedSumModel.
	 */
	WeightedSumModel* getWeightedSumModel();

	/*! \brief Return the pointer to the feasibility check model.
	 *
	 * This function returns a pointer to the feasibility check model.
	 * \return a pointer to a FeasibilityCheckModel.
	 */
	FeasibilityCheckModel* getFeasibilityCheckModel();

	/*! \brief Return the pointer to the dual benson model.
	 *
	 * This function returns a pointer to the dual benson model.
	 * \return a pointer to a DualBensonModel.
	 */
	DualBensonModel* getDualBensonModel();

	/*! \brief Return the pointer to the furthest feasible point model.
	 *
	 * This function returns a pointer to the further feasible point model.
	 * \return a pointer to a FurthestFeasiblePointModel.
	 */
	FurthestFeasiblePointModel* getFurthestFeasiblePointModel();

	/*! \brief Return true if the point should be checked.
	 *
	 * This function checks whether (1) the point has already been checked by Cplex and (2) if not, if it is weakly dominated by
	 * another existing point. If none of these conditions are met, the point should be checked for feasibility and the function
	 * returns true.
	 * \param y Point*. A pointer to the point being tested.
	 * \return true if the point should be checked.
	 */
	bool shouldBeChecked(Point* y);

	/*! \brief Return true if the point should be checked.
	 *
	 * This function checks whether (1) the point has already been checked by Cplex and (2) if not, if it is weakly dominated by
	 * another existing point. If none of these conditions are met, the point should be checked for feasibility and the function
	 * returns true.
	 * \param y Point*. A pointer to the point being tested.
	 * \return true if the point should be checked.
	 */
	bool shouldBeCheckedFull(Point* y);

	/*! \brief Clear the redundant facets.
	 */
	void clearRedundantFacets();

	/*! \brief Clear the dsicarded points.
	 *
	 * Note: new degenerate rays may be discovered
	 * \param D, list of Point*. A list containing the degenerate vertices.
	 */
	void clearDiscardedPoints();

	/*! \brief Updates the adjacency lists of the points contained in lists N and D
	 *
	 * Done using the double description method (vertex enumeration).
	 * \param N, list of Point*. The list of the polyhedron's new vertices.
	 * \param D, list of Point*. The list of the polyhedron's degenerate vertices.
	 */
	void updateAjacencyLists(std::list<Point*>& N, std::list<Point*>& D);

	/*! \brief Clear the status of all the concerned points before the next iteration.
	 *
	 * Note: D (degenerate points) is included in M (modified points).
	 * \param N, list of Point*. The list of the polyhedron's new vertices.
	 * \param M, list of Point*. The list of the polyhedron's degenerate vertices.
	 * \param R, list of Point*. The list of the polyhedron's new extreme rays
	 */
	void clearStatus(std::list<Point*>& N, std::list<Point*>& M, std::list<Point*>& R);

	/*! \brief Clear the status of all the concerned points before the next iteration.
	 *
	 * Note: D (degenerate points) is included in M (modified points).
	 * \param N, list of Point*. The list of the polyhedron's new vertices.
	 * \param M, list of Point*. The list of the polyhedron's degenerate vertices.
	 * \param R, list of Point*. The list of the polyhedron's new extreme rays
	 */
	void correctRedundancy(std::list<Point*>& N, std::list<Point*>& D);

	/*! \brief Generate an extreme ray in direction k at point w and filter redundant rays in R.
	 *
	 * \param w Point*. A pointer to the point attached to the ray computed.
	 * \param k int. The unbounded direction of the ray.
	 * \param R, list of Point*. The list of the polyhedron's new rays.
	 */
	void generateRay(Point* w, int k, std::list<Point*>& R);

	/* \brief Clear all the cuts generated in the model used to comptue the linear relaxation.
	 */
	void clearCuts();

	/*! \brief Add the cover cut described in cc to the model.
	 *
	 * \param cc vector of vector of int. Represent the cover cuts. Each row is a new cut, column 0 is the rhs, and the other columns are the indices of the variables to add to the cut.
	 */
	void applyCoverCuts(std::vector<std::vector<int>>& cc);

	/*! \brief Class the indices of the variables from the most to the least fractional among the extreme points of the lower bound set.
	 *
	 * \param cc vector int. Indices are stored in this vector.
	 */
	void computeFractionalProportion(std::vector<int>* I, SLUB& s);

	/*! \brief Check for each free variable wether each extreme value occurs among extreme points of the lower bound set in the region described in bd.
	 *
	 * \param bd BranchingDecisions*. A pointer to the branching decisions used for reference.
	 * \param h vector of vector of bool. If a variable take an extreme value, its corresponding cell is set to false.
	 */
	void checkFeasibleValues(BranchingDecisions* bd, std::vector<std::vector<bool>>* h);

	/* \brief Check whether new values for fixing variables should be tested afted variable i was fixed to v.
	 *
	 * \param toTest, vector of vector of bool. Where the values to check for each variables are recorded.
	 * \param int i. The index of the variable that was fixed.
	 * \param int v. The value to which x_i was fixed.
     */
	void updateFeasibleValues(std::vector<std::vector<bool>>* toTest, BranchingDecisions* bd, int i, int v);

	/* Get the minimal difference between the rhs of an hyperplane and the w.s. of local upper bound u using the normal vector of the hyperplane
	 * as the weigth vector.
	 *
	 * \param LocalUpperBound u. The local upper bound used for reference.
	 */
	double getMinFacetGap(LocalUpperBound& u);

	/* Retrieve the w.s. value of the extreme point with the smallest w.s. value, using weights l.
	 *
	 * \param vector of double, l. The weight vector.
	 */
	double getSmallestWeightedSumValue(std::vector<double>& l);

	/* Retrieve the variables that are not fixed and that takes at least once a fractional value among the extreme points of the LB set.
	 * If none takes a fractional value, returns the variables that take at least once different values.
	 * If no extreme point exists in the area defined by the subproblem, the set of all free variables is returned.
	 *
	 * \param isCandidate, vector of bool. True if variable i is candidate, false otherwise.
	 * \param BranchingDecision* bd. The branching decisions considered.
	 */
	void getCandidateSet(std::vector<bool>& isCandidate, BranchingDecisions* bd);

	/* Retrieve the set of candidates according the fractionnality of the extreme point of the LB set.
	 * Get the set of variables that is the most often fractional.
	 * If none is fractional, returns the set of variables taking different values.
	 * If no solution exist in the region, return the set of free variables.
	 */
	void getFractionalCandidateSet(std::vector<bool>& isCandidate, BranchingDecisions* bd);

	void checkViolatingCut();
	void checkViolatingCut(BranchingDecisions* bd);

	bool filterNdLubs(std::list<LocalUpperBound*>& lubDomi, Hyperplane* h);

	bool takeValue(BranchingDecisions* bd, int i, int v);

	std::vector<double> getAverageNormalVector(BranchingDecisions* bd);

	double getPercentageIntegrality();
	double getPercentageFeasibility();

	std::vector<double> getNadirPoint();

	void clear();

	void exportIteration();

	void exportProblem();

	void debug__resetAdjLists();

	int get_numberFacets();
	int get_numberVertices();
};