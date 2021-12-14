/**
* \author Nicolas Forget
*
* This class describe an extreme point of the lower bound set. Used in the linear relaxation (LP relax).
*/

#pragma once
#include <stdlib.h>
//#include "LinearProgram.h"
#include "Model.h"
#include <vector>
#include <list>
#include <math.h>

class Hyperplane;
class FeasibilityCheckModel;

class Point
{
private:
	std::vector<double> objVector; //!< objective vector of the point - coordinates in the objective space
	std::vector<double> preImage; //!< pre-image of the point - coordinates in the solution space. Note that if the point is on the bounding box, the pre-image is actually the pre-image of a feasible point that weakly dominates it.
	std::list<Point*> adjList; //!< pointers to the points adjacent to this point in polyhedron of the LP relax
	std::list<Hyperplane*> activeHyperplanes; //!< pointers to the hyperplanes active at this point - i.e. this point is located on these hyperplanes
	bool discarded; //!< true if the point becomes discarded by a new hyperplane added to the LB set.
	bool onBoundingBox; //!< true if it's located on the bounding box of the LP relax - i.e. share at least 1 component with the anti-ideal of the LP relax.
	bool feasible; //!< true if it has a feasible pre-image in the previous node (for the father LP relax)
	bool degenerate; //!< true if the point is degenerated, i.e. located on a new hyperplane added to the LB set.
	bool integer; //!< true if the point has only integer coordinates in the variable space.
	bool integratedInUB; //!< true if the point has already been integrated in the upper bound set. -> obsolete ??
	bool newPoint; //!< true if the point has been created in the current node, false if created in a father node
	bool checkpoint; //!< true if it has a feasible known pre-image in the current node
	bool isCplexChecked; //!< true if this point has already been checked by cplex in the current lp relax computation
	Point* copy; //!< used for copy purpose only. See copy constructor of LinearRelaxation for its purpose.
	bool isRay; //!< true if it is an extreme ray
	int unboundedDimension; //!< index of the unbounded dimension. -1 if the point is not a ray.
	//bool hasRay; //!< true if it is adjacent to an extreme ray
	std::vector<bool> hasRay; //!< one component for each dimension, component k is true if it is adjacent to an extreme ray in direction k.
	bool visited; //!< true if the point has been visited at a given iteration of benson's algorithm
	bool modified; //!< true if the point has been modified or used for computing a new point at a given of benson's algorithm
	std::list<Point*> C; //!< list of points related to the modifications (linked to modified)
	int domiSlub; //!< true if it is dominated by the current slub. -1 if unknown, 0 if non-dominated, 1 if dominated
	int refVarFix; //!< true if the point is a reference point for variable fixing, i.e. feasible for the current branching decisions.
	//int directionRay; //!< unbounded direction if the point is a ray

public:

	/*! \brief Constructor with a known objective vector.
	 *
	 * This function creates a new point in the objective space at the location specified by the parameter z.
	 * \param z std::vector of double. It defines the objective vector of the new point.
	 */
	Point(std::vector<double> const& z);

	/*! \brief Constructor of a ray with a known objective vector.
	 *
	 * This function creates a new point in the objective space at the location specified by the parameter z. The presence of the second
	 * parameter implies that the point is a ray. This second parameter describes which dimension is unbounded.
	 * \param z std::vector of double. It defines the objective vector of the new point.
	 * \param k int. Direction of the unbounded direction
	 */
	Point(std::vector<double> const& z, int k);

	/*! \brief Construct the point located at the intersection of the edge defined by u and v and hyperplane H.
	 *
	 * Note: lambda, the parameter that gives the actual location of the edge, is already computed previously to check redundancy.
	 * \param u Point*. A pointer to the first point defining the edge.
	 * \param v Point*. A pointer to the second point defining the edge.
	 * \param lambda double. Defines the location of the intesection on the edge.
	 * \param H Hyperplane*. A pointer to the hyperplane used to compute the intersection.
	 */
	Point(Point* u, Point* v, double lambda, Hyperplane* H);

	/*! \brief Construct the point located at the intersection of the edge defined by u and v and hyperplane H.
	 *
	 * Note: lambda, the parameter that gives the actual location of the edge, is already computed previously to check redundancy.
	 * \param u Point*. A pointer to the first point defining the edge.
	 * \param v Point*. A pointer to the second point defining the edge.
	 * \param lambda double. Defines the location of the intesection on the edge.
	 * \param H Hyperplane*. A pointer to the hyperplane used to compute the intersection.
	 */
	Point(Point* u, Point* v, std::vector<double> coord, Hyperplane* H);

	/*! \brief Constructor,
	 *
	 * This function creates a new point. Nothing is known about pre-image or objective vector yet.
	 */
	Point();

	bool operator==(const Point& b);

	/*! \brief Computes the intersection of the set of active hyperplanes of two points and return its size.
	 *
	 * This function calculates how many common active hyperplanes this point have with another point u.
	 * \param u Point. The second point used for the comparison.
	 * \return the number of common active hyperplanes
	 */
	int sizeIntersectionActiveHyperplanes(Point& u);

	/*! \brief Returns the value of a specific objective function.
	 *
	 * This function returns the value of the objective vector at coordinate obj. This correspond to the value of
	 * this point in objective obj.
	 * \param obj integer. The index of the objective to look at.
	 * \return the value of this objective, as a double.
	 */
	double get_objVector(int obj);

	/*! \brief Returns true if the point is on an already existing edge
	 *
	 * This is mostly a corrective patch to avoid linking non-adjacent points that passes through the algorithm.
	 * (Delete the point on the middle of the edge if such a case is discovered ??)
	 * \param y Point*. One of the three points used for edge comparisons
	 * \return true if it lies on an existing edge defined by two other points
	 */
	int isOnEdge(Point* y);

	/*! \brief Returns the value of a specific variable.
	 *
	 * This function returns the value of the solution vector (pre-image) at coordinate var. This correspond to the
	 * value of this point in variable var.
	 * \param var integer. The index of the variable to look at.
	 * \return the value of this variable, as a double.
	 */
	double get_preImage(int var);

	/*! \brief Returns a pointer to the objective vector.
	 *
	 * It is sometimes required to have the full objective vector at once. This function returns a pointer to this
	 * objective vector.
	 * \return a pointer to the objective vector of this point.
	 */
	std::vector<double>* get_objVector();

	/*! \brief Check whether the point is already discarded.
	 *
	 * This function checks whether the point is discarded.
	 * \return true, if the point is discarded
	 */
	bool isDiscarded();

	/*! \brief Check whether the point is degenerated.
	 *
	 * This function checks whether the point is degenerated.
	 * \return true, if the point is degenerated
	 */
	bool isDegenerate();

	/*! \brief Check whether the point is new, i.e. created in the current node.
	 *
	 * This function checks whether the point is new.
	 * \return true, if the point is new.
	 */
	bool isNew();

	/*! \brief Check whether the point is feasible in the current node.
	 *
	 * This function checks whether the point is feasible in the current node, by reading checkpoint.
	 * \return true, if the point is feasible.
	 */
	bool isCheckpoint();

	/*! \brief Returns a pointer to list of adjecent points.
	 *
	 * This function returns a pointer to the list of points adjecent to this point. This is used to go thought these
	 * adjecent points.
	 * \return a pointer to the list of adjacent points.
	 */
	std::list<Point*>* get_adjList();

	/*! \brief Returns a pointer to list of active hyperplanes.
	 *
	 * This function returns a pointer to the list of active hyperplanes. This is used to go thought these active
	 * hyperplanes.
	 * \return a pointer to the list of adjacent points.
	 */
	std::list<Hyperplane*>* get_activeHyperplanes();

	/*! \brief Returns a pointer to this point.
	 *
	 * This function returns a pointer to this point.
	 * \return a pointer to this point.
	 */
	Point* get_adress();

	/*! \brief Check whether the point is located above an hyperplane.
	 *
	 * This function checks if this point is located above the hyperplane H, by using its defining inequality.
	 * Note that this is not strictly above, but above or on. Hence, it returns true if $n^H.y \geq d$, where $n^H$
	 * is the normal vector of H, $y$ is this point, and $d$ is the right-hand side of H.
	 * \param H Hyperplane.
	 * \return true, if the point is located above H.
	 */
	bool above(Hyperplane& H);

	/*! \brief Check whether the point is located below an hyperplane.
	 *
	 * This function checks if this point is located below the hyperplane H, by using its defining inequality.
	 * Note that this is not strictly below, but above or on. Hence, it returns true if $n^H.y \leq d$, where $n^H$
	 * is the normal vector of H, $y$ is this point, and $d$ is the right-hand side of H.
	 * \param H Hyperplane.
	 * \return true, if the point is located below H.
	 */
	bool below(Hyperplane& H);

	/*! \brief Check whether the point is located below an hyperplane.
	 *
	 * This function checks if this point is located below the hyperplane H, by using its defining inequality.
	 * Note that this is not strictly below, but above or on. Hence, it returns true if $n^H.y \leq d$, where $n^H$
	 * is the normal vector of H, $y$ is this point, and $d$ is the right-hand side of H.
	 * CONSIDER EQUALITY AND AN EPSILON ABOVE AS BELOW.
	 * \param H Hyperplane.
	 * \return true, if the point is located below H.
	 */
	bool below2(Hyperplane& H);

	/*! \brief Check whether the point is located below an hyperplane.
	 *
	 * This function checks if this point is located below the hyperplane H, by using its defining inequality.
	 * Note that this is not strictly below, but above or on. Hence, it returns true if $n^H.y \leq d$, where $n^H$
	 * is the normal vector of H, $y$ is this point, and $d$ is the right-hand side of H.
	 * CONSIDER EQUALITY AND AN EPSILON ABOVE AS BELOW.
	 * \param H Hyperplane.
	 * \return true, if the point is located below H.
	 */
	bool isLocatedOn(Hyperplane& H);

	/*! \brief Check whether the point is located on an hyperplane.
	 *
	 * This function checks if this point is located on the hyperplane H, by using its defining inequality. Due to 
	 * floating point rounding, it happens that a point is detected as not on the hyperplane, while it is
	 * actually the case. Hence, a small epsilon is introduced. This function returns true if $|n^H.y - d| \leq \epsilon$
	 * where $n^H$ is the normal vector of H, $y$ is this point, and $d$ is the right-hand side of H. Currently,
	 * $\epsilon = 10^{-10}$.
	 * \param H Hyperplane.
	 * \return true, if the point is located on H.
	 */
	bool locatedOn(Hyperplane& H);

	double getEdgeCrossingValue(Point* u, Hyperplane* H);

	double distToHyperplane(Hyperplane* H);
	double distToHyperplane(Hyperplane* H, std::vector<double>& yAI);

	/*! \brief Compute the intersection point between an edge defined by this point, and an hyperplane.
	 *
	 * This function computes the intersection point between an edge, defined by this point and another point u, and an
	 * hyperplane H.
	 * \param u Point. The second point that defines the edge used for the computation.
	 * \param H Hyperplane. The hyperplane used for the computation.
	 * \return a std::vector of double that defines the coordinates of the intersection point.
	 */
	std::vector<double> edgeIntersection(Point& u, Hyperplane& H, int callDebug);

	/*! \brief Find the missing coordinate of this point, given that it is located on a specific hyperplane.
	 *
	 * This function computes the coordinate of this point in the objective space at index coord, given that all
	 * the other coordinates are known, and that it is located on the hyperplane H.
	 * \param H Hyperplane. The hyperplane used for the computation.
	 * \param coord integer. The index of the missing coordinate.
	 * \return the value of this coordinate after computation, as a double.
	 */
	double findMissingCoordinate(Hyperplane& H, int coord);

	/*! \brief Return the degree of bounding, i.e. the number of component shared with the anti-ideal points.
	 *
	 * \param antiIdealPoint, vector of doubles. The coordinates of the anti-ideal point.
	 * \return the degree of bounding, as an int.
	 */
	int degreeOfBounding(std::vector<double>& antiIdealPoint);

	/*! \brief Return the degree of bounding, i.e. the number of component shared with the anti-ideal points.
	 *
	 * \param antiIdealPoint, vector of doubles. The coordinates of the anti-ideal point.
	 * \return the degree of bounding, as an int.
	 */
	//int prelimAdjacencyCheck(std::vector<double>& antiIdealPoint, Point* y);

	/*! \brief Change the value of the objective vector at a specific coordinate.
	 *
	 * This function set the value of coordinate obj to value val.
	 * \param obj integer. The index of the coordinate to change.
	 * \param val double. The new value of the coordinate.
	 */
	void setObjVector(int obj, double val);

	/*! \brief Change the value of the pre-image vector at a specific coordinate.
	 *
	 * This function set the value of variable i to value val.
	 * \param i integer. The index of the variable to change.
	 * \param val double. The new value of the coordinate.
	 */
	void setPreImage(int i, double val);

	/*! \brief Initialize the pre-image vector (allocate n cells)
	 */
	void initPreImage(int s);

	/*! \brief Add a new point to the adjacency list.
	 *
	 * This function adds a new point to the adjacency list.
	 * \param adj Point*. A pointer to the new adjacent point.
	 */
	void addAdjacentPoint(Point* adj);

	/*! \brief Remove a specific point in the adjacency list.
	 *
	 * This function removes the point in adj in the adjacency list, if it is in the list.
	 * \param adj Point*. A pointer to the adjacent point to remove.
	 */
	void removeAdjacentPoint(Point* adj);

	/*! \brief Remove a specific point in the adjacency list.
	 *
	 * This function removes the point in adj in the adjacency list, if it is in the list.
	 * \param adj iterator of list of point*. The location of the point to delete in the list of adj points.
	 */
	void removeAdjacentPoint(std::list<Point*>::iterator adj);

	/*! \brief Add a new active hyperplane.
	 *
	 * This function adds a new hyperplane to the list of active hyperplanes.
	 * \param id Hyperplane*. A pointer to the new active hyperplane.
	 */
	void addActiveHyperplane(Hyperplane* H);

	/*! \brief Purge a degenerate vertex by deleting itself from the list of adjacent vertex, if that happened to be the case
	 *
	 * It also update the active hyperplanes.
	 * \param id Hyperplane*. A pointer to the new active hyperplane.
	 */
	void purge();

	/*! \brief Set this point as discarded.
	 *
	 * This function set this point as discarded.
	 */
	void becomesDiscarded();

	/*! \brief Set this point as non discarded.
	 *
	 * This function set this point as non discarded.
	 */
	void becomesNonDiscarded();

	/*! \brief Set this point as degenerated.
	 *
	 * This function set this point as degenerated.
	 */
	void becomesDegenerate();

	/*! \brief Set this point as non-degenerated.
	 *
	 * This function set this point as non-degenerated.
	 */
	void becomesNonDegenerate();

	/*! \brief Set this point as non-newPoint, i.e. as old.
	 *
	 * This function set this point as old, by setting newPoint to false.
	 */
	void becomesOld();

	/*! \brief Set this point as new.
	 *
	 * This function set this point as new, by setting newPoint to true.
	 */
	void becomesNew();

	/*! \brief Set this point as "on the bounding box".
	 *
	 * This function set this point as "on the bouding box", i.e. that it shares at least one component with the 
	 * anti-ideal point of the lower bound set.
	 */
	void becomesOnBoundingBox();

	/*! \brief Set this point as checkpoint.
	 *
	 * This function set this point as checkpoint (see LinearRelaxation).
	 */
	void becomesCheckpoint();

	/*! \brief Set this point as non-checkpoint.
	 *
	 * This function set this point as non-checkpoint & non-cplex-checked (see LinearRelaxation).
	 */
	void becomesNonCheckpoint();

	/*! \brief Note this point as checked by cplex.
	 *
	 * This function set this point as checked by cplex (see LinearRelaxation).
	 */
	void becomesCplexChecked();

	void becomesRefVarFix();
	void becomesNonRefVarFix();
	bool isRefVarFix();

	/*! \brief Check whether this point has been checked by cplex in this iteration.
	 *
	 * \return true if this point has been checked by cplex, false otherwise.
	 */
	bool isCheckedByCplex();

	/*! \brief Replace an adjacent vertex by a new one.
	 *
	 * This function removes an element (given by oldVertex) of the list of adjecent vertex and add a new one (newVertex)
	 * at the end of this list.
	 * \param oldVertex Point*. A pointer to the vertex to remove from the list.
	 * \param newVertex Point*. A pointer to the new vertex to add to the list.
	 */
	void replaceAdjVertex(Point* oldVertex, Point* newVertex); // can be opti in cpp if we can change value of an element in a list

	/*! \brief Check whether this point is adjacent to another point y.
	 *
	 * This is done by searching for y in the adjacency list of this point. If y is there, then the two points are adjacent.
	 * \param y Point. Pointer to the point used for comparison.
	 * \return true if the two points are adjacent.
	 */
	bool isAdjacent(Point* y);

	/*! \brief Update the set of active hyperplanes of this point.
	 *
	 * This function updates the set of active hyperplane of this point, given that it is located on the edge defined
	 * by the points u and v, and that it is located on the hyperplane H.
	 * It also notify the common hyperplanes that a new point defines them.
	 * \param u Point&. The first point that defines the corresponding edge. Passed by reference.
	 * \param v Point&. The second point that defines the corresponding edge. Passed by reference.
	 * \param H Hyperplane*. A pointer to the corresponding hyperplane.
	 * \return true if this point is actually an extreme point, false if it lies on a face only or in the interior of a facet.
	 */
	bool updateActiveHyperplanes(Point& u, Point& v, Hyperplane* H, Hyperplane* BB);
	//void updateActiveHyperplanes(Point& u, Point& v, Hyperplane* H);

	/*! \brief Update the set of active hyperplanes of this ray. Check whether it's a redundant ray at the sme time.
	 *
	 * This function updates the set of active hyperplane of a ray attached to point v, located on hyperplane H and unbouded
	 * in dimension unboundedDim.
	 * It also notify the common hyperplanes that a new point defines them.
	 * \param v Point&. The second point that defines the corresponding edge. Passed by reference.
	 * \param H Hyperplane*. A pointer to the corresponding hyperplane.
	 * \param unboundedDim int. The index of the unbounded dimension
	 * \return true if the ray is valid
	 */
	bool updateActiveHyperplanes(Point& v, Hyperplane* H, int unboundedDim);

	/*! \brief Update the set of active hyperplanes of this point.
	 *
	 * This function updates the set of active hyperplane of this point, given that it is located on the edge defined
	 * by the points u and v.
	 * \param u Point&. The first point that defines the corresponding edge. Passed by reference.
	 * \param v Point&. The second point that defines the corresponding edge. Passed by reference.
	 */
	void updateActiveHyperplanes(Point& u, Point& v);

	/*! \brief Check whether the point is feasible given the branching decisions made
	 *
	 * This function looks first at the decision space and check whether the last split variable is within the bounds defined
	 * in the branching decision. Then, it checks whether the point is in the objective branching cone defined in the branching
	 * decisions.
	 * \param BD BranchingDecisions. A pointer to the branching decisions used for comparison.
	 * \return true if the point is feasible given the branching decisions, false otherwise.
	 */
	bool satisfyBranchingDecisions(BranchingDecisions* BD);

	/* Retrieve the w.s. value of the point, with weights l.
	 *
	 * \param vector of double, l. The weight vector.
	 */
	double getWeightedSumValue(std::vector<double>& l);

	bool isNumericallyOnHyperplane(Hyperplane* H);

	bool isNumericallyRayOf(Point* y);

	void connect(Point* y);

	void disconnect(Point* y);

	void checkAdjacencyOld(Point* y, std::vector<Point*>& newVertices, std::vector<Point*>& degenerateVertices);

	void checkAdjacencyOld(Point* y, std::vector<Point*>& newVertices, std::list<Point*>& all);

	void debug__clearAdjList();

	/*! \brief Check whether the point has a known feasible pre-image.
	 *
	 * This function checks whether the point has a known feasible pre-image by looking at feasible.
	 * \return true, if the point has a known pre-image.
	 */
	bool isFeasible();

	/*! \brief Check whether the point is on the bounding box.
	 *
	 * This function checks whether the point is on the bounding box, by looking at onBoundingBox.
	 * \return true, if the point is on the bounding box.
	 */
	bool isOnBoundingBox();

	/*! \brief Check whether the point is on the bounding box.
	 *
	 * This function checks whether the point is on the bounding box, by looking at the common comonents between objVector
	 * and yI. If there is at least one, the point is on the bounding box.
	 * \param yI vector of double. The point defining the bounding box.
	 * \return true, if the point is on the bounding box.
	 */
	bool isOnBoundingBox(std::vector<double>& yI);

	/*! \brief Check whether the point is already in the upper bound set.
	 *
	 * This function checks whether the point is already in the upper bound set, by looking at integreatedInUB.
	 * \return true, if the point is in the upper bound set.
	 */
	bool isInUB();

	/*! \brief Set the point as "integrated in UB"
	 *
	 * This function checks set integratedInUB to true, meaning that the point has been added to the upper bound set if not
	 * dominated.
	 */
	void setAsIntegratedInUB();

	/*! \brief Check whether this point dominates another point y
	 *
	 * \param y Point*. A pointer to the other point used for comparison.
	 * \return true if this point dominates y.
	 */
	bool dominates(Point* y);

	/*! \brief Print this point.
	 *
	 * This function prints the point by showing its objective vector in a easily-readible way for the user.
	 */
	void print();

	/*! \brief Print this point.
	 *
	 * This function prints the point by showing its variable vector in a easily-readible way for the user.
	 */
	void printPreImage();

	/*! \brief Save the pre-image found for this point.
	 *
	 * This function set this point as feasible, and save the pre-image found for this point by asking the solution
	 * to the Cplex model M.
	 * \param M CplexModel*. A pointer to the Cplex model that check for the feasibility of a point.
	 */
	void isNowFeasible(FeasibilityCheckModel* M);

	/*! \brief Register the address of the last copy of this Point.
	 *
	 * This function register the address of the last copy of thus Point. This is done for copying LinearRelaxations
	 * and link the pointers properly.
	 * \param P Point*. The address of the copy.
	 */
	void setCopy(Point* P);

	/*! \brief Returns the address of the last copy of this point
	 *
	 * This function returns a pointer to the last copy of this point.
	 * \return a pointer to the copy.
	 */
	Point* get_copy();

	/*! \brief Returns the number of variables.
	 *
	 * This function returns the number of variables.
	 * \return the number of variables, as an int.
	 */
	int get_nbVar();

	/*! \brief Returns the number of objectives.
	 *
	 * This function returns the number of objectives.
	 * \return the number of objectives, as an int.
	 */
	int get_nbObj();

	/*! \brief Check whether the point is integer in the variable space.
	 *
	 * This function checks whether the point is integer in the variable space, i.e. whether its pre-image has only
	 * integer coordinates.
	 * \return true if the point is integer, false otherwise.
	 */
	bool isInteger();

	/*! \brief Check whether the point is significantly close to y.
	 *
	 * The two points are considered significantly close if their euclidian distance is smaller than a threashold thsld.
	 * This function computes the euclidian distance and compare it to the threashold value.
	 * \param Point* y. The point used for comparison with this.
	 * \return true if the two points are significantly close.
	 */
	bool isVeryCloseTo(Point* y);

	/*! \brief Check whether the point is significantly close to y.
	 *
	 * The two points are considered significantly close if their euclidian distance is smaller than a threashold thsld.
	 * This function computes the euclidian distance and compare it to the threashold value.
	 * \param Point* y. The point used for comparison with this.
	 * \return true if the two points are significantly close.
	 */
	bool isVeryCloseTo(Point* y, int callDebug);

	/*! \brief Check whether the point is significantly close to y.
	 *
	 * The two points are considered significantly close if their euclidian distance is smaller than a threashold thsld.
	 * This function computes the euclidian distance and compare it to the threashold value.
	 * \param Point* y. The point used for comparison with this.
	 * \return true if the two points are significantly close.
	 */
	bool isVeryCloseTo(std::vector<double> y);

	/*! \brief Merge the point with y.
	 *
	 * The functions replaces the coordinates of the current point by the coordinates of the point located at the middle
	 * of the edge defined by this point and y. It also merges the list of active hyperplanes and adjacent points.
	 * \param Point* y. The point used for comparison with this.
	 */
	void merge(Point* y);

	/*! \brief Merge the point with y.
	 *
	 * The functions replaces the coordinates of the current point by the coordinates of the point located at the middle
	 * of the edge defined by this point and y. It also merges the list of active hyperplanes and adjacent points.
	 * \param Point* y. The point used for comparison with this.
	 */
	void merge2(Point* y);

	/*! \brief Discard the hyperplane H from the list of active hyperplanes.
	 *
	 * The functions search for H in the list of active hyperplanes and removes it.
	 * \param Hyperplane* H. The hyperplane to discard.
	 */
	void discardHyperplane(Hyperplane* H);

	/*! \brief Check whether the point is an extreme ray.
	 * \return true if the point is a ray, false otherwise.
	 */
	bool is_ray();

	/*! \brief Label the point as extreme ray.
	 */
	void becomes_ray();
	
	/*! \brief Check whether the point is adjacent to an extreme ray in a given direction k.
	 *
	 * \param k, int. The index of the unbounded dimension being considered.
	 * \return true if the point is adjacent to an extreme ray unbounded in dimension k, false otherwise.
	 */
	bool has_ray(int k);
	
	/*! \brief Notify the point that it is now connected to an extreme ray unbounded in direction k.
	 *
	 * \param k, int. The unbounded direction of the ray.
	 */
	void receive_ray(int k);

	/*! \brief Notify the point that it is not longer connected to its extreme ray unbounded in direction k.
	 *
	 * \param k, int. The unbounded direction of the ray.
	 */
	void loose_ray(int k);

	/*! \brief Notify the adjacent points and active hyperplane that this point will be destroyed.
	 */
	void notifyDeletion();

	/*! \brief Returns the directions of the extreme rays adjacent to this vertex
	 *
	 * \param yAI, vector of double. The anti-ideal point, used to detect the direction of a ray.
	 * \return vector of int, containing the indices of the direction of the rays.
	 */
	std::vector<bool> getRayDirections(std::vector<double> yAI);

	/*! \brief Label the point as modified.
	 */
	void becomesModified();

	/*! \brief Label the point as non-modified.
	 */
	void becomesNonModified();

	/*! \brief Label the point as visited.
	 */
	void becomesVisited();

	/*! \brief Return true if the point is labelled as visited.
	 *
	 * \return true if the point is labelled as visited.
	 */
	bool isVisited();

	/*! \brief Return true if the point is labelled as modified.
	 *
	 * \return true if the point is labelled as modified.
	 */
	bool isModified();

	/*! \brief Label the point as non-visited.
	 */
	void becomesNonVisited();

	/*! \brief Add a point to the list of modifications (C).
	 *
	 * \param u Point*. A pointer to the point involved in the modification.
	 */
	void addModification(Point* u);

	/*! \brief Label the point as degenerate and correct possible computation made considering this point as non-degenerate.
	 *
	 * \param N std::list<Point*>*. The list of points computed so far at this iteration of benson's algorithm.
	 */
	void becomesDegenerate(std::list<Point*>* N, Hyperplane* H);

	/*! \brief Check whether the point is already discarded.
	 *
	 * This function checks whether the point is discarded.
	 * \return true, if the point is discarded
	 */
	bool isDiscarded(Hyperplane* H);

	/*! \brief Check whether the point is already discarded.
	 *
	 * This function checks whether the point is discarded.
	 * \return true, if the point is discarded
	 */
	bool isDegenerate(Hyperplane* H);

	/*! \brief Check adjacency of this point and y.
	 *
	 * This function checks whether this point and y are adjacent. The part of the polyhedron involved are described in N (new vertices)
	 * and D (degenerate vertices).
	 * \param y Point*. A pointer to the point used for adjacency test.
	 * \param N list of Point*. The list of the polyhedron's new vertices.
	 * \param D list of Point*. The list of the polyhedron's degenerate vertices.
	 */
	void checkAdjacency(Point* y, std::list<Point*>& N, std::list<Point*>& D, int iteration);

	/*! \brief Returns true if this point is included in the face defined by the hyperplanes included in F.
	 *
	 * \param F, list of hyperplanes. They describe the face.
	 * \return true if this point is included in the face.
	 */
	bool isIncludedInFace(std::list<Hyperplane*>& F);

	/*! \brief Clears the list of modifications (C)
	 */
	void clearModifications();

	/*! \brief Returns true if the set of active hyperplanes of this point is a subset of the set of active hyperplanes of y.
	 *
	 * If this is true, then this ray is not located on a minimal-dimensional face. In other words, this ray is redundant.
	 * \param y Point*. A pointer to the point used as reference.
	 * \return true if there is inclusion of the active hyperplane.
	 */
	bool hasSameFacets(Point* y);

	/*! \brief Notify the extreme rays attached to this point that they should be deleted in the concerned directions.
	 */
	void notifyRays(Hyperplane* H, std::list<Point*>& M);

	/* \brief Returns the unbounded dimension of the ray.
	 *
	 * \return the dimension, as an int.
	 */
	int getUnboundedDim();

	/* \brief Checks whether the point dominate the current
	 *
	 * -1 : unknown -> should be checked
	 * 0 : non-dominater -> should be skipped
	 * 1 : dominater -> should be considered
	 * \return status of the dominance
	 */
	int isDominater();
	void becomesDominater();
	void becomesNonDominater();
	void becomesUnknownDominater();
};
