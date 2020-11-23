/**
* \author Nicolas Forget
*
* This class describe an extreme point of the lower bound set. Used in the linear relaxation (LP relax).
*/

#pragma once
#include <stdlib.h>
//#include "LinearProgram.h"
#include "Model.h"
#include "Hyperplane.h"
#include <vector>
#include <list>
#include <math.h>

class Point
{
private:
	std::vector<double> objVector; //!< objective vector of the point - coordinates in the objective space
	std::vector<double> preImage; //!< pre-image of the point - coordinates in the solution space. Note that if the point is on the bounding box, the pre-image is actually the pre-image of a feasible point that weakly dominates it.
	std::list<Point*> adjList; //!< pointers to the points adjacent to this point in polyhedron of the LP relax
	std::list<Hyperplane*> activeHyperplanes; //!< pointers to the hyperplanes active at this point - i.e. this point is located on these hyperplanes
	bool discarded; //!< true if the point becomes discarded by a new hyperplane added to the LB set.
	bool onBoundingBox; //!< true if it's located on the bounding box of the LP relax - i.e. share at least 1 component with the anti-ideal of the LP relax.
	bool feasible; //!< true if it has a feasible pre-image (for the LP relax)
	bool degenerate; //!< true if the point is degenerated, i.e. located on a new hyperplane added to the LB set.
	bool integer; //!< true if the point has only integer coordinates in the variable space.
	bool integratedInUB; //!< true if the point has already been integrated in the upper bound set.
	Point* copy; //!< used for copy purpose only. See copy constructor of LinearRelaxation for its purpose.

public:

	/*! \brief Constructor with a known objective vector.
	 *
	 * This function creates a new point in the objective space at the location specified by the parameter z.
	 * \param z std::vector of double. It defines the objective vector of the new point.
	 */
	Point(std::vector<double> const& z);

	/*! \brief Constructor,
	 *
	 * This function creates a new point. Nothing is known about pre-image or objective vector yet.
	 */
	Point();

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

	/*! \brief Compute the intersection point between an edge defined by this point, and an hyperplane.
	 *
	 * This function computes the intersection point between an edge, defined by this point and another point u, and an
	 * hyperplane H.
	 * \param u Point. The second point that defines the edge used for the computation.
	 * \param H Hyperplane. The hyperplane used for the computation.
	 * \return a std::vector of double that defines the coordinates of the intersection point.
	 */
	std::vector<double> edgeIntersection(Point& u, Hyperplane& H);

	/*! \brief Find the missing coordinate of this point, given that it is located on a specific hyperplane.
	 *
	 * This function computes the coordinate of this point in the objective space at index coord, given that all
	 * the other coordinates are known, and that it is located on the hyperplane H.
	 * \param H Hyperplane. The hyperplane used for the computation.
	 * \param coord integer. The index of the missing coordinate.
	 * \return the value of this coordinate after computation, as a double.
	 */
	double findMissingCoordinate(Hyperplane& H, int coord);

	/*! \brief Change the value of the objective vector at a specific coordinate.
	 *
	 * This function set the value of coordinate obj to value val.
	 * \param obj integer. The index of the coordinate to change.
	 * \param val double. The new value of the coordinate.
	 */
	void setObjVector(int obj, double val);

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

	/*! \brief Add a new active hyperplane.
	 *
	 * This function adds a new hyperplane to the list of active hyperplanes.
	 * \param id Hyperplane*. A pointer to the new active hyperplane.
	 */
	void addActiveHyperplane(Hyperplane* H);

	/*! \brief Set this point as discarded.
	 *
	 * This function set this point as discarded.
	 */
	void becomesDiscarded();

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

	/*! \brief Set this point as "on the bounding box".
	 *
	 * This function set this point as "on the bouding box", i.e. that it shares at least one component with the 
	 * anti-ideal point of the lower bound set.
	 */
	void becomesOnBoundingBox();

	/*! \brief Replace an adjacent vertex by a new one.
	 *
	 * This function removes an element (given by oldVertex) of the list of adjecent vertex and add a new one (newVertex)
	 * at the end of this list.
	 * \param oldVertex Point*. A pointer to the vertex to remove from the list.
	 * \param newVertex Point*. A pointer to the new vertex to add to the list.
	 */
	void replaceAdjVertex(Point* oldVertex, Point* newVertex); // can be opti in cpp if we can change value of an element in a list

	/*! \brief Update the set of active hyperplanes of this point.
	 *
	 * This function updates the set of active hyperplane of this point, given that it is located on the edge defined
	 * by the points u and v, and that it is located on the hyperplane H.
	 * \param u Point&. The first point that defines the corresponding edge. Passed by reference.
	 * \param v Point&. The second point that defines the corresponding edge. Passed by reference.
	 * \param H Hyperplane*. A pointer to the corresponding hyperplane.
	 */
	void updateActiveHyperplanes(Point& u, Point& v, Hyperplane* H);

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

	/*! \brief Print this point.
	 *
	 * This function prints the point by showing its objective vector in a easily-readible way for the user.
	 */
	void print();

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
};
