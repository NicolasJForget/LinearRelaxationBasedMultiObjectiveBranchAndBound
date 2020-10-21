#include "Point.h"


/* ==========================================================
		Constructors
 ========================================================= */

 /*! \brief Constructor with a known objective vector.

  * \param z std::vector of double. It defines the objective vector of the new point.
  */
Point::Point(std::vector<double> const& z) : objVector(z), preImage(0), adjList(0), activeHyperplanes(0), discarded(false), onBoundingBox(false), feasible(false), degenerate(false), copy(NULL) {};

/*! \brief Constructor.
*/
Point::Point() : objVector(0), preImage(0), adjList(0), activeHyperplanes(0), discarded(false), onBoundingBox(false), feasible(false), degenerate(false), copy(NULL) {};

/* ==========================================================
		Regular methods
 ========================================================= */

 /*! \brief Computes the intersection of the set of active hyperplanes of two points and return its size.
  *
  * \param u Point. The second point used for the comparison.
  * \return the number of common active hyperplanes
  */
int Point::sizeIntersectionActiveHyperplanes(Point& u) {

	int s = 0;
	std::list<Hyperplane*>::iterator itThis;
	std::list<Hyperplane*>::iterator itu;
	std::list<Hyperplane*>* actHu = u.get_activeHyperplanes(); // pointer to list of active Hyperplane of u

	for (itThis = activeHyperplanes.begin(); itThis != activeHyperplanes.end(); ++itThis) { // for each active hyperplane of u
		for (itu = actHu->begin(); itu != actHu->end(); ++itu) { // for each active hyperplane of this point
			if (*itThis == *itu) { // if the two hyperplanes are the same
				++s;
			}
		}
	}

	return s;
}

/*! \brief Save the pre-image found for this point.
 *
 * \param M CplexModel*. A pointer to the Cplex model that check for the feasibility of a point.
 */
void Point::isNowFeasible(FeasibilityCheckModel* M) {

	feasible = true;
	preImage.resize(M->get_n());
	M->retrieveSolutionFeasibility(preImage);

}

/*! \brief Check whether the point is located above an hyperplane.
 *
 * \param H Hyperplane.
 * \return true, if the point is located above H.
 */
bool Point::above(Hyperplane& H) {

	double val = 0;
	for (int k = 0; k < H.get_dim() + 1; k++) {
		val += H.get_normalVector(k) * objVector[k];
	}

	return val >= H.get_rhs();

}

/*! \brief Check whether the point is located below an hyperplane.
 *
 * \param H Hyperplane.
 * \return true, if the point is located below H.
 */
bool Point::below(Hyperplane& H) {

	double val = 0;
	for (int k = 0; k < H.get_dim() + 1; k++) {
		val += H.get_normalVector(k) * objVector[k];
	}

	return val <= H.get_rhs();

}

/*! \brief Check whether the point is located on an hyperplane.
 *
 * \param H Hyperplane.
 * \return true, if the point is located on H.
 */
bool Point::locatedOn(Hyperplane& H) {

	double val = 0;
	for (int k = 0; k < H.get_dim() + 1; k++) {
		val += H.get_normalVector(k) * objVector[k];
	}

	return abs(val - H.get_rhs()) <= 0.0000000001; // the epsilon value

}

/*! \brief Compute the intersection point between an edge defined by this point, and an hyperplane.
 *
 * \param u Point. The second point that defines the edge used for the computation.
 * \param H Hyperplane. The hyperplane used for the computation.
 * \return a std::vector of double that defines the coordinates of the intersection point.
 */
std::vector<double> Point::edgeIntersection(Point& u, Hyperplane& H) {

	// compute the value of lambda, to find the point on the edge
	double denominator = 0;
	double lambda = H.get_rhs();
	for (int l = 0; l < H.get_dim() + 1; l++) {
		lambda -= H.get_normalVector(l) * objVector[l];
		denominator += H.get_normalVector(l) * (u.get_objVector(l) - objVector[l]);
	}
	lambda /= denominator;

	// compute the actual point
	std::vector<double> intersection(H.get_dim() + 1);
	for (int k = 0; k < H.get_dim() + 1; k++) {
		intersection[k] = lambda * u.get_objVector(k) + (1 - lambda) * objVector[k];
	}

	return intersection;

}

/*! \brief Find the missing coordinate of this point, given that it is located on a specific hyperplane.
 *
 * \param H Hyperplane. The hyperplane used for the computation.
 * \param coord integer. The index of the missing coordinate.
 * \return the value of this coordinate after computation, as a double.
 */
double Point::findMissingCoordinate(Hyperplane& H, int coord) {

	double val = H.get_rhs();
	for (int k = 0; k < H.get_dim() + 1; k++) {
		if (k != coord) {
			val -= H.get_normalVector(k) * objVector[k];
		}
	}
	val /= H.get_normalVector(coord);

	return val;

}

/*! \brief Change the value of the objective vector at a specific coordinate.
 *
 * \param obj integer. The index of the coordinate to change.
 * \param val double. The new value of the coordinate.
 */
void Point::setObjVector(int obj, double val) {
	objVector[obj] = val;
}

/*! \brief Add a new point to the adjacency list.
 *
 * \param adj Point*. A pointer to the new adjacent point.
 */
void Point::addAdjacentPoint(Point* adj) {
	adjList.push_back(adj);
}

/*! \brief Add a new active hyperplane.
 *
 * \param id Hyperplane*. A pointer to the new active hyperplane.
 */
void Point::addActiveHyperplane(Hyperplane* H) {
	activeHyperplanes.push_back(H);
}

/*! \brief Set this point as discarded.
 */
void Point::becomesDiscarded() {
	discarded = true;
}

/*! \brief Set this point as degenerated.
 */
void Point::becomesDegenerate() {
	degenerate = true;
}

/*! \brief Set this point as non-degenerated.
 */
void Point::becomesNonDegenerate() {
	degenerate = false;
}

/*! \brief Set this point as "on the bounding box".
 */
void Point::becomesOnBoundingBox() {
	onBoundingBox = true;
}

/*! \brief Replace an adjacent vertex by a new one.
 *
 * \param oldVertex Point*. A pointer to the vertex to remove from the list.
 * \param newVertex Point*. A pointer to the new vertex to add to the list.
 */
void Point::replaceAdjVertex(Point* oldVertex, Point* newVertex) {
	adjList.remove(oldVertex);
	adjList.push_back(newVertex);
}

/*! \brief Update the set of active hyperplanes of this point.
 *
 * \param u Point&. The first point that defines the corresponding edge. Passed by reference.
 * \param v Point&. The second point that defines the corresponding edge. Passed by reference.
 * \param H Hyperplane*. A pointer to the corresponding hyperplane.
 */
void Point::updateActiveHyperplanes(Point& u, Point& v, Hyperplane* H) {

	std::list<Hyperplane*>* hu = u.get_activeHyperplanes();
	std::list<Hyperplane*>::iterator itu;
	std::list<Hyperplane*>* hv = v.get_activeHyperplanes();
	std::list<Hyperplane*>::iterator itv;

	for (itu = hu->begin(); itu != hu->end(); ++itu) { // for each active hyperplane of u
		for (itv = hv->begin(); itv != hv->end(); ++itv) { // for each active hyperplane of v
			if (*itu == *itv) { // if there is a common hyperplane
				activeHyperplanes.push_back(*itu);
			}
		}
	}
	activeHyperplanes.push_back(H); // add the new hyperplane H as well
}

/*! \brief Check whether the point has a known feasible pre-image.
 *
 * \return true, if the point has a known pre-image.
 */
bool Point::isFeasible() {
	return feasible;
}

/*! \brief Check whether the point is on the bounding box.
 *
 * \return true, if the point is on the bounding box.
 */
bool Point::isOnBoundingBox() {
	return onBoundingBox;
}

/*! \brief Print this point.
 */
void Point::print() {
	int s = objVector.size();
	std::cout << " ( ";
	for (int k = 0; k < s - 1; k++) {
		std::cout << objVector[k] << " , ";
	}
	std::cout << objVector[s - 1] << " ) ";
}

/*! \brief Check whether the point is already discarded.
 *
 * \return true, if the point is discarded
 */
bool Point::isDiscarded() {
	return discarded;
}

/*! \brief Check whether the point is degenerated.
 *
 * \return true, if the point is degenerated
 */
bool Point::isDegenerate() {
	return degenerate;
}

/*! \brief Register the address of the last copy of this Point.
 *
 * \param P Point*. The address of the copy.
 */
void Point::setCopy(Point* P) {
    copy = P;
}

 /* ==========================================================
		 Getters
  ========================================================= */

/*! \brief Returns a pointer to the objective vector.
 *
 * \return a pointer to the objective vector of this point.
 */
double Point::get_objVector(int obj) {
	return objVector[obj];
}

/*! \brief Returns the value of a specific variable.
 *
 * \param var integer. The index of the variable to look at.
 * \return the value of this variable, as a double.
 */
double Point::get_preImage(int var) {
	return preImage[var];
}

/*! \brief Returns the value of a specific objective function.
 *
 * \param obj integer. The index of the objective to look at.
 * \return the value of this objective, as a double.
 */
std::vector<double>* Point::get_objVector() {
	return &objVector;
}

/*! \brief Returns a pointer to list of adjecent points.
 *
 * \return a pointer to the list of adjacent points.
 */
std::list<Point*>* Point::get_adjList() {
	return &adjList;
}

/*! \brief Returns a pointer to list of active hyperplanes.
 *
 * \return a pointer to the list of adjacent points.
 */
std::list<Hyperplane*>* Point::get_activeHyperplanes() {
	return &activeHyperplanes;
}

/*! \brief Returns a pointer to this point.
 *
 * \return a pointer to this point.
 */
Point* Point::get_adress() {
	return this;
}

/*! \brief Returns the address of the last copy of this point
 *
 * \return a pointer to the copy.
 */
Point* Point::get_copy() {
    return copy;
}
