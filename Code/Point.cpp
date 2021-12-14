#include "Point.h"
#include "Hyperplane.h"


/* ==========================================================
		Constructors
 ========================================================= */

 /*! \brief Constructor with a known objective vector.

  * \param z std::vector of double. It defines the objective vector of the new point.
  */
Point::Point(std::vector<double> const& z) : objVector(z), preImage(0), adjList(0), activeHyperplanes(0), discarded(false), onBoundingBox(false), feasible(false), degenerate(true), integer(false), integratedInUB(false), newPoint(true), checkpoint(false), isCplexChecked(false), copy(NULL), isRay(false), unboundedDimension(-1), hasRay(z.size(), false), visited(false), modified(false), C(0), domiSlub(-1), refVarFix(false) {}

/*! \brief Constructor with a known objective vector.
 *
 * \param z std::vector of double. It defines the objective vector of the new point.
 * \param ray bool. True if the point is an extreme ray, false otherwise.
 */
Point::Point(std::vector<double> const& z, int k) : objVector(z), preImage(0), adjList(0), activeHyperplanes(0), discarded(false), onBoundingBox(false), feasible(false), degenerate(true), integer(false), integratedInUB(false), newPoint(true), checkpoint(false), isCplexChecked(false), copy(NULL), isRay(true), unboundedDimension(k), hasRay(z.size(), false), visited(false), modified(false), C(0), domiSlub(-1), refVarFix(false) {}

/*! \brief Construct the point located at the intersection of the edge defined by u and v and hyperplane H.
 *
 * Note: lambda, the parameter that gives the actual location of the edge, is already computed previously to check redundancy.
 * Note: no ray is attached to the point yet.
 * \param u Point*. A pointer to the first point defining the edge (the infeasible one).
 * \param v Point*. A pointer to the second point defining the edge (the feasible one).
 * \param lambda double. Defines the location of the intesection on the edge.
 * \param H Hyperplane*. A pointer to the hyperplane used to compute the intersection.
 */
Point::Point(Point* u, Point* v, double lambda, Hyperplane* H) : objVector(u->get_nbObj()), preImage(0), adjList(0), activeHyperplanes(0), discarded(false), onBoundingBox(false), feasible(false), degenerate(true), integer(false), integratedInUB(false), newPoint(true), checkpoint(false), isCplexChecked(false), copy(NULL), isRay(false), unboundedDimension(-1), hasRay(v->get_nbObj(), false), visited(false), modified(false), C(0), domiSlub(-1), refVarFix(false) {

    // compute the coordinates

    for (int k = 0; k < objVector.size(); k++)
        objVector[k] = lambda * u->get_objVector(k) + (1 - lambda) * v->get_objVector(k);
    
    // update the connections of the vertices

    connect(v);
    u->disconnect(v);
    if (v->is_ray()) {
        hasRay[v->getUnboundedDim()] = true;
        //objVector[v->getUnboundedDim()] = v->get_objVector(v->getUnboundedDim());
    }

    // update the modification lists

    C.push_back(v); // we keep the feasible vertex as modification
    u->addModification(this);
    v->addModification(this);

    // update and connect the active hyperplanes

    std::list<Hyperplane*>* hppU = u->get_activeHyperplanes(), * hppV = v->get_activeHyperplanes();
    std::list<Hyperplane*>::iterator f, g;

    for (f = hppU->begin(); f != hppU->end(); f++) {
        for (g = hppV->begin(); g != hppV->end(); g++) {
            if (*f == *g) {
                activeHyperplanes.push_back(*f);
                (*f)->addVertex(this);
                break;
            }
        }
    }
    activeHyperplanes.push_back(H);
    H->addVertex(this);
}

/*! \brief Construct the point located at the intersection of the edge defined by u and v and hyperplane H.
 *
 * Note: lambda, the parameter that gives the actual location of the edge, is already computed previously to check redundancy.
 * Note: no ray is attached to the point yet.
 * \param u Point*. A pointer to the first point defining the edge (the infeasible one).
 * \param v Point*. A pointer to the second point defining the edge (the feasible one).
 * \param lambda double. Defines the location of the intesection on the edge.
 * \param H Hyperplane*. A pointer to the hyperplane used to compute the intersection.
 */
Point::Point(Point* u, Point* v, std::vector<double> coord, Hyperplane* H) : objVector(u->get_nbObj()), preImage(0), adjList(0), activeHyperplanes(0), discarded(false), onBoundingBox(false), feasible(false), degenerate(true), integer(false), integratedInUB(false), newPoint(true), checkpoint(false), isCplexChecked(false), copy(NULL), isRay(false), unboundedDimension(-1), hasRay(v->get_nbObj(), false), visited(false), modified(false), C(0), refVarFix(false) {

    // compute the coordinates

    for (int k = 0; k < objVector.size(); k++)
        objVector[k] = coord[k];

    // update the connections of the vertices

    connect(v);
    u->disconnect(v);
    if (v->is_ray()) {
        hasRay[v->getUnboundedDim()] = true;
        //objVector[v->getUnboundedDim()] = v->get_objVector(v->getUnboundedDim());
    }

    // update the modification lists

    C.push_back(v); // we keep the feasible vertex as modification
    u->addModification(this);
    v->addModification(this);

    // update and connect the active hyperplanes

    std::list<Hyperplane*>* hppU = u->get_activeHyperplanes(), * hppV = v->get_activeHyperplanes();
    std::list<Hyperplane*>::iterator f, g;

    for (f = hppU->begin(); f != hppU->end(); f++) {
        for (g = hppV->begin(); g != hppV->end(); g++) {
            if (*f == *g) {
                activeHyperplanes.push_back(*f);
                (*f)->addVertex(this);
                break;
            }
        }
    }
    activeHyperplanes.push_back(H);
    H->addVertex(this);
}

/*! \brief Constructor.
*/
Point::Point() : objVector(0), preImage(0), adjList(0), activeHyperplanes(0), discarded(false), onBoundingBox(false), feasible(false), degenerate(true), integer(false), integratedInUB(false), newPoint(true), checkpoint(false), isCplexChecked(false), copy(NULL), isRay(false), unboundedDimension(-1), hasRay(false), visited(false), modified(false), C(0), refVarFix(false) {}

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

/*! \brief Check whether the point is located below an hyperplane.
 *
 * CONSIDER EQUALITY AND AN EPSILON ABOVE AS BELOW.
 * \param H Hyperplane.
 * \return true, if the point is located below H.
 */
bool Point::below2(Hyperplane& H) {

    double val = 0;
    for (int k = 0; k < H.get_dim() + 1; k++) {
        val += H.get_normalVector(k) * objVector[k];
    }

    //std::cout << " dist : " << abs(val - H.get_rhs()) << "\n";

    return val <= H.get_rhs(); // epsilon value

}

/*! \brief Check whether the point is located below an hyperplane.
 *
 * \param H Hyperplane.
 * \return true, if the point is located below H.
 */
bool Point::isLocatedOn(Hyperplane& H) {

    double val = 0;
    for (int k = 0; k < H.get_dim() + 1; k++) {
        val += H.get_normalVector(k) * objVector[k];
    }

    //std::cout << " dist : " << abs(val - H.get_rhs()) << "\n";

    return abs(val - H.get_rhs()) <= EPS_PROXIMITY;
}

double Point::distToHyperplane(Hyperplane* H) {

    double val = 0;
    for (int k = 0; k < H->get_dim() + 1; k++) {
        val += H->get_normalVector(k) * objVector[k];
    }

    return val - H->get_rhs(); // negative => pts is below H
}

double Point::distToHyperplane(Hyperplane* H, std::vector<double>& yAI) {

    double val = 0;
    for (int k = 0; k < H->get_dim() + 1; k++) {
        // if it is the unbounded dimension of the ray
        if (isRay && objVector[k] > yAI[k]) {
            if (H->get_normalVector(k) != 0) // the facet is not unbounded in this direction, hence the ray necessarily above H
                val += M;
        }
        else {
            val += H->get_normalVector(k) * objVector[k];
        }
    }

    return val - H->get_rhs(); // negative => pts is below H
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

    int okpksvp = objVector.size() - 9;
    double v = pow(10, okpksvp); //objVector.size() - 9
    //return val >= H.get_rhs() && val <= H.get_rhs() + 0.00001;
    return abs(val - H.get_rhs()) <= 0.00001;//v; // the epsilon value, old is 10^-8

}

/*! \brief Compute the intersection point between an edge defined by this point, and an hyperplane.
 *
 * \param u Point. The second point that defines the edge used for the computation.
 * \param H Hyperplane. The hyperplane used for the computation.
 * \return a std::vector of double that defines the coordinates of the intersection point.
 */
std::vector<double> Point::edgeIntersection(Point& u, Hyperplane& H, int callDebug) {

    // check if edge lies on the hyperplane H
    int p = objVector.size();
    std::vector<double> edg(p);
    double orthogonalProd = 0;
    for (int k = 0; k < p; k++) {
        edg[k] = objVector[k] - u.get_objVector(k);
        orthogonalProd += edg[k] * H.get_normalVector(k);
    }
    //if (callDebug == 96)
        //std::cout << " otho prod : " << orthogonalProd << "\n";


    double lambda = 0;
    if (orthogonalProd <= 0.001) { // edge lies on H
        std::cout << " otho prod : " << orthogonalProd << "\n";
        lambda = 1;
        becomesDegenerate();
        u.becomesDegenerate();
    }
    else { // edge does not lies on H

        // compute the value of lambda, to find the point on the edge
        double denominator = 0;
        lambda = H.get_rhs();
        for (int l = 0; l < H.get_dim() + 1; l++) {
            lambda -= H.get_normalVector(l) * objVector[l];
            denominator += H.get_normalVector(l) * (u.get_objVector(l) - objVector[l]);
        }
        lambda /= denominator;

        /*if (callDebug == 164)
            std::cout << "           lambda = " << lambda << "\n";*/

        //if (lambda <= 0.001 || lambda >= 0.999) {
        //    std::cout << " lambda = " << lambda << "\n";// << " ... et pourtant, dist to H is " << u.distToHyperplane(&H) << "\n";
        //}

        // we can't divide by 0. Such case happen when both points lies on H and in particular have the same value for each
        // non-nul component of the normal vector of H. In this case, both points are degenerate, so we choose to set the new point
        // to one of the extreme points of the edges. Here lambda = 1 (the vertex considered as feasible) is choosen as it slightly
        // simplify one of the loops later in the algorithm.
        // Should be avoided because of the orthogonal product check, but that's just a security check.
        /*if (denominator == 0) {
            lambda = 1;
        }
        else {
            lambda /= denominator;
        }*/
    }

    double epsilon = 0;// 0.00001;
    if (lambda >= 1 - epsilon) {
        //std::cout << " lambda = " << lambda << "\n";
        lambda = 1;
        //u.becomesDegenerate();
    }
    else if (lambda <= 0 + epsilon) {
        //std::cout << " lambda = " << lambda << "\n";
        lambda = 0;
        //becomesDegenerate();
    }

	// compute the actual point
	std::vector<double> intersection(H.get_dim() + 1);
	for (int k = 0; k < H.get_dim() + 1; k++) {
        //if (objVector[k] == UNBOUNDED)
        //    intersection[k] = UNBOUNDED;
        //else
		    intersection[k] = lambda * u.get_objVector(k) + (1 - lambda) * objVector[k];
	}

	return intersection;

}

double Point::getEdgeCrossingValue(Point* u, Hyperplane* H) {

    double denominator = 0;
    double lambda = 0;
    lambda = H->get_rhs();
    for (int l = 0; l < H->get_dim() + 1; l++) {
        lambda -= H->get_normalVector(l) * objVector[l];
        denominator += H->get_normalVector(l) * (u->get_objVector(l) - objVector[l]);
    }
    lambda /= denominator;

    /*if (is_ray() && u->is_ray())
        std::cout << " double ray --> lambda = " << lambda << "\n";*/

    return lambda; // lambda is on u, the external vector
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

/*! \brief Change the value of the pre-image vector at a specific coordinate.
     *
     * This function set the value of variable i to value val.
     * \param i integer. The index of the variable to change.
     * \param val double. The new value of the coordinate.
     */
void Point::setPreImage(int i, double val) {
    preImage[i] = val;
}

/*! \brief Initialize the pre-image vector (allocate n cells)
 */
void Point::initPreImage(int s) {
    if (preImage.size() == 0) {
        preImage = std::vector<double>(s);
    }
}

/*! \brief Add a new point to the adjacency list.
 *
 * \param adj Point*. A pointer to the new adjacent point.
 */
void Point::addAdjacentPoint(Point* adj) {

    //std::list<Point*>::iterator vertex;
    //bool alreadyAdj = false;
    //for (vertex = adjList.begin(); vertex != adjList.end(); vertex++) {
    //    if (*vertex == adj) {
    //        //std::cout << " mer put1111111\n";
    //        alreadyAdj = true;
    //    }
    //}
    //if (adj == this) {
    //    alreadyAdj = true;
    //    //std::cout << "fdp\n";
    //}

    //if (!alreadyAdj) {
    //    adjList.push_back(adj);
    //}

    /*if (!isRay && adj->is_ray() && !adj->isNumericallyRayOf(this)) {
        std::cout << "*mumbling intensifies*\n";
    }*/

    //adjList.remove(adj);
    //if (adj != this) {
        adjList.push_back(adj);
    //}
}

/*! \brief Check whether this point is adjacent to another point y.
 *
 * \param y Point. Pointer to the point used for comparison.
 * \return true if the two points are adjacent.
 */
bool Point::isAdjacent(Point* y) {

    bool adj = false;
    std::list<Point*>::iterator vertex = adjList.begin();
    while (!adj && vertex != adjList.end()){
        if (*vertex == y) {
            adj = true;
        }
        vertex++;
    }

    return adj;
}


void Point::checkAdjacencyOld(Point* y, std::vector<Point*>& newVertices, std::vector<Point*>& degenerateVertices) {

    std::list<Hyperplane*> intersection(0);
    std::list<Hyperplane*>* hppOfY,* hppThird;
    std::list<Hyperplane*>::iterator f1, f2, f3;
    std::list<Point*>* adjY = y->get_adjList();
    std::list<Point*>::iterator v1 = adjList.begin() , v2;
    bool alreadyAdjacent = false;
    bool adjacent = true;

    // check if the two points are not already adjacent

    while (!alreadyAdjacent && v1 != adjList.end()) {
        if (*v1 == y) {
            std::cout << "we did avoid some redundant adjacency\n";
            alreadyAdjacent = true;
        }
        v1++;
    }
    if (alreadyAdjacent)
        adjacent = false;

    // compute minimal-dimensional face containing both this points and y

    if (!alreadyAdjacent) {
        hppOfY = y->get_activeHyperplanes();
        for (f1 = hppOfY->begin(); f1 != hppOfY->end(); f1++) {
            for (f2 = activeHyperplanes.begin(); f2 != activeHyperplanes.end(); f2++) {
                if (*f1 == *f2) {
                    intersection.push_back(*f1);
                }
            }
        }
    }
    

    // search for other newly generated points included in the face defined by intersection

    bool includedInIntersection;
    int nbIncludedInIntersection;
    int s = newVertices.size();
    int i = 0;

    // we start looping through each vertex
    while (adjacent && i < s) {
        if (this != newVertices[i] && y != newVertices[i]) { // we don't want to check this or y
            includedInIntersection = true;
            hppThird = newVertices[i]->get_activeHyperplanes();
            f1 = hppThird->begin();
            nbIncludedInIntersection = 0;
            // loop through the active hyperplanes of newVertices[i]
            while (adjacent && f1 != hppThird->end()) {
                f2 = intersection.begin();
                includedInIntersection = false;
                // loop through the intersection
                while (!includedInIntersection && f2 != intersection.end()) {
                    if (*f1 == *f2) { // this hyperplane of newVertices[i] is in the intersection !
                        ++nbIncludedInIntersection;
                        includedInIntersection = true;
                    }
                    f2++;
                }
                f1++;
            }
            // if all the hpp in the intersection are included in another vertex's actives hyperplanes list, this and y cannot be adjacent
            if (nbIncludedInIntersection == intersection.size()) {
                adjacent = false;
            }
        }
        i++;
    }

    // search for other degenerated points included in the face defined by intersection

    s = degenerateVertices.size();
    i = 0;

    // we start looping through each vertex
    while (adjacent && i < s) {
        if (this != degenerateVertices[i] && y != degenerateVertices[i]) { // we don't want to check this or y
            includedInIntersection = true;
            hppThird = degenerateVertices[i]->get_activeHyperplanes();
            f1 = hppThird->begin();
            nbIncludedInIntersection = 0;
            // loop through the active hyperplanes of newVertices[i]
            while (adjacent && f1 != hppThird->end()) {
                f2 = intersection.begin();
                includedInIntersection = false;
                // loop through the intersection
                while (!includedInIntersection && f2 != intersection.end()) {
                    if (*f1 == *f2) { // this hyperplane of newVertices[i] is in the intersection !
                        ++nbIncludedInIntersection;
                        includedInIntersection = true;
                    }
                    f2++;
                }
                f1++;
            }
            // if all the hpp in the intersection are included in another vertex's actives hyperplanes list, this and y cannot be adjacent
            if (nbIncludedInIntersection == intersection.size()) {
                adjacent = false;
            }
        }
        i++;
    }

    // if we conclude that this point of y are adjacent, we connect them

    if (adjacent) {
        addAdjacentPoint(y);
        y->addAdjacentPoint(this);
    }

}




void Point::checkAdjacencyOld(Point* y, std::vector<Point*>& newVertices, std::list<Point*>& all) {

    std::list<Hyperplane*> intersection(0);
    std::list<Hyperplane*>* hppOfY, * hppThird;
    std::list<Hyperplane*>::iterator f1, f2, f3;

    // compute minimal-dimensional face containing both this points and y

    hppOfY = y->get_activeHyperplanes();
    for (f1 = hppOfY->begin(); f1 != hppOfY->end(); f1++) {
        for (f2 = activeHyperplanes.begin(); f2 != activeHyperplanes.end(); f2++) {
            if (*f1 == *f2) {
                intersection.push_back(*f1);
            }
        }
    }

    // search for other newly generated points included in the face defined by intersection

    bool adjacent = true;
    bool includedInIntersection;
    int nbIncludedInIntersection;
    int s = newVertices.size();
    int i = 0;

    // we start looping through each vertex
    while (adjacent && i < s) {
        if (this != newVertices[i] && y != newVertices[i]) { // we don't want to check this or y
            includedInIntersection = true;
            hppThird = newVertices[i]->get_activeHyperplanes();
            f1 = hppThird->begin();
            nbIncludedInIntersection = 0;
            // loop through the active hyperplanes of newVertices[i]
            while (adjacent && f1 != hppThird->end()) {
                f2 = intersection.begin();
                includedInIntersection = false;
                // loop through the intersection
                while (!includedInIntersection && f2 != intersection.end()) {
                    if (*f1 == *f2) { // this hyperplane of newVertices[i] is in the intersection !
                        ++nbIncludedInIntersection;
                        includedInIntersection = true;
                    }
                    f2++;
                }
                f1++;
            }
            // if all the hpp in the intersection are included in another vertex's actives hyperplanes list, this and y cannot be adjacent
            if (nbIncludedInIntersection == intersection.size()) {
                adjacent = false;
            }
        }
        i++;
    }

    // search for other degenerated points included in the face defined by intersection

    std::list<Point*>::iterator vx = all.begin();
    // we start looping through each vertex
    while (adjacent && i < s) {
        if (this != *vx && y != *vx) { // we don't want to check this or y
            includedInIntersection = true;
            hppThird = (*vx)->get_activeHyperplanes();
            f1 = hppThird->begin();
            nbIncludedInIntersection = 0;
            // loop through the active hyperplanes of newVertices[i]
            while (adjacent && f1 != hppThird->end()) {
                f2 = intersection.begin();
                includedInIntersection = false;
                // loop through the intersection
                while (!includedInIntersection && f2 != intersection.end()) {
                    if (*f1 == *f2) { // this hyperplane of newVertices[i] is in the intersection !
                        ++nbIncludedInIntersection;
                        includedInIntersection = true;
                    }
                    f2++;
                }
                f1++;
            }
            // if all the hpp in the intersection are included in another vertex's actives hyperplanes list, this and y cannot be adjacent
            if (nbIncludedInIntersection == intersection.size()) {
                adjacent = false;
            }
        }
        vx++;
    }

    // if we conclude that this point of y are adjacent, we connect them

    if (adjacent) {
        addAdjacentPoint(y);
        y->addAdjacentPoint(this);
    }

}

/*! \brief Remove a specific point in the adjacency list.
 *
 * \param adj Point*. A pointer to the adjacent point to remove.
 */
void Point::removeAdjacentPoint(Point* adj) {
    adjList.remove(adj);
}

/*! \brief Remove a specific point in the adjacency list.
 *
 * \param adj Point*. A pointer to the adjacent point to remove.
 */
void Point::removeAdjacentPoint(std::list<Point*>::iterator adj) {
    adjList.erase(adj);
}

/*! \brief Add a new active hyperplane.
 *
 * \param id Hyperplane*. A pointer to the new active hyperplane.
 */
void Point::addActiveHyperplane(Hyperplane* H) {
    activeHyperplanes.remove(H);
	activeHyperplanes.push_back(H);
}

/*! \brief Set this point as discarded.
 */
void Point::becomesDiscarded() {
	discarded = true;
    //print();
    //std::cout << " BECOMES DISCARDED !! (unexpected condition)\n";
}

/*! \brief Set this point as discarded.
 */
void Point::becomesNonDiscarded() {
    discarded = false;
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

/*! \brief Set this point as non-newPoint, i.e. as old.
 */
void Point::becomesOld() {
    newPoint = false;
}

/*! \brief Set this point as new.
 */
void Point::becomesNew() {
    newPoint = true;
}

/*! \brief Set this point as "on the bounding box".
 */
void Point::becomesOnBoundingBox() {
	onBoundingBox = true;
}

/*! \brief Set this point as checkpoint.
 */
void Point::becomesCheckpoint() {
    checkpoint = true;
}

/*! \brief Set this point as non-checkpoint & non-cplex-checked.
 */
void Point::becomesNonCheckpoint() {
    checkpoint = false;
    isCplexChecked = false;
}

/*! \brief Note this point as checked by cplex.
 */
void Point::becomesCplexChecked() {
    isCplexChecked = true;
}

/*! \brief Check whether this point has been checked by cplex in this iteration.
 *
 * \return true if this point has been checked by cplex, false otherwise.
 */
bool Point::isCheckedByCplex() {
    return isCplexChecked;
}

/*! \brief Replace an adjacent vertex by a new one.
 *
 * \param oldVertex Point*. A pointer to the vertex to remove from the list.
 * \param newVertex Point*. A pointer to the new vertex to add to the list.
 */
void Point::replaceAdjVertex(Point* oldVertex, Point* newVertex) {

	adjList.remove(oldVertex);
    //oldVertex->removeAdjacentPoint(this);

    //adjList.remove(newVertex); // to avoid duplicates

    // this occurs when we merge a new duplicate point D with an already existing point V. At the creation,
    // D and V are said to be mutually adjacent. But when deleting D, we ask the adjacent vertex of D to consider
    // V instead of D. Because of that, we also ask V to be adjacent to itself, which creates duplicates in the
    // adjacency list of V and that we avoid by the following if statement.
    //if (newVertex != this) { 
        adjList.push_back(newVertex);
    //}
}

/*! \brief Update the set of active hyperplanes of this point.
 *
 * \param u Point&. The first point that defines the corresponding edge. Passed by reference.
 * \param v Point&. The second point that defines the corresponding edge. Passed by reference.
 * \param H Hyperplane*. A pointer to the corresponding hyperplane.
 */
bool Point::updateActiveHyperplanes(Point& u, Point& v, Hyperplane* H, Hyperplane* BB) {

	std::list<Hyperplane*>* hu = u.get_activeHyperplanes();
	std::list<Hyperplane*>::iterator itu;
	std::list<Hyperplane*>* hv = v.get_activeHyperplanes();
	std::list<Hyperplane*>::iterator itv;

    if (u.is_ray() && v.is_ray()) {
        isRay = true;
    }

    // the actual function

	for (itu = hu->begin(); itu != hu->end(); ++itu) { // for each active hyperplane of u
		for (itv = hv->begin(); itv != hv->end(); ++itv) { // for each active hyperplane of v
			if (*itu == *itv) { // if there is a common hyperplane
                /*if (isNumericallyOnHyperplane(*itu)) {
                    activeHyperplanes.push_back(*itu);
                    (*itu)->addVertex(this);
                }
                else {
                    std::cout << "bruh\n";
                }*/
                activeHyperplanes.push_back(*itu);
                (*itu)->addVertex(this);
			}
		}
	}
    //if (isNumericallyOnHyperplane(H)) {
    //    activeHyperplanes.push_back(H); // add the new hyperplane H as well
    //    H->addVertex(this);
    //}
    //else {
    //    std::cout << "big bruh\n";
    //}
    activeHyperplanes.push_back(H); // add the new hyperplane H as well
    H->addVertex(this);

    /*if (activeHyperplanes.size() <= objVector.size() - 1) {
        std::cout << "eh?\n";
    }*/

    /*if (!isRay) {
        for (itu = activeHyperplanes.begin(); itu != activeHyperplanes.end(); itu++) {
            if (*itu == BB) {
                std::cout << "mmmmmMMMMMMMM ?!?!?!\n";
            }
        }
    }*/

    /*if (u.is_ray() && v.is_ray()) {
        becomes_ray();
    }*/

    // check and correction prodecure

    //double val;
    //int nbHpp = activeHyperplanes.size();

    //// count how many of these hyperplanes actually defines this point
    //itv = activeHyperplanes.begin();
    //while (itv != activeHyperplanes.end()) {
    //    itu = itv;
    //    itv++;
    //    val = 0;
    //    for (int k = 0; k < objVector.size(); k++) {
    //        val += objVector[k] * (*itu)->get_normalVector(k);
    //    }
    //    if (abs(val - (*itu)->get_rhs()) >= 0.00001) { // the hyperplane *itu does not define this point !! We deconnect it from this point
    //        nbHpp--;
    //        (*itu)->removeVertex(this);
    //        activeHyperplanes.erase(itu);
    //    }
    //}

    //// if there are less active hyperplane than dimension, then this point is not an extreme point and lies on a face or in the interior 
    //// of a facet. This imply that u and v should not be adjacent and consequently, we disconnect them.
    //int ray = 0;
    //if (v.is_ray()) ray = 1;
    //if (nbHpp < objVector.size() - ray) {

    //    std::cout << " cheh !\n";

    //    /*u.removeAdjacentPoint(&v);
    //    v.removeAdjacentPoint(&u);

    //    for (itu = activeHyperplanes.begin(); itu != activeHyperplanes.end(); itu++) {
    //        (*itu)->print();
    //        (*itu)->removeVertex(this);
    //    }
    //    std::cout << " cheh !\n";*/
    //}

    return true; // nbHpp >= objVector.size();
}

bool Point::isNumericallyOnHyperplane(Hyperplane* H) {
    double val = 0;
    for (int k = 0; k < objVector.size(); k++) {
        val += objVector[k] * H->get_normalVector(k);
        //std::cout << " val = " << val << "\n";
    }
    //std::cout << " diff val : " << abs(val - H->get_rhs()) << "\n";
    return abs(val - H->get_rhs()) <= 0.0001;
}

bool Point::isNumericallyRayOf(Point* y) {

    int nbDivergingObjectives = objVector.size();
    for (int k = 0; k < objVector.size(); k++) {
        if (abs(objVector[k] - y->get_objVector(k)) <= 0.00000001) {
            --nbDivergingObjectives;
        }
    }

    return nbDivergingObjectives > 1;
}

/*! \brief Returns true if the point is on an already existing edge
 *
 * \param y Point*. One of the three points used for edge comparisons
 * \return true if it lies on an existing edge defined by two other points
 */
int Point::isOnEdge(Point* y) {

    int p = objVector.size();
    std::list<Point*>::iterator a = adjList.begin();
    double lambda = 0;
    double hypotheticalValue = 0;
    int k = 0;
    int l = 0;
    bool isOnEdge = false;
    bool lambdaFound = false;
    int cas = 0;

    while (!isOnEdge && a != adjList.end()) {
        if (*a != y) {
            lambdaFound = false;
            while (!lambdaFound && l < p) {
                if (abs((*a)->get_objVector(l) - y->get_objVector(l)) >= 0.00001) {
                    lambda = (objVector[l] - y->get_objVector(l)) / ((*a)->get_objVector(l) - y->get_objVector(l));
                    lambdaFound = true;
                }
                else {
                    l++;
                }
            }
            isOnEdge = true;
            k = 0;
            while (isOnEdge && k < p) {
                if (k != l) {
                    hypotheticalValue = lambda * (*a)->get_objVector(k) + (1 - lambda) * y->get_objVector(k);
                    if (abs(hypotheticalValue - objVector[k]) >= 0.00001) { // significant difference, the point is not on edge.
                        isOnEdge = false;
                    }
                }
                k++;
            }
        }
        if (!isOnEdge) {
            a++;
        }
    }

    /*a--;
    if (isOnEdge) {
        std::cout << "oof! with lambda = " << lambda << "\n";
        std::cout << "t = ";
        print();
        std::cout << "\na = ";
        (*a)->print();
        std::cout << "\ny = ";
        y->print();
        std::cout << "\n";
    }*/

    std::list<Hyperplane*>* listHpp;
    std::list<Hyperplane*>::iterator f;
    std::list<Point*>* adjacentVertex;
    std::list<Point*>::iterator vertex3;
    //a--;
    if (isOnEdge) {
        if (lambda <= -0.0000001) { // adjacent point is on the interior of the edge
            cas = 1;
            (*a)->becomesDiscarded();
        }
        else if (lambda >= 1.0000001) { // y is in the interior of the edge
            cas = 3;
            y->becomesDiscarded();
        }
        else { // this* is in the interior of the edge
            cas = 2;
            becomesDiscarded();
        }
        //std::cout << "out !\n";
    }
    else {
        cas = 0;
    }

    return cas;
}

/*! \brief Update the set of active hyperplanes of this point.
 *
 * \param u Point&. The first point that defines the corresponding edge. Passed by reference.
 * \param v Point&. The second point that defines the corresponding edge. Passed by reference.
 */
void Point::updateActiveHyperplanes(Point& u, Point& v) {

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
    //activeHyperplanes.push_back(H); // add the new hyperplane H as well
}

/*! \brief Return the degree of bounding, i.e. the number of component shared with the anti-ideal points.
 *
 * \param antiIdealPoint, vector of doubles. The coordinates of the anti-ideal point.
 * \return the degree of bounding, as an int.
 */
int Point::degreeOfBounding(std::vector<double>& antiIdealPoint) {
    
    int degree = 0;
    for (int k = 0; k < objVector.size(); k++) {
        if (objVector[k] == antiIdealPoint[k]) {
            degree++;
        }
    }

    return degree;
}

/*! \brief Return the degree of bounding, i.e. the number of component shared with the anti-ideal points.
     *
     * \param antiIdealPoint, vector of doubles. The coordinates of the anti-ideal point.
     * \return the degree of bounding, as an int.
     */
//bool Point::prelimAdjacencyCheck(std::vector<double>& antiIdealPoint, Point* y) {
//
//    int degree = 0;
//    bool canBeAdj = true;
//    for (int k = 0; k < objVector.size(); k++) {
//        if (!(abs(objVector[k] - y->get_objVector(k)) <= 0.000001)) { // if not equal
//            if (objVector[k] == antiIdealPoint[k] || y->get_objVector(k) == antiIdealPoint[k]) {
//                degree++;
//            }
//        }
//    }
//
//}

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

/*! \brief Check whether the point is on the bounding box.
 *
 * \param yI vector of double. The point defining the bounding box.
 * \return true, if the point is on the bounding box.
 */
bool Point::isOnBoundingBox(std::vector<double>& yI) {

    for (int k = 0; k < objVector.size(); k++) {
        if (yI[k] == objVector[k]) {
            onBoundingBox = true;
        }
    }

    return onBoundingBox;
}

/*! \brief Print this point.
 */
void Point::print() {
	int s = objVector.size();
    //std::cout << " is discarded: " << discarded;
	std::cout << " ( ";
	for (int k = 0; k < s - 1; k++) {
		std::cout << objVector[k] << " , ";
	}
	std::cout << objVector[s - 1] << " ) ";
}

/*! \brief Print this point.
 */
void Point::printPreImage() {
    int s = preImage.size();
    //std::cout << " is discarded: " << discarded;
    std::cout << " ( ";
    for (int i = 0; i < s - 1; i++) {
        std::cout << preImage[i] << " , ";
    }
    std::cout << preImage[s - 1] << " )\n";
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

/*! \brief Check whether the point is new, i.e. created in the current node.
 *
 * This function checks whether the point is new.
 * \return true, if the point is new.
 */
bool Point::isNew() {
    return newPoint;
}

/*! \brief Check whether the point is integer in the variable space.
 *
 * \return true if the point is integer, false otherwise.
 */
bool Point::isInteger() {

    if (!integer) {
        integer = true;
        //int i = 0;
        int n = preImage.size();
        if (n == 0) integer = false;
        //while (integer && i < n) {
        //    if (preImage[i] - floor(preImage[i] + 0.00000000001) >= 0.00000000001) { // epsilons for numerical instabilities
        //        integer = false;
        //    }
        //    ++i;
        //}
        for (int i = 0; i < n; i++) {
            if (preImage[i] - floor(preImage[i]) <= EPS_INT) { // epsilons for numerical instabilities // 0.00000001
                //integer = false;
                preImage[i] = floor(preImage[i]);
            }
            else if (preImage[i] - floor(preImage[i]) >= 1 - EPS_INT) { // 0.999999999
                //integer = false;
                preImage[i] = floor(preImage[i]) + 1;
            }
            else {
                integer = false;
            }
        }
    }

    return integer;
}

/*! \brief Check whether the point is already in the upper bound set.
 *
 * \return true, if the point is in the upper bound set.
 */
bool Point::isInUB() {
    return integratedInUB;
}

/*! \brief Check whether the point is feasible in the current node.
 *
 * \return true, if the point is feasible.
 */
bool Point::isCheckpoint() {
    return checkpoint;
}

/*! \brief Set the point as "integrated in UB"
 *
 * This function checks set integratedInUB to true, meaning that the point has been added to the upper bound set if not
 * dominated.
 */
void Point::setAsIntegratedInUB() {
    integratedInUB = true;
}

/*! \brief Register the address of the last copy of this Point.
 *
 * \param P Point*. The address of the copy.
 */
void Point::setCopy(Point* P) {
    copy = P;
}

/*! \brief Check whether this point dominates another point y
 *
 * \param y Point*. A pointer to the other point used for comparison.
 * \return true if this point dominates y.
 */
bool Point::dominates(Point* y) {

    bool domi = true;
    //bool allEqual = true;
    int k = 0;
    int p = objVector.size();

    while (domi && k < p) {
        if (objVector[k] > y->get_objVector(k) + 0.00001) {
            domi = false;
            //allEqual = false;
        }
        //else if (objVector[k] != y->get_objVector(k)) {
        //    allEqual = false;
        //}
        k++;
    }

    return domi ;
}

/*! \brief Check whether the point is feasible given the branching decisions made
 *
 * \param BD branchingDecisions. A pointer to the branching decisions used for comparison.
 * \return true if the point is feasible given the branching decisions, false otherwise.
 */
bool Point::satisfyBranchingDecisions(BranchingDecisions* BD) {
    
    bool fsb = true;

    if (preImage[BD->lastSplittedIndex] <= BD->ub[BD->lastSplittedIndex] && preImage[BD->lastSplittedIndex] >= BD->lb[BD->lastSplittedIndex]) { // decision space
        for (int k = 0; k != BD->slub.size(); k++) {
            if (!(objVector[k] <= BD->slub[k] - 1)) {
                fsb = false;
            }
        }
    }
    else {
        fsb = false;
    }

    return fsb;
}

/* Retrieve the w.s. value of the point, with weights l.
 *
 * \param vector of double, l. The weight vector.
 */
double Point::getWeightedSumValue(std::vector<double>& l) {

    double val = 0;
    for (int k = 0; k < objVector.size(); k++) {
        val += l[k] * objVector[k];
    }

    return val;
}

/*! \brief Check whether the point is significantly close to y.
 *
 * \param Point* y. The point used for comparison with this.
 * \return true if the two points are significantly close.
 */
bool Point::isVeryCloseTo(Point* y, int callDebug) {

    int dim = objVector.size();
    // 0.00000001 * dim // 0.001
    double trshld = 0.0001; // to fine tune -> 0.001 is the dist allowed for each component approx. // sqrt(pow(0.001,2) * dim)
    double dist = 0;
    bool closeOnEachDimension = true;

    for (int k = 0; k < dim; k++) {
        dist += pow(objVector[k] - y->get_objVector(k) , 2);
        if (abs(objVector[k] - y->get_objVector(k)) >= trshld) {
            closeOnEachDimension = false;
        }
    }
    dist = sqrt(dist);
    /*if (callDebug == 164)
        std::cout << "       dist = " << dist << "\n";*/

    //if (dist <= EPS_PROXIMITY)
        //std::cout << " dist : " << dist << "\n";

    return  dist <= EPS_PROXIMITY; // closeOnEachDimension;
}

/*! \brief Check whether the point is significantly close to y.
 *
 * \param Point* y. The point used for comparison with this.
 * \return true if the two points are significantly close.
 */
bool Point::isVeryCloseTo(Point* y) {

    int dim = objVector.size();
    double dist = 0;
    bool closeOnEachDimension = true;

    for (int k = 0; k < dim; k++) {
        dist += pow(objVector[k] - y->get_objVector(k), 2);
    }
    dist = sqrt(dist);

    return  dist <= EPS_PROXIMITY * 100;
}

/*! \brief Check whether the point is significantly close to y.
 *
 * \param Point* y. The point used for comparison with this.
 * \return true if the two points are significantly close.
 */
bool Point::isVeryCloseTo(std::vector<double> y) {

    int dim = objVector.size();
    double dist = 0;

    for (int k = 0; k < dim; k++) {
        dist += pow(objVector[k] - y[k], 2);
    }
    dist = sqrt(dist);

    return  dist <= EPS_PROXIMITY; // closeOnEachDimension;
}

/*! \brief Merge the point with y.
 *
 * \param Point* y. The point used for comparison with this.
 */
void Point::merge(Point* y) {

    // check if the point is close to an integer point
    std::vector<double> intCoord(objVector.size(), 0);
    int k = 0;
    bool isCloseToInt = true;
    while (isCloseToInt && k < objVector.size()) {
        if (y->get_objVector(k) - floor(y->get_objVector(k)) <= 0.001) {
            intCoord[k] = floor(y->get_objVector(k));
        }
        else if (y->get_objVector(k) - floor(y->get_objVector(k)) >= 0.999) {
            intCoord[k] = floor(y->get_objVector(k)) + 1;
        }
        else {
            isCloseToInt = false;
        }
        k++;
    }

    // if not discarded, we keep the coordinates of the old degenerate vertex. Otherwise we take the coord of the new pts.
    if (isCloseToInt) {
        /*std::cout << "it happens lol with :\n";
        print();
        std::cout << " and ";
        y->print();
        std::cout << "\n";*/
        for (int l = 0; l < objVector.size(); l++) {
            objVector[l] = intCoord[l];
        }
    }
    else if (isDiscarded()) {
        for (int k = 0; k < objVector.size(); k++) {
            objVector[k] = y->get_objVector(k);//max(y->get_objVector(k),objVector[k]);
            //max(objVector[k],y->get_objVector(k));
            //(objVector[k] + y->get_objVector(k)) / 2;
            // y->get_objVector(k);
        }
    }
    

    // merge the adjacent points
    std::list<Point*>* adjListY = y->get_adjList();
    std::list<Point*>::iterator adjY;
    std::list<Point*>::iterator adjThis;
    bool isAlreadyAdjacent;

    for (adjY = adjListY->begin(); adjY != adjListY->end(); adjY++) { // for each vertex adjY adjacent to y
        if (*adjY != this) { // that case happens when the degenerate point this is considered as feasible.
            isAlreadyAdjacent = false;
            adjThis = adjList.begin();
            while (!isAlreadyAdjacent && adjThis != adjList.end()) { // search for adjY in the adjacency list of this point
                if (*adjY == *adjThis) { // if it's already there, then we don't add the
                    isAlreadyAdjacent = true;
                }
                adjThis++;
            }
            if (!isAlreadyAdjacent) {
                adjList.push_back(*adjY);
            }
        }
    }

    // merge the active hyperplanes

    std::list<Hyperplane*>* activeHppY = y->get_activeHyperplanes();
    std::list<Hyperplane*>::iterator hppY;
    std::list<Hyperplane*>::iterator hppThis;
    bool isAlreadyActive;

    for (hppY = activeHppY->begin(); hppY != activeHppY->end(); hppY++) {
        isAlreadyActive = false;
        hppThis = activeHyperplanes.begin();
        while (!isAlreadyActive && hppThis != activeHyperplanes.end()) {
            if (*hppY == *hppThis) { // search for potential hyperplane not included in this point's list
                isAlreadyActive = true;
            }
            hppThis++;
        }
        if (!isAlreadyActive) { // if this hyperplane is not already in the list, add it
            activeHyperplanes.push_back(*hppY);
            (*hppY)->addVertex(this); // added after -> is this true ?
        }
    }
}



/*! \brief Merge the point with y.
 *
 * \param Point* y. The point used for comparison with this.
 */
void Point::merge2(Point* y) {

    // merge the adjacent points
    std::list<Point*>* adjListY = y->get_adjList();
    std::list<Point*>::iterator adjY;
    std::list<Point*>::iterator adjThis;
    bool isAlreadyAdjacent;

    /*std::cout << "\n";
    for (adjThis = adjList.begin(); adjThis != adjList.end(); adjThis++) {
        std::cout << "   -> " << *adjThis << " : ";
        (*adjThis)->print();
        std::cout << "\n";
    }*/

    for (adjY = adjListY->begin(); adjY != adjListY->end(); adjY++) { // for each vertex adjY adjacent to y
        if (*adjY != this) { // that case happens when the degenerate point this is considered as feasible.
            isAlreadyAdjacent = false;
            adjThis = adjList.begin();
            while (!isAlreadyAdjacent && adjThis != adjList.end()) { // search for adjY in the adjacency list of this point
                if (*adjY == *adjThis) { // if it's already there, then we don't add the
                    isAlreadyAdjacent = true;
                }
                adjThis++;
            }
            if (!isAlreadyAdjacent) {
                //adjList.push_back(*adjY);
                connect(*adjY);
            }
        }
    }

    // merge the active hyperplanes

    std::list<Hyperplane*>* activeHppY = y->get_activeHyperplanes();
    std::list<Hyperplane*>::iterator hppY;
    std::list<Hyperplane*>::iterator hppThis;
    bool isAlreadyActive;

    for (hppY = activeHppY->begin(); hppY != activeHppY->end(); hppY++) {
        isAlreadyActive = false;
        hppThis = activeHyperplanes.begin();
        while (!isAlreadyActive && hppThis != activeHyperplanes.end()) {
            if (*hppY == *hppThis) { // search for potential hyperplane not included in this point's list
                isAlreadyActive = true;
            }
            hppThis++;
        }
        if (!isAlreadyActive) { // if this hyperplane is not already in the list, add it
            activeHyperplanes.push_back(*hppY);
            (*hppY)->addVertex(this); // added after -> is this true ?
        }
    }
}

/*! \brief Purge a degenerate vertex by deleting itself from the list of adjacent vertex, if that happened to be the case
 *
 * \param id Hyperplane*. A pointer to the new active hyperplane.
 */
void Point::purge() {

    std::list<Point*>::iterator adjVtx;
    std::list<Hyperplane*>::iterator hpp;
    for (adjVtx = adjList.begin(); adjVtx != adjList.end(); adjVtx++) {
        if (*adjVtx == this) {
            std::cout << "We have a problem Houston.\n";
            for (hpp = activeHyperplanes.begin(); hpp != activeHyperplanes.end(); hpp++) {
                (*hpp)->removeVertex(this);
            }
        }
    }
    adjList.remove(this);
}

/*! \brief Discard the hyperplane H from the list of active hyperplanes.
 *
 * \param Hyperplane* H. The hyperplane to discard.
 */
void Point::discardHyperplane(Hyperplane* H) {

    activeHyperplanes.remove(H);
    /*std::list<Hyperplane*>::iterator hpp = activeHyperplanes.begin();
    bool isFound = false;
    while (!isFound) {
        if (*hpp == H) {
            activeHyperplanes.erase(hpp);
            isFound = true;
        }
        else {
            hpp++;
            if (hpp == activeHyperplanes.end()) {
                throw std::string("Error: Trying to discard a non-active hyperplane.");
            }
        }
    }*/
}

/*! \brief Check whether the point is an extreme ray.
 * \return true if the point is a ray, false otherwise.
 */
bool Point::is_ray() {
    return isRay;
}

/*! \brief Check whether the point is adjacent to an extreme ray.
 * \return true if the point is adjacent to an extreme ray, false otherwise.
 */
bool Point::has_ray(int k) {
    return hasRay[k];
}

/*! \brief Notify the point that it is now connected to an extreme ray.
 */
void Point::receive_ray(int k) {
    hasRay[k] = true;
}

/*! \brief Notify the point that it is no longer connected on of its extreme ray.
 */
void Point::loose_ray(int k) {
    hasRay[k] = false;
}

/*! \brief Notify the adjacent points and active hyperplane that this point will be destroyed.
 */
void Point::notifyDeletion() {

    std::list<Point*>* lalist;
    std::list<Point*>::iterator vertex, vertex2, vtx;
    std::list<Hyperplane*>::iterator f;
    bool ondebugdesrass;

    // notify adjacent points
    if (adjList.size() >= 1) {
        vertex2 = adjList.begin();
        while (vertex2 != adjList.end()) {
            vertex = vertex2;
            vertex2++;
            if (isRay)
                (*vertex)->loose_ray(unboundedDimension);
            (*vertex)->removeAdjacentPoint(this);
        }
    }
    

    // notify active hyperplanes
    for (f = activeHyperplanes.begin(); f != activeHyperplanes.end(); f++) {
        (*f)->removeVertex(this);
    }
}

void Point::connect(Point* y) {

    std::list<Point*>* adjY = y->get_adjList();;
    std::list<Point*>::iterator v;

    for (v = adjList.begin(); v != adjList.end(); v++) {
        if (*v == y)
            std::cout << " AH !\n";
    }

    for (v = adjY->begin(); v != adjY->end(); v++) {
        if (*v == this)
            std::cout << " AH !... aaaaaaaaa.. FHFSOSSHOF\n";
    }

    y->addAdjacentPoint(this);
    addAdjacentPoint(y);
}

void Point::disconnect(Point* y) {
    y->removeAdjacentPoint(this);
    removeAdjacentPoint(y);
}

/*! \brief Returns the directions of the extreme rays adjacent to this vertex
     *
     * \param yAI, vector of double. The anti-ideal point, used to detect the direction of a ray.
     * \return vector of int, containing the indices of the direction of the rays.
     */
std::vector<bool> Point::getRayDirections(std::vector<double> yAI) {

    std::vector<bool> directions(objVector.size());
    std::list<Point*>::iterator vtx;
    bool objTreated;
    
    for (int k = 0; k < objVector.size(); k++) {
        objTreated = false;
        vtx = adjList.begin();
        while (!objTreated && vtx != adjList.end()) {
            if ((*vtx)->is_ray() && (*vtx)->get_objVector(k) == yAI[k]) {
                objTreated = true;
            }
            vtx++;
        }
        if (!objTreated) directions[k] = false;
        else directions[k] = true;
    }

    return directions;
}

/*! \brief Update the set of active hyperplanes of this point.
 *
 * \param v Point&. The second point that defines the corresponding edge. Passed by reference.
 * \param H Hyperplane*. A pointer to the corresponding hyperplane.
 * \param unboundedDim int. The index of the unbounded dimension
 */
bool Point::updateActiveHyperplanes(Point& v, Hyperplane* H, int unboundedDim) {

    std::list<Hyperplane*>* hpp = v.get_activeHyperplanes();
    std::list<Hyperplane*>::iterator h;
    bool isValid = true;
    double val = 0;
    int nbDefiningHpp = hpp->size() + 1;
    std::vector<Hyperplane*> validHpp(0);

    for (h = hpp->begin(); h != hpp->end(); h++) {

        /*val = 0;
        for (int k = 0; k < objVector.size(); k++) {
            val += v.get_objVector(k) * (*h)->get_normalVector(k);
        }*/

        if (!isNumericallyOnHyperplane(*h)) { //abs(val - (*h)->get_rhs()) >= 0.0001
            --nbDefiningHpp;
        }
        else {
            validHpp.push_back(*h);
        }

        // old
        /*if ((*h)->get_normalVector(unboundedDim) == 0) {
            activeHyperplanes.push_back(*h);
            (*h)->addVertex(this);
        }*/
    }

    // eRay[k]->findMissingCoordinate(*hpp, k)
    if (false) { // less defining hyperplanes than the dimension of the space - 1, because it's an extreme ray
        // nbDefiningHpp < objVector.size() - 2
        isValid = false;
    }
    else {
        for (int i = 0; i < validHpp.size(); i++) {
            activeHyperplanes.push_back(validHpp[i]);
            validHpp[i]->addVertex(this);
        }
        //activeHyperplanes.push_back(H);
        //H->addVertex(this);
        //std::cout << "   ... and we add a new vertex to ";
        //H->print();
    }

    return isValid;
}

/*! \brief Label the point as extreme ray.
 */
void Point::becomes_ray() {
    isRay = true;
}

/*! \brief Label the point as modified.
     */
void Point::becomesModified() {
    modified = true;
}

/*! \brief Label the point as non-modified.
 */
void Point::becomesNonModified() {
    modified = false;
}

/*! \brief Label the point as visited.
 */
void Point::becomesVisited() {
    visited = true;
}

/*! \brief Label the point as non-visited.
 */
void Point::becomesNonVisited() {
    visited = false;
}

/*! \brief Add a point to the list of modifications (C).
 *
 * \param u Point*. A pointer to the point involved in the modification.
 */
void Point::addModification(Point* u) {
    C.push_back(u);
}

void Point::debug__clearAdjList() {
    adjList.clear();
}

/*! \brief Label the point as degenerate and correct possible computation made considering this point as non-degenerate.
 *
 * \param N std::list<Point*>*. The list of points computed so far at this iteration of benson's algorithm.
 */
void Point::becomesDegenerate(std::list<Point*>* N, Hyperplane* H) {
    
    std::list<Point*>::iterator s, t, a;
    std::list<Point*>* adjVtx;

    // undo previous connections related to this

    for (s = C.begin(); s != C.end(); s++) {
        adjVtx = (*s)->get_adjList();
        a = adjVtx->begin();
        while (a != adjVtx->end()) {
            t = a;
            a++;
            (*t)->disconnect(*s);
            (*t)->connect(this);
        }
        N->remove(*s);
    }

    // delete the points created with this

    for (s = C.begin(); s != C.end(); s++) {
        (*s)->notifyDeletion();
        delete* s;
    }

    // update status of this
    
    C.clear();
    discarded = false;
    degenerate = true;
    activeHyperplanes.push_back(H);
    H->addVertex(this);
}

/*! \brief Returns true if this point is included in the face defined by the hyperplanes included in F.
 *
 * \param F, list of hyperplanes. They describe the face.
 * \return true if this point is included in the face.
 */
bool Point::hasSameFacets(Point* y) {

    std::list<Hyperplane*>* hppY;
    std::list<Hyperplane*>::iterator f = activeHyperplanes.begin() , g;
    bool flag = true;

    while (flag && f != activeHyperplanes.end()) {
        flag = false;
        hppY = y->get_activeHyperplanes();
        for (g = hppY->begin(); g != hppY->end(); g++) {
            if (*f == *g) {
                flag = true;
                break;
            }
        }
        f++;
    }

    return flag;
}

/*! \brief Check whether the point is already discarded.
 *
 * This function checks whether the point is discarded.
 * \return true, if the point is discarded
 */
bool Point::isDiscarded(Hyperplane* H) {

    if (!visited) {
        if (isRay) {
            if (H->get_normalVector(unboundedDimension) == 0)
                discarded = true;
        }
        else if (below2(*H))
            discarded = true;
    }

    return discarded;
}

/*! \brief Check whether the point is already discarded.
 *
 * This function checks whether the point is discarded.
 * \return true, if the point is discarded
 */
bool Point::isDegenerate(Hyperplane* H) {

    //std::cout << " yoooodayo\n";
    if (!visited) {

        if (isRay) {
            if (H->get_normalVector(unboundedDimension) != 0)
                degenerate = false;
            else if (isLocatedOn(*H))
                degenerate = true;
        }
        else if (isLocatedOn(*H))
            degenerate = true;

        if (degenerate) {
            addActiveHyperplane(H);
            H->addVertex(this);
        }
    }

    return degenerate;
}

/*! \brief Notify the extreme rays attached to this point that they should be deleted in the concerned directions.
 */
void Point::notifyRays(Hyperplane* H, std::list<Point*>& M) {

    std::list<Point*>::iterator v;
    for (v = adjList.begin(); v != adjList.end(); v++) {
        if ((*v)->is_ray()) {
            if (H->get_normalVector((*v)->getUnboundedDim()) == 0) { // this ray is degenerated
                (*v)->becomesDegenerate();
                receive_ray((*v)->getUnboundedDim());
                (*v)->addActiveHyperplane(H);
                H->addVertex(*v);
                M.push_back(*v); // this ray is modified and should be status-reseted at the end of the iteration
                //(*v)->becomesDiscarded();
            }
            else { // this ray is not degenerated
                (*v)->becomesNonDegenerate();
                receive_ray((*v)->getUnboundedDim());
            }
        }
    }
}

/*! \brief Check adjacency of this point and y.
 *
 * This function checks whether this point and y are adjacent. The part of the polyhedron involved are described in N (new vertices)
 * and D (degenerate vertices).
 * \param y Point*. A pointer to the point used for adjacency test.
 * \param N list of Point*. The list of the polyhedron's new vertices.
 * \param D list of Point*. The list of the polyhedron's degenerate vertices.
 */
void Point::checkAdjacency(Point* y, std::list<Point*>& N, std::list<Point*>& D, int iteration) {

    std::list<Hyperplane*> F(0); // description of the face of minimal dimension that contains both this and y
    std::list<Hyperplane*>* hppY = y->get_activeHyperplanes();
    std::list<Hyperplane*>::iterator f, g;
    std::list<Point*>::iterator v = adjList.begin();
    bool flag = false;

    // check if the two points are not already adjacent

    while (flag == false && v != adjList.end()) {
        if (*v == y) {
            //std::cout << "we did avoid some redundant adjacency\n";
            flag = true;
        }
        v++;
    }

    // compute F, as the intersection of actives hyperplanes of this and y

    for (f = activeHyperplanes.begin(); f != activeHyperplanes.end(); f++) {
        for (g = hppY->begin(); g != hppY->end(); g++) {
            if (*f == *g) {
                F.push_back(*f);
                break;
            }
        }
    }

    // the actual test

    // a weaker condition is tested before doing the stronger test:
    if (F.size() >= objVector.size() - 1) { // if there is less than p - 1 hpp in the intersection, this point and y are not adjacent // objVector.size() - 1

        v = N.begin();
        while (flag == false && v != N.end()) {
            if (*v != this && *v != y && (*v)->isIncludedInFace(F)) { //  && !(*v)->isVeryCloseTo(this) && !(*v)->isVeryCloseTo(y)
                /*if (iteration == 4143 && abs(get_objVector(0) + 203.7) <= 0.1) {
                    std::cout << " adjacency cancelled by ";
                    (*v)->print();
                    std::cout << "\n";
                }*/
                flag = true;
            }
            v++;
        }

        v = D.begin();
        while (flag == false && v != D.end()) {
            if (*v != this && *v != y && (*v)->isIncludedInFace(F)) { //  && !(*v)->isVeryCloseTo(this) && !(*v)->isVeryCloseTo(y)
                /*if (iteration == 4143 && abs(get_objVector(0) + 203.7) <= 0.1) {
                    std::cout << " adjacency cancelled by ";
                    (*v)->print();
                    std::cout << "\n";
                }*/
                flag = true;
                if (!(*v)->isDegenerate())
                    std::cout << "bruh\n";
            }
            v++;
        }

        if (flag == false) {
            /*if (iteration == 4143 && abs(get_objVector(0) + 203.7) <= 0.1) {
                std::cout << " ARE ADJACENT !!\n";
            }*/
            connect(y);
        }
    }
    /*else {
        if (iteration == 4143 && abs(get_objVector(0) + 203.7) <= 0.1) {
            std::cout << " is not adj due to intersection < p\n";
        }
    }*/
}

/*! \brief Returns true if this point is included in the face defined by the hyperplanes included in F.
 *
 * \param F, list of hyperplanes. They describe the face.
 * \return true if this point is included in the face.
 */
bool Point::isIncludedInFace(std::list<Hyperplane*>& F) {

    std::list<Hyperplane*>::iterator f = F.begin(), g;
    bool flag = true;

    while (flag == true && f != F.end()) {
        flag = false;
        for (g = activeHyperplanes.begin(); g != activeHyperplanes.end(); g++) {
            if (*f == *g) {
                flag = true;
                break;
            }
        }
        f++;
    }

    return flag == true;
}

/*! \brief Clears the list of modifications (C)
 */
void Point::clearModifications() {
    C.clear();
    modified = false;
}

/* \brief Checks whether the point dominate the current
 *
 * -1 : unknown -> should be checked
 * 0 : non-dominater -> should be skipped
 * 1 : dominater -> should be considered
 * \return status of the dominance
 */
int Point::isDominater() {
    return domiSlub;
}

void Point::becomesDominater() {
    domiSlub = 1;
}
void Point::becomesNonDominater() {
    domiSlub = 0;
}

void Point::becomesUnknownDominater() {
    domiSlub = -1;
}

void Point::becomesRefVarFix() {
    refVarFix = true;
}

void Point::becomesNonRefVarFix() {
    refVarFix = false;
}

bool Point::isRefVarFix() {
    return refVarFix;
}

 /* ==========================================================
		 Getters
  ========================================================= */

  /*! \brief Return true if the point is labelled as visited.
   *
   * \return true if the point is labelled as visited.
   */
bool Point::isVisited() {
    return visited;
}


/*! \brief Return true if the point is labelled as modified.
 *
 * \return true if the point is labelled as modified.
 */
bool Point::isModified() {
    return modified;
}

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

bool Point::operator==(const Point& b) {
    bool equal = true;
    int p = objVector.size();
    int k = 0;
    while (equal && k < p) {
        if (objVector[k] != b.objVector[k]) {
            equal = false;
        }
        k++;
    }
    return equal;
}

/*! \brief Returns the number of variables.
 *
 * \return the number of variables, as an int.
 */
int Point::get_nbVar() {
    return preImage.size();
}

/*! \brief Returns the number of objectives.
 *
 * \return the number of objectives, as an int.
 */
int Point::get_nbObj() {
    return objVector.size();
}

/* \brief Returns the unbounded dimension of the ray.
 *
 * \return the dimension, as an int.
 */
int Point::getUnboundedDim() {
    return unboundedDimension;
}
