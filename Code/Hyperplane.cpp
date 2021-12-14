#include "Hyperplane.h"
#include "Point.h"


/* ==========================================================
		Constructors
 ========================================================= */

/*! \brief Constructor with normal vector and point
 *
 * \param nv std::vector of double. It defines the normal vector of the hyperplane.
 * \param pt std::vector of double. It defines the coordinates of a point of the hyperplane.
 */
Hyperplane::Hyperplane(std::vector<double> const& nv, std::vector<double> const& pt) : normalVector(nv), nbDefiningPts(0), defPts(0), dim(nv.size() - 1), copy(NULL), redundant(false), isNew(true) {

	// compute the right-hand side given the normal vector and a point
	rhs = 0;
	for (int i = 0; i < dim + 1; i++) {
		rhs += nv[i] * pt[i];
	}
}

/*! \brief Constructor with normal vector and right-hand side
 *
 * \param nv std::vector of double. It defines the normal vector of the hyperplane.
 * \param rhs double. It defines the right-hand side of the hyperplane's equation.
 */
Hyperplane::Hyperplane(std::vector<double> const& nv, double const& rhs) : normalVector(nv), nbDefiningPts(0), defPts(0), dim(nv.size() - 1), rhs(rhs), copy(NULL), redundant(false), isNew(true) {

}


/* ==========================================================
		Regular Methods
 ========================================================= */

 /*! \brief Check whether the hyperplane is redundant using a quick but approximate method.
  *
  * \return true if it is redundant.
  */
bool Hyperplane::quickCheckRedundancy() {
    return defPts.size() <= dim;
}

 /*! \brief Check whether the hyperplane is redundant
  *
  * \return true if the hyperplane is redundant in the representation of the lower bound set , false otherwise.
  */
void Hyperplane::checkRedundancy() {

    // if less points than dimension, we are left with only a face (dim between 0 and p - 2) and so it is redundant
    if (defPts.size() <= dim) {
        if (normalVector[1] == 1) {
            //std::cout << "Deleted by defPts condition : ";
            //print();
        }
        redundant = true;
    }
    else { // if (dim >= 3) {

        std::list<Point*>::iterator vertex = defPts.begin();
        redundant = true;

        while (redundant && vertex != defPts.end()) {
            if (!(*vertex)->isDegenerate()) {
                redundant = false;
            }
            vertex++;
        }
        if (redundant && normalVector[1] == 1) {
            //std::cout << "Deleted by total degenerancy condition : ";
            //print();
        }
    }
}

/*! \brief Compute the intersection point between an edge defined by this point, and an hyperplane.
 *
 * \param u Point*. A pointer to the first point that defines the edge used for the computation.
 * \param v Point*. A pointer to the second point that defines the edge used for the computation.
 * \return lambda, int, s.t. the new point y = lambda * u + (1 - lambda) * v
 */
double Hyperplane::edgeIntersection(Point* u, Point* v) {

    double denominator = 0;
    double lambda = rhs;

    for (int l = 0; l < dim + 1; l++) {
        lambda -= normalVector[l] * v->get_objVector(l);
        denominator += normalVector[l] * (u->get_objVector(l) - v->get_objVector(l));
    }
    lambda /= denominator;
    
    return lambda; // y = lamda * u + (1 - lambda) * v
}

// formula :
// lamda * u + (1 - lambda) * v = y AND n . y = d, where y is the new point
// n . ( l * u + (1 - l) * v ) = d
// n . ( l u + v - l v ) = d
// n . ( l (u - v) + v ) = d
// n . ( l (u - v) ) + n . v = d
// l * n . (u - v) = d - n . v
// l = (d - n . v) / [n . (u - v)]

/*! \brief Compute the intersection point between an edge defined by this point, and an hyperplane.
 *
 * \param u Point*. A pointer to the first point that defines the edge used for the computation.
 * \param v Point*. A pointer to the second point that defines the edge used for the computation.
 * \return lambda, int, s.t. the new point y = lambda * u + (1 - lambda) * v
 */
std::vector<double> Hyperplane::edgeIntersection2(Point* u, Point* v) {

    std::vector<double> res(u->get_nbObj());
    double denominator = 0;
    double lambda = rhs;

    for (int l = 0; l < dim + 1; l++) {
        lambda -= normalVector[l] * v->get_objVector(l);
        denominator += normalVector[l] * (u->get_objVector(l) - v->get_objVector(l));
    }
    lambda /= denominator;

    for (int k = 0; k < u->get_nbObj(); k++)
        res[k] = lambda * u->get_objVector(k) + (1 - lambda) * v->get_objVector(k);

    return res; // y = lamda * u + (1 - lambda) * v
}

/*! \brief Check whether the hyperplane is redundant
 *
 * This functions states whether the hyperplane is redundant, by returning the value of redundant.
 * \return true if the hyperplane is redundant.
 */
bool Hyperplane::isRedundant() {
    return redundant;
}

/*! \brief Prints the hyperplane
 */
void Hyperplane::print() {

	std::cout << "( ";
	for (int k = 0; k < dim; k++) {
		std::cout << normalVector[k] << " , ";
	}
	std::cout << normalVector[dim] << " ) " << " = " << rhs << std::endl;

}

/*! \brief Counts that a new vertex defines this hyperplane
 */
void Hyperplane::addVertex() {
	nbDefiningPts++;
}

/*! \brief Add a new defining vertex to this hyperplane
 * \param vertex Point*. Pointer to the vertex added to the list of defining points.
 */
void Hyperplane::addVertex(Point* vertex) {
    //defPts.remove(vertex);
    defPts.push_back(vertex);
    nbDefiningPts++;
}

/*! \brief Counts that a vertex defining this hyperplane doesn't exists anymore.
 */
void Hyperplane::removeVertex() {
	nbDefiningPts--;
}

/*! \brief Remove the vertex from the list of defining points of this hyperplane
 * \param vertex Point*. Pointer to the vertex removed from the list of defining points.
 */
void Hyperplane::removeVertex(Point* vertex) {

    nbDefiningPts--;
    defPts.remove(vertex);

    /*std::list<Point*>::iterator pts = defPts.begin();
    bool erased = false;

    while (!erased && pts != defPts.end()) {
        if (*pts == vertex) {B
            defPts.erase(pts);
            erased = true;
            nbDefiningPts--;
        }
        else {
            pts++;
        }
    }*/
    //do {
    //    if (*pts == vertex) {
    //        defPts.erase(pts);
    //        erased = true;
    //    }
    //    else {
    //        pts++;
    //        if (pts == defPts.end()) {
    //            erased = true;
    //            //throw std::string("delete a non-defining point");
    //        }
    //    }
    //} while (!erased); //pts != defPts.end() || 
    //nbDefiningPts--;
}

/*! \brief Register the address of the last copy of this hyperplane.
 *
 * \param H Hyperplane*. The address of the copy.
 */
void Hyperplane::setCopy(Hyperplane* H) {
    copy = H;
}

/*! \brief Check whether the point y is dominated by the hyperplane.
 *
 * \param y vector of int. Represents the point used for the comparison.
 * \return true if the point is dominated by the hyperplane.
 */
bool Hyperplane::dominates(std::vector<int>& y) {

    double val = 0;
    for (int k = 0; k <= dim; k++) {
        val += y[k] * normalVector[k];
    }

    return val >= rhs;
}

/*! \brief Notify the defining points that this hyperplane will be deleted.
 */
void Hyperplane::notifyDeletion() {

    std::list<Point*>::iterator vertex;
    for (vertex = defPts.begin(); vertex != defPts.end(); vertex++) {
        (*vertex)->discardHyperplane(this);
    }
}

bool Hyperplane::is_new() {
    return isNew;
}

void Hyperplane::becomes_new() {
    isNew = true;
}
void Hyperplane::becomes_old() {
    isNew = false;
}

/* ==========================================================
		Getters
 ========================================================= */

 /*! \brief Returns the right-hand side of the hyperplane
  *
  * \return the right-hand side of the hyperplane.
  */
double Hyperplane::get_rhs() {
	return rhs;
}

/*! \brief Returns the dimension of the hyperplane
 *
 * \return the dimension of the hyperplane.
 */
int Hyperplane::get_dim() {
	return dim;
}

/*! \brief Returns a specific coordinate of the normal vector of the hyperplane.
 *
 * \param coord integer. Defines the index of the coordinate we want to look at.
 * \return the coordinate of the normal vector at index coord.
 */
double Hyperplane::get_normalVector(int coord) {
	return normalVector[coord];
}

/*! \brief Returns the number of points defining the hyperplane.
 *
 * \return the number of points defining the hyperplane.
 */
int Hyperplane::get_nbDefiningPts() {
	return nbDefiningPts;
}

/*! \brief Returns the address of the last copy of this hyperplane.
 *
 * \return a pointer to the copy.
 */
Hyperplane* Hyperplane::get_copy() {
    return copy;
}

/*! \brief Returns the address of the list of defining points
 *
 * \return a pointer to the list
 */
std::list<Point*>* Hyperplane::get_defPts() {
    return &defPts;
}

std::vector<double>* Hyperplane::getNormalVector() {
    return &normalVector;
}