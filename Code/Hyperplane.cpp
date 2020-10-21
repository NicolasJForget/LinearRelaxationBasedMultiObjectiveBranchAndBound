#include "Hyperplane.h"


/* ==========================================================
		Constructors
 ========================================================= */

/*! \brief Constructor with normal vector and point
 *
 * \param nv std::vector of double. It defines the normal vector of the hyperplane.
 * \param pt std::vector of double. It defines the coordinates of a point of the hyperplane.
 */
Hyperplane::Hyperplane(std::vector<double> const& nv, std::vector<double> const& pt) : normalVector(nv), nbDefiningPts(0), dim(nv.size() - 1), copy(NULL) {

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
Hyperplane::Hyperplane(std::vector<double> const& nv, double const& rhs) : normalVector(nv), nbDefiningPts(0), dim(nv.size() - 1), rhs(rhs), copy(NULL) {

}


/* ==========================================================
		Regular Methods
 ========================================================= */

 /*! \brief Check whether the hyperplane is redundant
  *
  * \return true if the hyperplane is redundant in the representation of the lower bound set , false otherwise.
  */
bool Hyperplane::isRedundant() {
	return nbDefiningPts <= dim;
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

/*! \brief Counts that a vertex defining this hyperplane doesn't exists anymore.
 */
void Hyperplane::removeVertex() {
	nbDefiningPts--;
}

/*! \brief Register the address of the last copy of this hyperplane.
 *
 * \param H Hyperplane*. The address of the copy.
 */
void Hyperplane::setCopy(Hyperplane* H) {
    copy = H;
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