/**
* \author Nicolas Forget
*
* This class describe an hyperplane. Used in the linear-relaxation-based lower bound sets.
*/

#pragma once
#include <iostream>
#include <vector>
#include<list>

class Point;

class Hyperplane
{
private:
	std::vector<double> normalVector; //!< normal vector of the hyperplane
	int nbDefiningPts; //!< number of points known on the hyperplane - used if the hyperplane defines a facet of a polyhedron
	std::list<Point*> defPts; //!< defining points of the facet
	int dim; //!< dimension of the hyperplane (p - 1)
	double rhs; //!< a constant defining the right-hand side of the hyperplane's equation
	Hyperplane* copy; //!< used for copy purpose only. See copy constructor of LinearRelaxation for its purpose.
	bool redundant; //!< true if the hyperplane is redundant, i.e. if it is a face but not a facet of the lower bound set.
	bool isNew; //!< true if computed at the current iteration of the BB

public:

	/*! \brief Constructor with normal vector and point
	 *
	 * This function builds the an hyperplane when its normal vector and a point of this hyperplane are given.
	 * \param nv std::vector of double. It defines the normal vector of the hyperplane.
	 * \param pt std::vector of double. It defines the coordinates of a point of the hyperplane.
	 */
	Hyperplane(std::vector<double> const& nv, std::vector<double> const& pt);

	/*! \brief Constructor with normal vector and right-hand side
	 *
	 * This function builds the an hyperplane when its normal vector and its right-hand side are given.
	 * \param nv std::vector of double. It defines the normal vector of the hyperplane.
	 * \param rhs double. It defines the right-hand side of the hyperplane's equation.
	 */
	Hyperplane(std::vector<double> const& nv, double const& rhs);

	/*! \brief Check whether the hyperplane is redundant
	 *
	 * This function checks whether the hyperplane is redundant in the current description of the lower bound set.
	 * This is done by checking whether the number of extreme points of the lower bound set located on this hyperplane
	 * is lower than its dimension (dim); or if its defining points belongs to a single face of the lower bound set.
	 * If this is true, the hyperplane defines a face (not a facet) of the lower
	 * bound set and is thus redundant. In this case, redundant is set to true.
	 * \param int p. Dimension of the objective space.
	 */
	void checkRedundancy();

	/*! \brief Check whether the hyperplane is redundant using a quick but approximate method.
	 *
	 * This function quickly checks whether the hyperplane is redundant in the current description of the lower bound set.
	 * This is done by checking whether the number of extreme points of the lower bound set located on this hyperplane
	 * is lower than its dimension (dim). 
	 * \return true if it is redundant.
	 */
	bool quickCheckRedundancy();

	/*! \brief Check whether the hyperplane is redundant
	 *
	 * This functions states whether the hyperplane is redundant, by returning the value of redundant.
	 * \return true if the hyperplane is redundant.
	 */
	bool isRedundant();

	/*! \brief Prints the hyperplane
	 *
	 * This function prints the hyperplane in the command line by giving its normal vector and its right-hand side in
	 * an easily-readible way for the user.
	 */
	void print();

	/*! \brief Check whether the point y is dominated by the hyperplane.
	 *
	 * This function checks if the point y is located above the hyperplane, i.e. if there is at least one point of the
	 * hyperplane that dominates it.
	 * \param y vector of int. Represents the point used for the comparison.
	 * \return true if the point is dominated by the hyperplane.
	 */
	bool dominates(std::vector<int>& y);

	/*! \brief Returns the right-hand side of the hyperplane
	 *
	 * This function returns the right-hand side of the hyperplane.
	 * \return the right-hand side of the hyperplane.
	 */
	double get_rhs();

	/*! \brief Returns the dimension of the hyperplane
	 *
	 * This function returns the dimension of the hyperplane.
	 * \return the dimension of the hyperplane.
	 */
	int get_dim();

	/*! \brief Returns a specific coordinate of the normal vector of the hyperplane.
	 *
	 * This function returns the coordinate of the normal vector at index coord.
	 * \param coord integer. Defines the index of the coordinate we want to look at.
	 * \return the coordinate of the normal vector at index coord.
	 */
	double get_normalVector(int coord);

	/*! \brief Returns the number of points defining the hyperplane.
	 *
	 * This function returns the number of extreme points of the current lower bound set that are located on this
	 * hyperplane.
	 * \return the number of points defining the hyperplane.
	 */
	int get_nbDefiningPts();

	/*! \brief Counts that a new vertex defines this hyperplane
	 *
	 * This function increases by one the number of extreme points of the lower bound set that are located on this
	 * hyperplane.
	 */
	void addVertex();

	/*! \brief Counts that a new vertex defines this hyperplane
	 *
	 * This function increases by one the number of extreme points of the lower bound set that are located on this
	 * hyperplane and add the vertex to the list of defining vertex.
	 * \param vertex Point*. Pointer to the vertex added to the list of defining points.
	 */
	void addVertex(Point* vertex);

	/*! \brief Counts that a vertex defining this hyperplane doesn't exists anymore.
	 *
	 * This function decreases by one the number of extreme points of the lower bound set that are located on this
	 * hyperplane.
	 */
	void removeVertex();

	/*! \brief Counts that a vertex defining this hyperplane doesn't exists anymore.
	 *
	 * This function decreases by one the number of extreme points of the lower bound set that are located on this
	 * hyperplane and remove the vertex to the list of defining vertex.
	 * \param vertex Point*. Pointer to the vertex removed from the list of defining points.
	 */
	void removeVertex(Point* vertex);

	/*! \brief Register the address of the last copy of this hyperplane.
	 *
	 * This function register the address of the last copy of thus hyperplane. This is done for copying LinearRelaxations
	 * and link the pointers properly.
	 * \param H Hyperplane*. The address of the copy.
	 */
	void setCopy(Hyperplane* H);

	/*! \brief Returns the address of the last copy of this hyperplane.
	 *
	 * This function returns a pointer to the last copy of this hyperplane.
	 * \return a pointer to the copy.
	 */
	Hyperplane* get_copy();

	/*! \brief Returns the address of the list of defining points
	 *
	 * This function returns a pointer to the list of defining points
	 * \return a pointer to the list
	 */
	std::list<Point*>* get_defPts();

	/*! \brief Notify the defining points that this hyperplane will be deleted.
	 */
	void notifyDeletion();

	/*! \brief Compute the intersection point between an edge defined by this point, and an hyperplane.
	*
	* \param u Point*. A pointer to the first point that defines the edge used for the computation.
	* \param v Point*. A pointer to the second point that defines the edge used for the computation.
	* \return lambda, int, s.t. the new point y = lambda * u + (1 - lambda) * v
	*/
	double edgeIntersection(Point* u, Point* v);

	/*! \brief Compute the intersection point between an edge defined by this point, and an hyperplane.
	*
	* \param u Point*. A pointer to the first point that defines the edge used for the computation.
	* \param v Point*. A pointer to the second point that defines the edge used for the computation.
	* \return lambda, int, s.t. the new point y = lambda * u + (1 - lambda) * v
	*/
	std::vector<double> edgeIntersection2(Point* u, Point* v);

	bool is_new();
	void becomes_new();
	void becomes_old();
	std::vector<double>* getNormalVector();
};

