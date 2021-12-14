#include "Lub.h"

// ===============================================================================================================================
//							LocalUpperBound
// ===============================================================================================================================

/* ==========================================================
        Constructors
 ========================================================= */

 /*! \brief Default constructor of a local upper bound.S
  */
LocalUpperBound::LocalUpperBound() : coordinates(0), definingPoints(0), id(-1), discarded(false) {}

 /*! \brief Default constructor of a local upper bound.
  *
  * \param lp MathematicalModel. Used to get the number of objectives.
  */
LocalUpperBound::LocalUpperBound(MathematicalModel& lp) : coordinates(lp.get_p(), INT_MAX), id(0), discarded(false) {
    definingPoints = std::vector<std::list<Solution*>>(lp.get_p());
    for (int k = 0; k < lp.get_p(); k++) {
        definingPoints[k] = std::list<Solution*>(0);
    }
}

/*! \brief Constructor of a local upper bound defined by an old one u and a new point y on objective j.
     *
     * This function creates a local upper bounds with the same coordinates and defining points as u, except on objective
     * j, where z is used as defining point and defines the jth coordinates too.
     * \param y Solution. The new solution that defines the local upper bound.
     * \param u LocalUpperBound. The old local upper boud used to defines the other p - 1 coordinates.
     * \param j int. The objective defined by the new solution y.
     */
LocalUpperBound::LocalUpperBound(Solution* y, LocalUpperBound& u, int j, int id) : id(id), discarded(false) {

    //std::cout << "building new lub...\n";
    int p = u.get_dim();
    std::list<Solution*>* defPts;
    std::list<Solution*>::iterator yDef;
    
    coordinates = std::vector<int>(p);
    definingPoints = std::vector<std::list<Solution*>>(p);

    for (int k = 0; k < p; k++) {
        definingPoints[k] = std::list<Solution*>(0);
        if (k == j) {
            coordinates[k] = y->get_objVector(k);
            definingPoints[k].push_back(y);
            //std::cout << "init def pts!\n";
            //std::cout << " new def coord: from address " << y;
            //y->print();
        }
        else {
            coordinates[k] = u.get_coordinate(k);
            defPts = u.get_definingPoints(k);
            //std::cout << "copy def pts!\n";
            for (yDef = defPts->begin(); yDef != defPts->end(); yDef++) {
                //std::cout << " new def pts: ";
                //(*yDef)->print();
                if ( !(*yDef)->isDiscarded() && (*yDef)->get_objVector(j) < y->get_objVector(j)) {
                    definingPoints[k].push_back(*yDef);
                }
            }
        }
    //std::cout << "defining pts on component " << k << ": " << definingPoints[k].size() << std::endl;
    }
}

/* ==========================================================
        Regular Methods
 ========================================================= */

 /*! \brief Prints the local upper bound.
  */
void LocalUpperBound::print() {
    std::cout << "( ";
    for (int i = 0; i < coordinates.size() - 1; i++) {
        std::cout << coordinates[i] << " , ";
    }
    std::cout << coordinates[coordinates.size() - 1] << " )\n";
}

/*! \brief Check whether this Local Upper Bound is dominated by a Solution y.
 *
 * \param y Solution. The Solution used for comparison.
 * \return true if this local upper bound is strictly dominated by y
 */
bool LocalUpperBound::isStrictlyDominated(Solution* y) {
    bool dominated = true;
    int k = 0;
    int p = coordinates.size();
    while (dominated && k < p) {
        if (coordinates[k] <= y->get_objVector(k)) {
            dominated = false;
        }
        ++k;
    }

    return dominated;
}

/*! \brief Check whether this Local Upper Bound is defined by a Solution y on a given component k.
 *
 * \param y Solution. The Solution used for comparison.
 * \param k int. The index of the coordinate.
 * \return true if this local upper bound is defined by y
 */
bool LocalUpperBound::isDefinedBy(Solution* y, int k) {

    //if (k == 1 && y->get_objVector(0) == -29 && y->get_objVector(1) == -21 && y->get_objVector(2) == -30) {
    //    std::cout << "on passe ici\n";
    //    y->print();
    //    std::cout << " vs ";
    //    print();
    //}

    bool definedBy = false;
    //std::cout << " COMP " << k << " : ";
    //print();
    //std::cout << " vs ";
    //y->print();

    if (coordinates[k] == y->get_objVector(k)) {
        //std::cout << "true, equal... ";
        definedBy = true;
        //if (k == 1 && y->get_objVector(0) == -29 && y->get_objVector(1) == -21 && y->get_objVector(2) == -30) {
        //    std::cout << "check -> ";
        //    y->print();
        //    std::cout << " vs ";
        //    print();
        //}
        for (int i = 0; i < coordinates.size(); i++) {
            //if (k == 1 && y->get_objVector(0) == -29 && y->get_objVector(1) == -21 && y->get_objVector(2) == -30) {
            //    std::cout << "\n\n\n\n\n objective " << i << " -> coord = " << coordinates[i] << " , solution = " << y->get_objVector(i);
            //}
            if (i != k && coordinates[i] <= y->get_objVector(i)) { // strict ?
                definedBy = false;
            }
            //else {
            //    std::cout << i << " true, defined ... ";
            //}
        }
    }
    //std::cout << "\n";

    return definedBy;
}

/*! \brief Add a new defining component of this local upper bound on coordinate k
 *
 * This function adds the solution y as new defining component of this local upper bound on objective k.
 * \param y Solution. The Solution used for comparison.
 * \param k int. The index of the objective.
 */
void LocalUpperBound::addDefiningComponent(Solution* y, int k) {
    definingPoints[k].push_back(y);
}

/*! \brief Filter the discarded defining solutions on coordinate k
 *
 * This function remove from the list of defining solutions of coordinate k the solution that are discarded, i.e. dominated.
 * \param k int. The index of the objective.
 */
void LocalUpperBound::filterDefiningComponent(int k) {

    std::list<Solution*>::iterator y;
    std::list<Solution*>::iterator yNext = definingPoints[k].begin();

    while (yNext != definingPoints[k].end()) {
        y = yNext;
        yNext++;
        if ((*y)->isDiscarded()) {
            definingPoints[k].erase(y);
        }
    }
}

/*! \brief Compute the critical value z^{max}_j in the computation of the local upper bounds
 * 
 * \param j int. The index of the coordinate we are looking at.
 * \return the value of the critical value.
 */
int LocalUpperBound::computeCriticalValue(int j) {

    int p = coordinates.size();
    int max;
    int min;
    std::list<Solution*>::iterator y;
    std::list<Solution*>::iterator nextY;
    std::list<Solution*>::iterator zbleh;

    max = INT_MIN;
    for (int k = 0; k < p; k++) {
        //std::cout << " coord " << k << ", nb of defining points: " << definingPoints[k].size() << std::endl;
        if (k != j && definingPoints[k].size() >= 1) {
            min = INT_MAX;
            nextY = definingPoints[k].begin();
            //for (y = definingPoints[k].begin(); y != definingPoints[k].end(); y++) {
            while (nextY != definingPoints[k].end()) {
                y = nextY;
                nextY++;
                //else {
                    //std::cout << " -> obj " << k << ", pt is: ";
                    //(*y)->print();
                if ((*y)->get_objVector(j) < min) {
                    //(*y)->print();
                    min = (*y)->get_objVector(j);
                }
                //}
                //if ((*y)->isDiscarded()) {
                //    definingPoints[k].erase(y);
                    //if (definingPoints[k].empty()) {
                    //    min = INT_MIN; // case where in fact, the list is empty, i.e. we should not have entered this loop and thus set min to INT_MAX
                    //}
                //}
            }
            if (min > max) {
                max = min;
            }
        }
    }

    return max;
}

/*! \brief Set this local upper bound as discarded.
 */
void LocalUpperBound::becomesDiscarded() {
    discarded = true;
}

/*! \brief Check whether the local upper bound is already discarded.
 *
 * \return true if the local upper bound is already discarded.
 */
bool LocalUpperBound::isDiscarded() {
    return discarded;
}

/*! \brief Check whether a local upper bound is located above the hyperplane.
 *
 * \param H Hyperplane. A pointer to the hyperplane used for comparison.
 * \param P Parameters. A pointer to the parameters of the algorithm, to access the GCD of each objective.
 * \return a pointer to the list of defining Solutions.
 */
bool LocalUpperBound::above(Hyperplane* H, Parameters* P) {

    bool above = false;
    double lhs = 0;
    double epsilon = 0.001; // 0.001;

    for (int k = 0; k <= H->get_dim(); k++) {
        lhs += H->get_normalVector(k) * (coordinates[k] - 1); //P->GCD[k]
    }

    if (lhs + epsilon >= H->get_rhs()) {
        above = true;
    }

    return above;
}

bool LocalUpperBound::above(Hyperplane* H) {

    bool above = false;
    double lhs = 0;
    double epsilon = 0.001; // 0.001;

    for (int k = 0; k <= H->get_dim(); k++) {
        lhs += H->get_normalVector(k) * (coordinates[k] - 1); //P->GCD[k]
    }

    if (lhs + epsilon >= H->get_rhs()) {
        above = true;
    }

    return above;
}

/*! \brief Compute the weigthed sum of the components in the objective space.
 *
 * \param vector of double l. The weight vector.
 */
double LocalUpperBound::getWeightedSum(std::vector<double> l) {

    double ws = 0;
    for (int k = 0; k < coordinates.size(); k++) {
        ws += l[k] * coordinates[k];
    }

    return ws;
}

bool LocalUpperBound::isRedundant(std::list<LocalUpperBound>& NU) {

    std::list<LocalUpperBound>::iterator u;
    bool isDomi = true;

    for (u = NU.begin(); u != NU.end(); u++) {
        if (&(*u) != this) {
            isDomi = true;
            for (int k = 0; k != coordinates.size(); k++) {
                if (coordinates[k] < u->get_coordinate(k)) {
                    isDomi = false;
                    break;
                }
            }
            if (isDomi) {
                std::cout << " AH ! ";
                print();
                std::cout << " is domi by ";
                u->print();
                std::cout << "\n";
            }
        }   
    }

    return isDomi;
}

/* ==========================================================
        Getters
 ========================================================= */

/*! \brief Return the number of objectives.
 */
int LocalUpperBound::get_dim() {
    return coordinates.size();
}

/*! \brief Return the objective vector as a pointer.
 *
 * This function returns a pointer to the coordinate vector.
 * \return a pointer to the vector of coordinates
 */
int LocalUpperBound::get_coordinate(int k) {
    return coordinates[k];
}

/*! \brief Return the defining-points vector as a pointer.
*
* This function returns a pointer to the defining-points vector.
* \return a pointer to the vector of defining-points
*/
std::list<Solution*>* LocalUpperBound::get_definingPoints(int k) {
    return &definingPoints[k];
}

/*! \brief Return the id of the local upper bound.
 *
 * This function returns the id of this local upper bounx
 * \return the value of the id, as an int
 */
int LocalUpperBound::get_id() {
    return id;
}