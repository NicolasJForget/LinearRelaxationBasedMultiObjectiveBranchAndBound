#include "solution.h"

/* ==========================================================
        Constructors
 ========================================================= */

/*! \brief Default constructor of a solution.
 */
Solution::Solution() : variableVector(0), objectiveVector(0), dominated(false) {}

/*! \brief Construct a solution from an extreme point of the LinearRelaxation
 *
 * This function creates a solution from an extreme point of the LinearRelaxation, by extracting its pre-image and
 * objective vector.
 * \param pts Point. The point the solution is created from.
 */
Solution::Solution(Point& pts) : dominated(false) {

    // fill solution vector
    int n = pts.get_nbVar();
    variableVector = std::vector<int>(n);
    for (int i = 0; i < n; i++) {
        variableVector[i] = (int)round(pts.get_preImage(i));
    }
    
    // fill objective vector
    int p = pts.get_nbObj();
    objectiveVector = std::vector<int>(p);
    for (int k = 0; k < p; k++) {
        objectiveVector[k] = (int)round(pts.get_objVector(k));
    }
}

/* ==========================================================
        Regular Methods
 ========================================================= */

 /*! \brief Check whether this Solution dominates another Solution y.
  *
  * \param y Solution. The Solution used for comparison.
  */
bool Solution::dominate(Solution& y) {
    
    bool domi = true;
    int p = objectiveVector.size();
    int i = 0;
    
    do {
        if (y.get_objVector(i) < objectiveVector[i]) {
            domi = false;
        }
        i++;
    } while (domi && i < p);
        
    return domi;
}

/*! \brief Prints the solution.
 */
void Solution::print() {

    int s = objectiveVector.size();
    int s2 = objectiveVector.size() - 1;

    if (s2 >= 0) { // objectiveVector.size() - 1
        std::cout << " ( ";
        for (int i = 0; i < objectiveVector.size() - 1; i++) {
            std::cout << objectiveVector[i] << " , ";
        }
        std::cout << objectiveVector[objectiveVector.size() - 1] << " )\n";
    }
    else {
        std::cout << " xxx \n";
    }
}

/*! \brief Discard the solution.
 *
 * Discard this solution by setting dominated to true.
 */
void Solution::discard() {
    dominated = true;
}

/* ==========================================================
        Getters
 ========================================================= */

/*! \brief Returns the value of a specific objective function.
     *
     * This function returns the value of the objective vector at coordinate obj. This correspond to the value of
     * this point in objective obj.
     * \param obj integer. The index of the objective to look at.
     * \return the value of this objective, as a double.
     */
int Solution::get_objVector(int obj) {
    return objectiveVector[obj];
}

/*! \brief Checks whether the solution is discarded.
 *
 * \return true if the solution is discarded, false otherwise.
 */
bool Solution::isDiscarded() {
    return dominated;
}