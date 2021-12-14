#include "solution.h"

/* ==========================================================
        Constructors
 ========================================================= */

/*! \brief Default constructor of a solution.
 */
Solution::Solution() : variableVector(0), objectiveVector(0), dominated(false), cpuCreation(0) {}

Solution::Solution(std::vector<int>& y) : variableVector(0), objectiveVector(y), dominated(false), cpuCreation(0) {}

/*! \brief Construct a solution from an extreme point of the LinearRelaxation
 *
 * This function creates a solution from an extreme point of the LinearRelaxation, by extracting its pre-image and
 * objective vector.
 * \param pts Point. The point the solution is created from.
 */
Solution::Solution(Point& pts, MathematicalModel* lp) : dominated(false), cpuCreation(0) {

    // fill solution vector
    //int n = pts.get_nbVar();
    int n = lp->get_n();
    variableVector = std::vector<int>(n);
    for (int i = 0; i < n; i++) {
        variableVector[i] = (int)round(pts.get_preImage(i));
    }
    
    // fill objective vector
    int p = pts.get_nbObj();
    objectiveVector = std::vector<int>(p);
    for (int k = 0; k < p; k++) {
        for (int i = 0; i < n; i++) {
            objectiveVector[k] += (int)round(pts.get_preImage(i) * lp->get_objective(k,i));;
        }
        //objectiveVector[k] = (int)round(objectiveVector[k]);
    }

  /*  std::cout << "initial point was ";
    pts.print();
    std::cout << " and current point is ";
    print();
    std::cout << "\n";*/
}

Solution::Solution(BranchingDecisions* bd, MathematicalModel* lp) : dominated(false), cpuCreation(0) {

    // fill solution vector
    int n = lp->get_n();
    variableVector = std::vector<int>(n);
    for (int i = 0; i < n; i++) {
        variableVector[i] = bd->lb[i];
    }

    // fill objective vector
    int p = lp->get_p();
    objectiveVector = std::vector<int>(p);
    for (int k = 0; k < p; k++) {
        for (int i = 0; i < n; i++) {
            objectiveVector[k] += bd->lb[i] * lp->get_objective(k, i);//(int)round(pts.get_preImage(i) * lp->get_objective(k, i));;
        }
    }
}

Solution::Solution(std::vector<double>& y, MathematicalModel* lp) : dominated(false), cpuCreation(0) {

    // fill solution vector
    int n = lp->get_n();
    variableVector = std::vector<int>(n);
    for (int i = 0; i < n; i++) {
        variableVector[i] = y[i];
    }

    // fill objective vector
    int p = lp->get_p();
    objectiveVector = std::vector<int>(p);
    for (int k = 0; k < p; k++) {
        for (int i = 0; i < n; i++) {
            objectiveVector[k] += y[i] * lp->get_objective(k, i);//(int)round(pts.get_preImage(i) * lp->get_objective(k, i));;
        }
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

        // pre-img
        if (variableVector.size() != 0) {
            std::cout << "\n   -> ( ";
            for (int i = 0; i < variableVector.size() - 1; i++) {
                std::cout << variableVector[i] << " , ";
            }
            std::cout << variableVector[variableVector.size() - 1] << " )\n";
        }
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

bool Solution::isFeasible(MathematicalModel* lp) {

    bool feasible = true;
    int j = 0, lhs;

    while (feasible && j < lp->get_m()) {
        
        // compute left-hand side of cte

        lhs = 0;
        for (int i = 0; i < lp->get_n(); i++) {
            lhs += variableVector[i] * lp->get_constraint(j, i);
        }

        // check with the rhs of the constraint

        if (lp->get_signCte(j) == 0) { // >= cte
            if (lhs < lp->get_rhs(j)) feasible = false;
        }
        else if (lp->get_signCte(j) == 1) { // <= cte
            if (lhs > lp->get_rhs(j)) feasible = false;
        }
        else if (lp->get_signCte(j) == 2) { // == cte
            if (lhs != lp->get_rhs(j)) feasible = false;
        }

        j++;
    }

    return feasible;
}



void Solution::setCpu(double t) {
    cpuCreation = t;
}

bool Solution::isEqual(Solution* y) {

    bool equal = true;
    for (int k = 0; k < objectiveVector.size(); k++) {
        if (y->get_objVector(k) != objectiveVector[k]) {
            equal = false;
        }
    }

    return equal;
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

/*! \brief Returns the value of variable x_i.
 *
 * \param i integer. The index of the variable to look at.
 * \return the value of this variable, as a double.
 */
int Solution::get_preImage(int i) {
    return variableVector[i];
}

/*! \brief Checks whether the solution is discarded.
 *
 * \return true if the solution is discarded, false otherwise.
 */
bool Solution::isDiscarded() {
    return dominated;
}


double Solution::getCpu() {
    return cpuCreation;
}