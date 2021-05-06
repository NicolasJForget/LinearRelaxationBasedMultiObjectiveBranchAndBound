#include "SLUB.h"

/* ==========================================================
        Constructors
 ========================================================= */

/*! \brief Constructor of slub given branching decisions.
	 *
	 * This function creates a slub, by assigning the slub field of the BranchingDecision to coord.
	 * \param branchDec BranchingDecision. The branching decision from which we read the slub field.
	 */
SLUB::SLUB(BranchingDecisions& branchDec, Parameters* param) : dim(branchDec.slub.size()) {
    coord = std::vector<int>(dim);
    for (int k = 0; k < dim; k++) {
        coord[k] = branchDec.slub[k];
        //if (param->objectiveBranching != NO_OBJECTIVE_BRANCHING)
            coord[k] -= param->GCD[k];
    }
}

/*! \brief Constructor of slub given a local upper bound.
 *
 * This function creates a slub, by assigning the coordinates of the local upper bound to the coord field.
 * \param lub LocalUpperBound*. A pointer to the local upper bound from which we read the coordinates.
 */
SLUB::SLUB(LocalUpperBound& lub) : dim(lub.get_dim()) {

    coord = std::vector<int>(lub.get_dim());
    for (int k = 0; k < lub.get_dim(); k++) {
        coord[k] = lub.get_coordinate(k);
    }
}

/*! \brief Constructor for an empty slub in dimension p.
 * \param p int. The dimension of the slub.
 */
SLUB::SLUB(int p) : coord(p,INT_MIN) , dim(p) {
}


/* ==========================================================
        Regular Methods
 ========================================================= */

/*! \brief Merges two slub.
 *
 * This function merges this slub with the slub given in parameter. It replaces coord of this slub by the nadir point of
 * this slub and the other.
 * \param slub SLUB. The slub that has to be merged with this one.
 */
void SLUB::merge(SLUB& slub) {

    for (int k = 0; k < dim; k++) {
        if (slub.get_coordinate(k) > coord[k]){
            coord[k] = slub.get_coordinate(k);
        }
    }
}

/*! \brief Merges a slub and a lub.
 *
 * This function merges this slub with the lub given in parameter. It replaces coord of this slub by the nadir point of
 * this slub and the lub.
 * \param lub LocalUpperBound. The lub that has to be merged with this slub.
 */
void SLUB::merge(LocalUpperBound& lub) {

    for (int k = 0; k < dim; k++) {
        if (lub.get_coordinate(k) > coord[k]) {
            coord[k] = lub.get_coordinate(k);
        }
    }
}

/*! \brief Check whether this SLUB is dominated the point pts.
 *
 * \param pts Point*. A pointer to the point used for comparison.
 * \return true if this slub is dominated by the point pts, i.e. pts is located in the cone defined by this SLUB.
 */
bool SLUB::dominated(Point* pts) {

    bool dominated = true;
    if (pts->isDominater() == 1) {
        dominated = true;
    }
    else if (pts->isDominater() == 0) {
        dominated = false;
    }
    else {
        int k = 0;
        while (dominated && k < dim) {
            if (pts->get_objVector(k) > coord[k]) {
                dominated = false;
            }
            k++;
        }
    }

    return dominated;
}

/*! \brief Return the value of a specific coordinate.
 *
 * \param k int. The index of the coordinate.
 * \return the value of the coordinate
 */
int SLUB::get_coordinate(int k) {
    return coord[k];
}

/* \brief Print the slub
 */
void SLUB::print() {
    std::cout << "( ";
    for (int k = 0; k < dim - 1; k++) {
        std::cout << coord[k] << " , ";
    }
    std::cout << coord[dim - 1] << " )\n";
}