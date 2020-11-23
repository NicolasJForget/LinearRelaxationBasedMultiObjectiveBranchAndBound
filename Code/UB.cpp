#include "UB.h"

/* ==========================================================
        Constructors
 ========================================================= */

 /*! \brief Default constructor of an upper bound set.
  *
  * \param lp MathematicalModel. Used to get the number of objectives.
  */
UpperBoundSet::UpperBoundSet(MathematicalModel& lp) : incumbentSet(0), lastGeneratedId(0), lp(&lp) {
    lubSet = std::list<LocalUpperBound>(0);
    lubSet.push_back(LocalUpperBound(lp));
}

 /* ==========================================================
         Regular Methods
  ========================================================= */

  /*! \brief Generate a new unique id for a local upper bound.
   */
int UpperBoundSet::newId() {
    ++lastGeneratedId;
    return lastGeneratedId;
}

/*! \brief Update the upper bound set with a new point.
 *
 * \param P Point. The new point added to the upper bound set. Assumed to be integer.
 */
void UpperBoundSet::updateUB(Point& P) {

    //printLub();

    Solution* newSol = new Solution(P);
    bool nonDominated = true; //!< true if the new Solution sol is not dominated by any points of the incumbentSet so far
    std::list<Solution*>::iterator y = incumbentSet.begin();
    std::list<Solution*>::iterator nextY = incumbentSet.begin();
    
    //int cpt = 0;
    while (nonDominated && y != incumbentSet.end()) {
        if ((*y)->dominate(*newSol)) {
            nonDominated = false;
        }
        y++;
    }

    if (nonDominated) {

        std::cout << "\n\n new sol: ";
        newSol->print();

        // search for dominated points
        while (nextY != incumbentSet.end()) {
            y = nextY;
            nextY++;
            if (newSol->dominate(**y)) {
                (*y)->discard();
                (*y)->print();
                std::cout << " is discarded\n";
            }
        }

        // update local upper bounds
        incumbentSet.push_back(newSol);
        updateLUB(newSol);

        // delete the dominated points
        nextY = incumbentSet.begin();
        while (nextY != incumbentSet.end()) {
            y = nextY;
            nextY++;
            if ((*y)->isDiscarded()) {
                delete* y;
                incumbentSet.erase(y);
            }
        }
    }
    else {
        delete newSol;
    }
}

/*! \brief Prints the upper bound set.
 *
 * This function prints the objective vector of each solution of the upper bound set.
 */
void UpperBoundSet::print() {

    std::cout << "\n\n --- UB set ---\n";
    std::cout << "   -> " << incumbentSet.size() << " non-dominated points\n";
    std::list<Solution*>::iterator y;
    for (y = incumbentSet.begin(); y != incumbentSet.end(); y++) {
        (*y)->print();
    }
}

/*! \brief Prints the set of local upper bounds.
 *
 * This function prints the set of local upper bounds, by showing their coordindates in the objective space
 */
void UpperBoundSet::printLub() {
    std::cout << "\n\n ============= Local Upper Bounds =============\n";
    std::cout << "   -> " << lubSet.size() << " local upper bounds, corresponding to " << incumbentSet.size() << " non-dominated points\n\n";
    std::list<LocalUpperBound>::iterator u;
    for (u = lubSet.begin(); u != lubSet.end(); u++) {
        u->print();
    }
}

/*! \brief Update the set of local upper bounds
 *
 * This function updates the set of local upper bounds, given the new Solution s added to the upper bound set.
 * \param s Solution. The new solution added to the upper bound set.
 */
void UpperBoundSet::updateLUB(Solution* s) {

    std::list<LocalUpperBound>::iterator u;
    std::list<LocalUpperBound>::iterator uNext;
    std::list<LocalUpperBound*> A; // list of lub that should be discarded because strictly dominated by s
    std::list<LocalUpperBound*>::iterator uDiscarded;
    int zMax;

    // Update defining lists, in case we detect this new point defines already existing local upper bounds

    for (u = lubSet.begin(); u != lubSet.end(); u++) {
        if (u->isStrictlyDominated(s)) {
            u->becomesDiscarded();
            A.push_back(&(*u));
            std::cout << "  -> new lub discarded: ";
            u->print();
        }
        else {
            for (int k = 0; k < lp->get_p(); k++) {
                if (u->isDefinedBy(s, k)) {
                    std::cout << "  -> lub defined by new sol on component " << k << " : ";
                    u->print();
                    u->filterDefiningComponent(k);
                    u->addDefiningComponent(s, k);
                }
            }
        }
    }

    // Compute the new local upper bounds that have to be computed

    for (uDiscarded = A.begin(); uDiscarded != A.end(); uDiscarded++) {
        //std::cout << " ===> new lubs with ";
        //(*uDiscarded)->print();
        //(*uDiscarded)->print();
        for (int j = 0; j < lp->get_p(); j++) {
            zMax = (*uDiscarded)->computeCriticalValue(j);
            //std::cout << "    zMax = " << zMax << std::endl;
            if (s->get_objVector(j) > zMax) {
                lubSet.push_back(LocalUpperBound(s, **uDiscarded, j, genNewId()));
                //for (int i = 0; i < lp->get_p(); i++) std::cout << lubSet.back().get_definingPoints(i)->size();
                std::cout << "     -> on component " << j << ", new lub: ";
                lubSet.back().print();
            }
        }
    }

    // Discarded the old local upper bounds

    uNext = lubSet.begin();
    do
    {
        u = uNext;
        uNext++;
        if (u->isDiscarded()) {
            //std::cout << "size UB is " << lubSet.size() << std::endl;
            lubSet.erase(u);
            //std::cout << "size UB is " << lubSet.size() << std::endl;
        }
    } while (uNext != lubSet.end());
}

/*! \brief Generate a new unique id for a new local upper bound
 *
 * \return the id as an int.
 */
int UpperBoundSet::genNewId() {
    ++lastGeneratedId;
    return lastGeneratedId;
}

/*! \brief Provides a pointer to the list of local upper bounds.
 *
 * This function returns a pointer to the list of local upper bounds of this upper bound set.
 * \return a pointer to the list of local upper bounds.
 */
std::list<LocalUpperBound>* UpperBoundSet::getLubs() {
    return &lubSet;
}

/*! \brief Returns the number of non-dominated points
 *
 * This function returns the number of non-dominated points, by returning the size of incumbentSet.
 * \return number of non-dominated points, as an int.
 */
int UpperBoundSet::getNbNonDominatedPoints() {
    return incumbentSet.size();
}

/*! \brief Reset the UpperBoundSet.
 *
 * This function reset the upper bound set by deleting all its components.
 */
void UpperBoundSet::reset() {

    lubSet.clear();
    lubSet.push_back(LocalUpperBound(*lp));

    std::list<Solution*>::iterator sol;
    for (sol = incumbentSet.begin(); sol != incumbentSet.end(); sol++) {
        delete *sol;
    }
    incumbentSet.clear();

    lastGeneratedId = 0;
}