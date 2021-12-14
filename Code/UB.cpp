#include "UB.h"

/* ==========================================================
        Constructors
 ========================================================= */

 /*! \brief Default constructor of an upper bound set.
  *
  * \param lp MathematicalModel. Used to get the number of objectives.
  */
UpperBoundSet::UpperBoundSet(MathematicalModel& lp) : incumbentSet(0), lastGeneratedId(0), lp(&lp), cpuSol() {
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
    /*std::vector<int> v = { 2919, 3342, 1230 };
    Solution sRef = Solution(v);*/

    Solution* newSol = new Solution(P, lp);
    //newSol->print();
    /*if (newSol->isEqual(&sRef)) {
        std::cout << " We found the missing point : ";
        newSol->print();
    }*/
    cpuSol.StopTimer();
    newSol->setCpu(cpuSol.CumulativeTime("sec"));
    cpuSol.StartTimer();
    bool nonDominated = true; //!< true if the new Solution sol is not dominated by any points of the incumbentSet so far
    std::list<Solution*>::iterator y = incumbentSet.begin();
    std::list<Solution*>::iterator nextY = incumbentSet.begin();
    
    //int cpt = 0;
    while (nonDominated && y != incumbentSet.end()) {
        if ((*y)->dominate(*newSol)) {
            nonDominated = false;

            /*if (newSol->isEqual(&sRef)) {
                (*y)->print();
                std::cout << " dominates (newSol) ";
                newSol->print();
            }*/

        }
        y++;
    }

    if (nonDominated) {

        //std::cout << "\n\n new sol: ";
        //newSol->print();

        // search for dominated points
        while (nextY != incumbentSet.end()) {
            y = nextY;
            nextY++;
            if (newSol->dominate(**y)) {
                (*y)->discard();
                //(*y)->print();
                //std::cout << " is discarded\n";
                /*if ((*y)->isEqual(&sRef)) {
                    newSol->print();
                    std::cout << " (newSol) dominates ";
                    (*y)->print();
                }*/
            }
        }

        // update local upper bounds
        /*if (newSol->isEqual(&sRef)) {
            std::cout << " We add the missing point...\n";
        }*/
        incumbentSet.push_back(newSol);
        updateLUB(newSol);

        // delete the dominated points
        nextY = incumbentSet.begin();
        while (nextY != incumbentSet.end()) {
            y = nextY;
            nextY++;
            if ((*y)->isDiscarded()) {
                /*if ((*y)->isEqual(&sRef)) {
                    std::cout << " We delete the missing point...\n";
                }*/
                delete* y;
                incumbentSet.erase(y);
            }
        }
        if (ACTIVATE_NO_GOOD_CUTS) {
            fcm->generateNoGoodConstraint(*newSol);
            pm->generateNoGoodConstraint(*newSol);
        }
    }
    else {
        delete newSol;
    }
}

/*! \brief Update the upper bound set with a new point.
 *
 * \param P Point. The new point added to the upper bound set. Assumed to be integer.
 */
void UpperBoundSet::updateUB(BranchingDecisions* bd) {

    //printLub();
    /*std::vector<int> v = { 2919, 3342, 1230 };
    Solution sRef = Solution(v);*/

    Solution* newSol = new Solution(bd, lp);
    //std::cout << "   -> new sol tested : ";
    //newSol->print();
    cpuSol.StopTimer();
    newSol->setCpu(cpuSol.CumulativeTime("sec"));
    cpuSol.StartTimer();
    bool nonDominated = true; //!< true if the new Solution sol is not dominated by any points of the incumbentSet so far
    std::list<Solution*>::iterator y = incumbentSet.begin();
    std::list<Solution*>::iterator nextY = incumbentSet.begin();

    // test for feasibility
    nonDominated = newSol->isFeasible(lp); // note: nonDominated has the same role as feasible

    //int cpt = 0;
    while (nonDominated && y != incumbentSet.end()) {
        if ((*y)->dominate(*newSol)) {
            nonDominated = false;
            /*if (newSol->isEqual(&sRef)) {
                std::cout << " We found the missing point : ";
                newSol->print();
            }*/
        }
        y++;
    }

    if (nonDominated) {

        //std::cout << "\n\n new sol: ";
        //newSol->print();

        // search for dominated points
        while (nextY != incumbentSet.end()) {
            y = nextY;
            nextY++;
            if (newSol->dominate(**y)) {
                (*y)->discard();
                //(*y)->print();
                //std::cout << " is discarded\n";
                /*if ((*y)->isEqual(&sRef)) {
                    newSol->print();
                    std::cout << " (newSol) dominates ";
                    (*y)->print();
                }*/
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
                /*if ((*y)->isEqual(&sRef)) {
                    std::cout << " We delete the missing point...\n";
                }*/
                delete* y;
                incumbentSet.erase(y);
            }
        }
        if (ACTIVATE_NO_GOOD_CUTS) {
            fcm->generateNoGoodConstraint(*newSol);
            pm->generateNoGoodConstraint(*newSol);
        }
    }
    else {
        delete newSol;
    }
}

/*! \brief Update the upper bound set with a new point.
 *
 * \param P Point. The new point added to the upper bound set. Assumed to be integer.
 */
void UpperBoundSet::updateUB(std::vector<double>& s) {

    bool isInteger = true;

    for (int i = 0; i < s.size(); i++) {
        if (s[i] - floor(s[i]) >= EPS_INT) {
            isInteger = false;
        }
    }

    if (isInteger) {

        Solution* newSol = new Solution(s, lp);
        //std::cout << "   -> new sol tested : ";
        //newSol->print();
        cpuSol.StopTimer();
        newSol->setCpu(cpuSol.CumulativeTime("sec"));
        cpuSol.StartTimer();
        bool nonDominated = true; //!< true if the new Solution sol is not dominated by any points of the incumbentSet so far
        std::list<Solution*>::iterator y = incumbentSet.begin();
        std::list<Solution*>::iterator nextY = incumbentSet.begin();

        // test for feasibility
        nonDominated = newSol->isFeasible(lp); // note: nonDominated has the same role as feasible

        //int cpt = 0;
        while (nonDominated && y != incumbentSet.end()) {
            if ((*y)->dominate(*newSol)) {
                nonDominated = false;
                /*if (newSol->isEqual(&sRef)) {
                    std::cout << " We found the missing point : ";
                    newSol->print();
                }*/
            }
            y++;
        }

        if (nonDominated) {

            //std::cout << "\n\n new sol: ";
            //newSol->print();

            // search for dominated points
            while (nextY != incumbentSet.end()) {
                y = nextY;
                nextY++;
                if (newSol->dominate(**y)) {
                    (*y)->discard();
                    //(*y)->print();
                    //std::cout << " is discarded\n";
                    /*if ((*y)->isEqual(&sRef)) {
                        newSol->print();
                        std::cout << " (newSol) dominates ";
                        (*y)->print();
                    }*/
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
                    /*if ((*y)->isEqual(&sRef)) {
                        std::cout << " We delete the missing point...\n";
                    }*/
                    delete* y;
                    incumbentSet.erase(y);
                }
            }
            if (ACTIVATE_NO_GOOD_CUTS) {
                fcm->generateNoGoodConstraint(*newSol);
                pm->generateNoGoodConstraint(*newSol);
            }
        }
        else {
            delete newSol;
        }
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
            //std::cout << "  -> new lub discarded: ";
            //u->print();
        }
        else {
            //u->print();
            for (int k = 0; k < lp->get_p(); k++) {
                if (u->isDefinedBy(s, k)) {
                    //std::cout << "  -> lub defined by new sol on component " << k << " : ";
                    //u->print();
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
                //std::cout << "     -> on component " << j << ", new lub: ";
                //lubSet.back().print();
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

    // correction check

    /*for (u = lubSet.begin(); u != lubSet.end(); u++) {
        if (u->isRedundant(lubSet)) {

        }
    }*/
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

/*! \brief Check the largest difference between the ws value w, and the weighted sum of a lub in cone s with weights l.
 *
 * \param double w. The weighted sum value.
 * \param vector of double l. The weight vector of the objective functions.
 * \param branchingDecisions bd. Include the slub s.
 */
double UpperBoundSet::getLargestGap(double w, std::vector<double>& l, BranchingDecisions* bd) {

    double largestGap = -1;
    double gap = -1;
    std::list<LocalUpperBound>::iterator u;
    SLUB s = SLUB(bd->slub);

    for (u = lubSet.begin(); u != lubSet.end(); u++) {
        if (s.dominated(*u)) {
            gap = u->getWeightedSum(l) - w;
            if (gap > largestGap) {
                largestGap = gap;
            }
        }
    }

    return -largestGap; // return the opposite because of min heap data structure in node selection
}

/* Test whether the weighted sum with weights l of the lubs included in the OB cone defined in bd is lower or greater
 * than value ws. Return false if at least one of them has a greater weighted sum value, true otherwise.
 *
 * \param BranchingDecisions* bd. The branching decisions that define the cone we are located in in the objective space.
 * \param vector of double l. The weight vector of the weighted sums considered.
 * \param double ws. The value used for comparison of the weighted sums.
 */
bool UpperBoundSet::testWeightedSumValue(BranchingDecisions* bd, std::vector<double>& l, double w) {

    SLUB s = SLUB(bd->slub);
    bool atLeastOneDominated = false;
    double wsValue;
    std::list<LocalUpperBound>::iterator u = lubSet.begin();

    if (lubSet.size() == 1) {
        atLeastOneDominated = true;
    }
    else {
        while (!atLeastOneDominated && u != lubSet.end()) {
            if (s.dominated(*u)) { // if the lub u is located in C(s)
                wsValue = u->getWeightedSum(l);
                //std::cout << wsValue << " vs " << w << "\n";
                if (wsValue >= w) {
                    atLeastOneDominated = true;
                }
            }
            u++;
        }
    }

    return !atLeastOneDominated;
}

/* Get the largest value for a local upper bound with weight l and in the slub defined in bd.
 *
 * \param BranchingDecisions* bd. The branching decisions that define the cone we are located in in the objective space.
 * \param vector of double l. The weight vector of the weighted sums considered.
 */
double UpperBoundSet::getLargestWeightedSumValue(BranchingDecisions* bd, std::vector<double>& l) {

    SLUB s = SLUB(bd->slub);
    double wsValue, wsMax = -1000000;
    std::list<LocalUpperBound>::iterator u;
    //std::cout << "lub size: " << lubSet.size() << "\n";

    if (lubSet.size() == 1) {
        wsMax = 10000000;
    }
    else {
        for (u = lubSet.begin(); u != lubSet.end(); u++) {
            if (s.dominated(*u)) { // if the lub u is located in C(s)
                wsValue = u->getWeightedSum(l);
                //std::cout << wsValue << " vs " << wsMax << "\n";
                if (wsValue >= wsMax) {
                    wsMax = wsValue;
                }
            }
        }
    }

    return wsMax;
}


std::vector<double> UpperBoundSet::getLubWeigth(BranchingDecisions* bd) {

    SLUB s(bd->slub);
    int p = lubSet.begin()->get_dim();
    std::vector<int> idealPts(p, 1000000);
    std::list<LocalUpperBound>::iterator u;
    std::vector<std::list<LocalUpperBound>::iterator> lexiSol(p);
    std::vector<double> w(p);

    for (u = lubSet.begin(); u != lubSet.end(); u++) {
        if (s.dominated(*u)) {
            for (int k = 0; k < p; k++) {
                if (u->get_coordinate(k) < idealPts[k]) {
                    idealPts[k] = u->get_coordinate(k);
                    lexiSol[k] = u;
                }
            }
        }
    }

    // NOT FINISHED !!

    return w;
}


/*! \brief Provides a pointer to the list of local upper bounds.
 *
 * \return a pointer to the list of local upper bounds.
 */
std::list<Solution*>* UpperBoundSet::getYN() {
    return &incumbentSet;
}

void UpperBoundSet::addSol(Solution* s) {
    incumbentSet.push_back(s);
    updateLUB(s);
}

void UpperBoundSet::setCplexModels(FeasibilityCheckModel* ff, ProbingModel* pp) {
    fcm = ff;
    pm = pp;
}



void UpperBoundSet::startTimerUB() {
    cpuSol.StartTimer();
}

void UpperBoundSet::stopTimerUB() {
    cpuSol.StopTimer();
}

std::list<Solution*>* UpperBoundSet::getIncumbentSet() {
    return &incumbentSet;
}