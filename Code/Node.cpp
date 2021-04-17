#include "Node.h"

/* ==========================================================
        Constructors
 ========================================================= */

/*! \brief Default constructor of a node.
 */
Node::Node() : P(nullptr), param(nullptr), LB(), score(-1), branchingDecision(), splittingIndex(-1), ndLub(0), timeLBComputation(), timeDominanceTest(), timeUpdateUB(), stat(nullptr), depth(0) {}; //, status(UNSOLVED)

/* \brief Destructor of the node.
 *
 */
Node::~Node() {
    if (param->LBset == LP_RELAX || param->LBset == WARMSTARTED_LP_RELAX) {
        delete dynamic_cast<LinearRelaxation*>(LB);
    }
    else {
        std::cout << "Warning: Trying to destroy an undefined lower bound set.";
        //throw std::string("Error: Trying to destroy an undefined lower bound set.");
    }
}

/*! \brief Default constructor of the root node.
 *
 * \param lp MathematicalModel*. A pointer to the initial problem.
 */
Node::Node(MathematicalModel* lp, Parameters* par, Statistics* stat) : P(lp), param(par), score(-1), branchingDecision(lp->get_n(),lp->get_p()), splittingIndex(-1), ndLub(0), timeLBComputation(), timeDominanceTest(), timeUpdateUB(), stat(stat), depth(0) { //, status(UNSOLVED)

    // lb & ub variables & ob
    if (lp->asBoundsGiven()) {
        for (int i = 0; i < lp->get_n(); i++) {
            branchingDecision.lb[i] = lp->getLb(i); // change to lb defined in file
            branchingDecision.ub[i] = lp->getUb(i); // change to lb defined in file
        }
    }
    else {
        for (int i = 0; i < lp->get_n(); i++) {
            branchingDecision.lb[i] = 0; // change to lb defined in file
            branchingDecision.ub[i] = 1; // change to lb defined in file
        }
    }

    for (int k = 0; k < lp->get_p(); k++) {
        for (int i = 0; i < lp->get_n(); i++) {
            if (lp->get_objective(k, i) >= 0) {
                branchingDecision.slub[k] += lp->get_objective(k, i) * lp->getUb(i);
            }
        }
        branchingDecision.slub[k] += 2;
    }
    //for (int k = 0; k < lp->get_p(); k++) {
    //    branchingDecision.slub[k] = INT_MAX;
    //}

    // define appropriate lb -> switch lb case ??
    if (par->LBset == LP_RELAX || par->LBset == WARMSTARTED_LP_RELAX) {
        LB = new LinearRelaxation(lp,&branchingDecision,stat);
    }
    else {
        //std::cout << "  !!! Undefined lower bound set !!!\n";
        throw std::string("Error: Undefined lower bound set");
    }
}

/*! \brief Creates a node given a splitting index and a bound.
 *
 * \param nd Node*. A pointer to the node it is created from.
 * \param index int. The index of the variable to split.
 * \param bound int. States the sign of the new constraint: are we adding an upper or a lower bound value on a variable?
 * \param val int. The value of the bound added on the variable, i.e. the right-hand side of the new branching constraint.
 */
Node::Node(Node& nd, int index, int bound, int val, SLUB& slub) : P(nd.P), param(nd.param), score(-1), splittingIndex(index), timeLBComputation(), timeDominanceTest(), timeUpdateUB(), ndLub(nd.ndLub), stat(nd.stat), depth(nd.depth + 1) {

    branchingDecision = BranchingDecisions(nd.branchingDecision);
    branchingDecision.depth = depth;
    branchingDecision.lastSplittedIndex = index;
    if (bound == IS_LB) {
        branchingDecision.lb[index] = val;
    }
    else if (bound == IS_UB) {
        branchingDecision.ub[index] = val;
    }
    else {
        throw std::string("Error: constraint sign for branching decision is not valid\n");
    }
    for (int k = 0; k < P->get_p(); k++) {
        branchingDecision.slub[k] = slub.get_coordinate(k);
    }

    if (nd.param->LBset == LP_RELAX) {
        LinearRelaxation* cst = dynamic_cast<LinearRelaxation*>(nd.LB);
        LB = new LinearRelaxation(nd.P,cst->getWeightedSumModel(),cst->getFeasibilityCheckModel(),cst->getDualBensonModel(),cst->getFurthestFeasiblePointModel(),&branchingDecision,stat);
    }
    else if (nd.param->LBset == WARMSTARTED_LP_RELAX) {
        if (depth % CORRECTION_WARMSTART == 0) {
            LinearRelaxation* cst = dynamic_cast<LinearRelaxation*>(nd.LB);
            LB = new LinearRelaxation(nd.P, cst->getWeightedSumModel(), cst->getFeasibilityCheckModel(), cst->getDualBensonModel(), cst->getFurthestFeasiblePointModel(), &branchingDecision, stat);
        }
        else
            LB = new LinearRelaxation(dynamic_cast<LinearRelaxation*>(nd.LB),&branchingDecision);
    }
    else {
        throw std::string("Error: Lower bound set not supported for node splitting\n");
    }
}

/* ==========================================================
        Regular Methods
 ========================================================= */

 /*! \brief Process the node
  */
void Node::process(UpperBoundSet& U, int iteration) {

    //print();
    
    stat->timeLBComputation.StartTimer();
    LB->setIteration(iteration);
    LB->applyBranchingDecisions();
    //std::cout << "nik ";
    LB->compute();
    //std::cout << "patamer\n";
    stat->timeLBComputation.StopTimer();

    if (LB->getStatus() == UNSOLVED) {
        throw std::string("Error: lower bound set unsolved");
    }
    else if (LB->getStatus() != INFEASIBLE) {

        stat->timeUpdateUB.StartTimer();
        LB->gatherIntegerSolutions(U);
        stat->timeUpdateUB.StopTimer();

        if (LB->getStatus() != OPTIMAL) {
            stat->timeDominanceTest.StartTimer();
            LB->applyDominanceTest(U,param,ndLub);
            stat->timeDominanceTest.StopTimer();
        }
    }

}

/*! \brief Check whether the node is fathomed.
 *
 * \return true is the node is fathomed, false otherwise.
 */
bool Node::isFathomed() {

    bool fathomed = true;

    if (LB->getStatus() == INFEASIBLE) stat->nbFathomedInfeasibility++;
    else if (LB->getStatus() == OPTIMAL) stat->nbFathomedOptimality++;
    else if (LB->getStatus() == DOMINATED) stat->nbFathomedDominance++;
    else fathomed = false;

    return fathomed;//(LB->getStatus() == INFEASIBLE || LB->getStatus() == OPTIMAL || LB->getStatus() == DOMINATED);
}

/*! \brief Split the problem in the objective space
 *
 * \param Q pointer to a list of pointers of Nodes. The list in which the new nodes created are stored.
 */
void Node::splitOS(std::list<Node*>* Q, UpperBoundSet* U, int iteration) {

    if (param->objectiveBranching == NO_OBJECTIVE_BRANCHING) {
        stat->timeComputeOB.StartTimer();
        SLUB slub = SLUB(branchingDecision, param);
        stat->timeComputeOB.StopTimer();
        splitVS(Q, slub, iteration);
    }
    else if (param->objectiveBranching == CONE_OBJECTIVE_BRANCHING) {

        stat->timeComputeOB.StartTimer();
        SLUB slub = SLUB(P->get_p());
        //slub.print();
        std::list<LocalUpperBound>* NU = U->getLubs();
        std::list<LocalUpperBound>::iterator u;
        std::list<int>::iterator nextId = ndLub.begin();
        bool dominated;
        // search for dominated lubs to compute OB. Note that all non-existing lubs have been deleted beforehand
        // in the dominance test, and thus there is no need to take care of this case.
        for (u = NU->begin(); u != NU->end(); u++) {
            dominated = true;
            if (nextId != ndLub.end() && *nextId == u->get_id()) {
                dominated = false;
                nextId++;
            }
            if (dominated) {
                slub.merge(*u);
            }
            //slub.print();
        }
        stat->timeComputeOB.StopTimer();
        splitVS(Q, slub, iteration);
    }
    else if (param->objectiveBranching == FULL_OBJECTIVE_BRANCHING) {

        stat->timeComputeOB.StartTimer();
        // building the initial set of SLUBs S.
        std::list<SLUB*> S(0);
        std::list<LocalUpperBound>* lubs = U->getLubs();
        std::list<LocalUpperBound>::iterator u;
        std::list<int>::iterator nextId = ndLub.begin();
        for (u = lubs->begin(); u != lubs->end(); u++) { // take only non-dominated ones !!!
            if (nextId != ndLub.end() && *nextId == u->get_id()) {
                nextId++;
            }
            else {
                S.push_back(new SLUB(*u));
            }
        }

        // merging the SLUBs if necessary
        std::vector<int> intersectionPoint(P->get_p(),0);
        std::list<SLUB*>::iterator s1 = S.begin(), s2, sTempo;
        bool aMergingOperationOccured = true;
        bool jobFinished = false;
        while (aMergingOperationOccured) { // it can be improved, e.g. by taking care of a list of merged slubs and checking those only
            s1 = S.begin();
            aMergingOperationOccured = false;
            while (s1 != S.end()) { // !jobFinished
                s2 = s1;
                s2++;
                while (s2 != S.end()) {
                    for (int k = 0; k < P->get_p(); k++) {
                        intersectionPoint[k] = min((*s1)->get_coordinate(k), (*s2)->get_coordinate(k));
                    }
                    if (LB->dominates(intersectionPoint)) {
                        (*s1)->merge(**s2);
                        delete *s2;
                        S.erase(s2);
                        s2 = s1;
                        aMergingOperationOccured = true;
                    }
                    s2++;
                }
                s1++;
            }
        }

        stat->timeComputeOB.StopTimer();
        //std::cout << " sub-pb created: " << S.size() << std::endl;
        for (s1 = S.begin(); s1 != S.end(); s1++) {
            //(*s1)->print();
            splitVS(Q, **s1, iteration);
        }

    }
    else {
        throw std::string("Error: objective branching parameter not supported\n");
    }

}

/*! \brief Split the problem in the variable space
 *
 * \param Q pointer to a list of pointers of Nodes. The list in which the new nodes created are stored.
 */
void Node::splitVS(std::list<Node*>* Q, SLUB& slub, int iteration) {

    stat->timeVarSel.StartTimer();
    //std::vector<int> index(3,0); // 0: index splitted, 1: -ub 1st node, as a negative number, 2: lb 2nd node
    if (param->variableSelection == FIRST_INDEX) {
        if (P->isBinary()) {
            Q->push_back(new Node(*this, splittingIndex + 1, IS_UB, 0, slub));
            Q->push_back(new Node(*this, splittingIndex + 1, IS_LB, 1, slub));
        }
        else {
            throw std::string("Error: integer programs not supported by the current variable selection.");
        }
    }
    else if (param->variableSelection == MOST_OFTEN_FRACTIONAL) {
        if (param->LBset == LP_RELAX || param->LBset == WARMSTARTED_LP_RELAX) {
            int index = dynamic_cast<LinearRelaxation*>(LB)->computeMostOftenFractionalIndex(slub);
            if (P->isBinary()) {
                Q->push_back(new Node(*this, index, IS_UB, 0, slub));
                Q->push_back(new Node(*this, index, IS_LB, 1, slub));
            }
            else {
                int splittingValue = dynamic_cast<LinearRelaxation*>(LB)->computeMedianSplittingValue(slub, index);
                Q->push_back(new Node(*this, index, IS_UB, splittingValue, slub));
                Q->push_back(new Node(*this, index, IS_LB, splittingValue + 1, slub));
                //if (iteration == 1436)
                    //std::cout << " split value is " << splittingValue << " on variable " << splittingIndex << "\n";
            }
        }
        else {
            throw std::string("Error: variable selection parameter not supported with this LB set\n");
        }
    }
    else {
        throw std::string("Error: variable selection parameter not supported\n");
    }
    stat->timeVarSel.StopTimer();
}

/*! \brief Prints the node
 */
void Node::print() {

    std::cout << "\n\n\n ============= New node ===============\n\n";

    for (int i = 0; i < P->get_n(); i++) {
        if (branchingDecision.lb[i] >= 1) { //change for integer var + update lb & ub only if relevant when branching?
            std::cout << " x[" << i << "] >= " << branchingDecision.lb[i] << std::endl;
        }
        if (branchingDecision.ub[i] <= 0) { //change for integer var + update lb & ub only if relevant when branching?
            std::cout << " x[" << i << "] <= " << branchingDecision.ub[i] << std::endl;
        }
    }
    std::cout << "\n   -> last split index : " << branchingDecision.lastSplittedIndex;
    std::cout << "\n\n";
    for (int k = 0; k < P->get_p(); k++) {
        std::cout << " objective " << k << " : " << branchingDecision.slub[k] << "\n";
    }
}

/*! \brief Prints the status of the node.
 *
 * This function prints the status of the node, by reading the status of its lower bound set.
 */
void Node::showStatus() {
    if (LB->getStatus() == SOLVED) {
        std::cout << "\n   -> The node is split";
    }
    else if (LB->getStatus() == INFEASIBLE) {
        std::cout << "\n   -> The node is fathomed by infeasibility";
    }
    else if (LB->getStatus() == OPTIMAL) {
        std::cout << "\n   -> The node is fathomed by optimality";
    }
    else if (LB->getStatus() == DOMINATED) {
        std::cout << "\n   -> The node is fathomed by dominance";
    }
    else {
        std::cout << "\n   -> The node is fathomed by not solved... why ??";
    }
}

/*! \brief Rreturn the status of the node
 *
 * \return status of the node, as an integer
 */
int Node::getStatus() {
    return LB->getStatus();
}

/*! \brief Write statistics of this node
 *
 * \param stat Statistics*. Pointer to the struct where statistics are recorded.
 */
//void Node::getStatistics(Statistics* stat) {
//
//    stat->nbNodes++;
//
//    // status
//    if (LB->getStatus() == INFEASIBLE) {
//        stat->nbFathomedInfeasibility++;
//    }
//    else if (LB->getStatus() == DOMINATED) {
//        stat->nbFathomedDominance++;
//    }
//    else if (LB->getStatus() == OPTIMAL) {
//        stat->nbFathomedOptimality++;
//    }
//
//    // timers
//    stat->timeDominanceTest += timeDominanceTest.ElapsedTime("sec");
//    stat->timeLBComputation += timeLBComputation.ElapsedTime("sec");
//    stat->timeUpdateUB += timeUpdateUB.ElapsedTime("sec");
//
//    // LP relax
//    if (param->LBset == LP_RELAX || param->LBset == WARMSTARTED_LP_RELAX) {
//        dynamic_cast<LinearRelaxation*>(LB)->getStatistics(stat);
//    }
//}

/*! \brief Rreturn the split index.
 * \return the index, as an integer
 */
int Node::getSplittingIndex() {
    return splittingIndex;
}

/*! \brief Rreturn the depth of the node.
 *
 * \return the depth, as an integer
 */
int Node::getDepth() {
    return depth;
}

/*! \brief Print the LB of the node.
 */
void Node::showLB() {
    LB->print();
}

bool Node::isOurCulprit() {

    if (
        branchingDecision.lb[5] == 2 &&
        branchingDecision.lb[8] == 1 &&
        branchingDecision.lb[9] == 1 &&
        branchingDecision.ub[1] == 0 &&
        branchingDecision.ub[2] == 0 &&
        branchingDecision.ub[3] == 0 &&
        branchingDecision.ub[4] == 2 &&
        branchingDecision.ub[5] == 2 &&
        branchingDecision.ub[6] == 0 &&
        branchingDecision.ub[7] == 0 &&
        branchingDecision.ub[8] == 1 &&
        branchingDecision.ub[9] == 1
        )
        return true;
    else
        return false;

}