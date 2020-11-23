#include "Node.h"

/* ==========================================================
        Constructors
 ========================================================= */

/*! \brief Default constructor of a node.
 */
Node::Node() : P(nullptr), param(nullptr), LB(), score(-1), branchingDecision(), splittingIndex(-1), timeLBComputation(), timeDominanceTest(), timeUpdateUB() {}; //, status(UNSOLVED)

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
Node::Node(MathematicalModel* lp, Parameters* par) : P(lp), param(par), score(-1), branchingDecision(lp->get_n(),lp->get_p()), splittingIndex(-1), timeLBComputation(), timeDominanceTest(), timeUpdateUB() { //, status(UNSOLVED)

    // lb & ub variables & ob
    for (int i = 0; i < lp->get_n(); i++) {
        branchingDecision.lb[i] = 0; // change to lb defined in file
        branchingDecision.ub[i] = 1; // change to lb defined in file
    }

    for (int k = 0; k < lp->get_p(); k++) {
        branchingDecision.slub[k] = INT_MAX;
    }

    // define appropriate lb -> switch lb case ??
    if (par->LBset == LP_RELAX || par->LBset == WARMSTARTED_LP_RELAX) {
        LB = new LinearRelaxation(lp,&branchingDecision);
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
Node::Node(Node& nd, int index, int bound, int val) : P(nd.P), param(nd.param), score(-1), splittingIndex(index), timeLBComputation(), timeDominanceTest(), timeUpdateUB() {

    branchingDecision = BranchingDecisions(nd.branchingDecision);
    if (bound == IS_LB) {
        branchingDecision.lb[index] = val;
    }
    else if (bound == IS_UB) {
        branchingDecision.ub[index] = val;
    }
    else {
        throw std::string("Error: constraint sign for branching decision is not valid\n");
    }

    if (nd.param->LBset == LP_RELAX) {
        LinearRelaxation* cst = dynamic_cast<LinearRelaxation*>(nd.LB);
        LB = new LinearRelaxation(nd.P,cst->getWeightedSumModel(),cst->getFeasibilityCheckModel(),cst->getDualBensonModel(),cst->getFurthestFeasiblePointModel(),&branchingDecision);
    }
    else if (nd.param->LBset == WARMSTARTED_LP_RELAX) {
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
void Node::process(UpperBoundSet& U) {

    timeLBComputation.StartTimer();
    LB->applyBranchingDecisions();
    LB->compute();
    timeLBComputation.StopTimer();

    if (LB->getStatus() == UNSOLVED) {
        throw std::string("Error: lower bound set unsolved");
    }
    else if (LB->getStatus() != INFEASIBLE) {

        timeUpdateUB.StartTimer();
        LB->gatherIntegerSolutions(U);
        timeUpdateUB.StopTimer();

        if (LB->getStatus() != OPTIMAL) {
            timeDominanceTest.StartTimer();
            LB->applyDominanceTest(U,param);
            timeDominanceTest.StopTimer();
        }
    }

}

/*! \brief Check whether the node is fathomed.
 *
 * \return true is the node is fathomed, false otherwise.
 */
bool Node::isFathomed() {
    return (LB->getStatus() == INFEASIBLE || LB->getStatus() == OPTIMAL || LB->getStatus() == DOMINATED);
}

/*! \brief Split the problem in the objective space
 *
 * \param Q pointer to a list of pointers of Nodes. The list in which the new nodes created are stored.
 */
void Node::splitOS(std::list<Node*>* Q) {

    if (param->objectiveBranching == NO_OBJECTIVE_BRANCHING) {
        splitVS(Q);
    }
    else {
        throw std::string("Error: objective branching parameter not supported\n");
    }

}

/*! \brief Split the problem in the variable space
 *
 * \param Q pointer to a list of pointers of Nodes. The list in which the new nodes created are stored.
 */
void Node::splitVS(std::list<Node*>* Q) {

    std::vector<int> index(3,0); // 0: index splitted, 1: -ub 1st node, as a negative number, 2: lb 2nd node
    if (param->variableSelection == FIRST_INDEX) {
        Q->push_back(new Node(*this, splittingIndex + 1, IS_UB, 0));
        Q->push_back(new Node(*this, splittingIndex + 1, IS_LB, 1));
    }
    else {
        throw std::string("Error: variable selection parameter not supported\n");
    }
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
    std::cout << "\n";
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
void Node::getStatistics(Statistics* stat) {

    stat->nbNodes++;

    // status
    if (LB->getStatus() == INFEASIBLE) {
        stat->nbFathomedInfeasibility++;
    }
    else if (LB->getStatus() == DOMINATED) {
        stat->nbFathomedDominance++;
    }
    else if (LB->getStatus() == OPTIMAL) {
        stat->nbFathomedOptimality++;
    }

    // timers
    stat->timeDominanceTest += timeDominanceTest.ElapsedTime("sec");
    stat->timeLBComputation += timeLBComputation.ElapsedTime("sec");
    stat->timeUpdateUB += timeUpdateUB.ElapsedTime("sec");
}