#include "BranchAndBound.h"

/* ==========================================================
        Constructors
 ========================================================= */

 /*! \brief Default constructor of a branch-and-bound
  *
  * \param instance string. The path to the file that contains the instance solved by this branch and bound
  */
BranchAndBound::BranchAndBound(std::string instance) : lp(instance), queue(0), P(), U(lp), stat() {}

 /* ==========================================================
         Regular Methods
  ========================================================= */

  /*! \brief Runs the branch-and-bound algorithm with the given parameters.
   *
   * \param lb int. The identifier of the lower bound set.
   * \param nodeSel int. The identifier of the node selection strategy.
   * \param varSel int. The identifier of the variable selection strategy.
   * \param ob int. The identifier of the objective branching strategy.
   */
void BranchAndBound::run(int lb, int nodeSel, int varSel, int ob) {

    try {
        reset();

        Timer cpuTime;
        cpuTime.StartTimer();

        // set parameters
        P.LBset = lb;
        P.nodeSelection = nodeSel;
        P.variableSelection = varSel;
        P.objectiveBranching = ob;
        P.GCD = std::vector<int>(lp.get_p(), 1); // compute with euclidean algorithm for each objective?

        // initialization of the B&B
        queue.push_back(new Node(&lp,&P)); // creation of the root node
        Node* currentNode = nullptr;

        // explore the tree
        while (queue.size() != 0) {
            currentNode = selectNode();
            //currentNode->print();
            currentNode->process(U);
            //U.print();
            //currentNode->showStatus();
            if (!currentNode->isFathomed()) {
                //split
                currentNode->splitOS(&queue);
            }
            currentNode->getStatistics(&stat);
            delete currentNode;
        }

        cpuTime.StopTimer();
        stat.totalTime = cpuTime.ElapsedTime("sec");
    }
    catch (std::string msg) {
        std::cout << msg << std::endl;
        exit(1);
    }
}

/*! \brief Select the next node to be explored in the tree
 *
 * \return a pointer to a Node.
 */
Node* BranchAndBound::selectNode() {

    Node* selectedNode = nullptr;

    if (P.nodeSelection == DEPTH_FIRST) {
        selectedNode = queue.back();
        queue.pop_back();
    }
    else if (P.nodeSelection == BREADTH_FIRST) {
        selectedNode = queue.front();
        queue.pop_front();
    }
    else {
        //std::cout << "  !!! Undefined node selection !!!\n";
        throw std::string("Error: Undefined node selection");
    }

    return selectedNode;
}

/*! \brief Print out YN
 */
void BranchAndBound::printYN() {
    std::cout << "\n\n\n ================================================== \n";
    if (queue.size() >= 1) {
        std::cout << " !!! There are still some unexplored nodes, this is thus an approximation of YN !!!";
    }
    U.print();
}

/*! \brief Reset the branch-and-bound
 *
 * This function reset the branch-and-bound by deleting all its elements except the linear program.
 */
void BranchAndBound::reset() {

    std::list<Node*>::iterator nd;
    for (nd = queue.begin(); nd != queue.end(); nd++) {
        delete *nd;
    }
    queue.clear();

    U.reset();
    P = Parameters();
    stat = Statistics();
}

/*! \brief Print out statistics
 *
 * This function prints various statistics about this branch-and-bound run.
 */
void BranchAndBound::printStatistics() {

    std::cout << "\n\n =====================================================================\n";
    std::cout << "       Node selection | ";
    if (P.nodeSelection == BREADTH_FIRST) std::cout << "Breadth first\n";
    else if (P.nodeSelection == DEPTH_FIRST) std::cout << "Depth first\n";

    std::cout << "   Variable selection | ";
    if (P.variableSelection == FIRST_INDEX) std::cout << "First free index\n";
    else if (P.variableSelection == MOST_OFTEN_FRACTIONAL) std::cout << "Most often fractional\n";

    std::cout << "  Objective Branching | ";
    if (P.objectiveBranching == NO_OBJECTIVE_BRANCHING) std::cout << "None\n";
    else if (P.objectiveBranching == FULL_OBJECTIVE_BRANCHING) std::cout << "Full\n";
    else if (P.objectiveBranching == CONE_OBJECTIVE_BRANCHING) std::cout << "Cone\n";

    std::cout << "      Lower Bound Set | ";
    if (P.LBset == LP_RELAX) std::cout << "LP relax (no warmstart)\n";
    else if (P.LBset == WARMSTARTED_LP_RELAX) std::cout << "Warmstarted LP relax\n";

    std::cout << "\n Non-dominated points found: " << U.getNbNonDominatedPoints() << "\n";

    double cputOther = stat.totalTime - stat.timeLBComputation - stat.timeUpdateUB - stat.timeDominanceTest;
    std::cout << "\n Elapsed time: " << stat.totalTime << " sec\n";
    std::cout << "    -> LB computation: " << stat.timeLBComputation << " sec ( " << 100 * stat.timeLBComputation / stat.totalTime << " % )\n";
    std::cout << "    -> Dominance test: " << stat.timeDominanceTest << " sec ( " << 100 * stat.timeDominanceTest / stat.totalTime << " % )\n";
    std::cout << "    -> Update UB: " << stat.timeUpdateUB << " sec ( " << 100 * stat.timeUpdateUB / stat.totalTime << " % )\n";
    std::cout << "    -> others: " << cputOther << " sec ( " << 100 * cputOther / stat.totalTime << " % )\n";

    int nbFathomed = stat.nbFathomedInfeasibility + stat.nbFathomedDominance + stat.nbFathomedOptimality;
    std::cout << "\n Node explored: " << stat.nbNodes << "\n";
}