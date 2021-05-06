#include "BranchAndBound.h"

/* ==========================================================
        Constructors
 ========================================================= */

 /*! \brief Default constructor of a branch-and-bound
  *
  * \param instance string. The path to the file that contains the instance solved by this branch and bound
  */
BranchAndBound::BranchAndBound(std::string instance) : lp(instance), queue(0), P(), U(lp), stat() {

    inst = instance;
    inst.erase(0, 20);

    /*std::stringstream ss(instance);
    std::vector<std::string> result;
    while (ss.good())
    {
        std::string substr;
        getline(ss, substr, '/');
        result.push_back(substr);
    }
    inst = result[result.size() - 1];*/
}

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
void BranchAndBound::run(int lb, int nodeSel, int varSel, int ob, int valBranch, int timeout) {

    try {
        reset();
        int iteration = 0;

        //lp.printObjective();
        //lp.printConstraints();

        Timer timeLimit;
        //cpuTime.StartTimer();
        stat.totalTime.StartTimer();

        // set parameters
        P.LBset = lb;
        P.nodeSelection = nodeSel;
        P.variableSelection = varSel;
        P.objectiveBranching = ob;
        P.branchingValueSelection = valBranch;
        P.timeOut = timeout;
        P.GCD = std::vector<int>(lp.get_p(), 1); // compute with euclidean algorithm for each objective?
        //if (lp.get_p() == 3)
        //    P.objectiveBranching = CONE_OBJECTIVE_BRANCHING;
        //else
        //    P.objectiveBranching = NO_OBJECTIVE_BRANCHING;

        // initialization of the B&B
        queue.push_back(new Node(&lp,&P,&stat)); // creation of the root node
        Node* currentNode = nullptr;

        std::string go;
        // explore the tree
        while (queue.size() != 0 && timeLimit.CumulativeTime("sec") < P.timeOut) { // && iteration < 64    4 // TIME_OUT
            timeLimit.StartTimer();
            stat.timeNodeSel.StartTimer();
            currentNode = selectNode();
            stat.timeNodeSel.StopTimer();
            stat.nbNodes++;

            if (iteration % 1000 == 0) { // 8771 // int(DEBUG_IT)
                /*std::cout << "\n\n Start next iteration...";
                std::cin >> go;*/
                //currentNode->print();
                std::cout << "\niteration: " << iteration << "\n";
                //currentNode->showLB();
                //currentNode->showStatus();
                //if (iteration == 94) throw std::string("Debug\n");
            }
            currentNode->process(U,iteration);
            //currentNode->showLB();
            /*if (currentNode->isOurCulprit()) {
                currentNode->print();
                currentNode->showLB();
                std::cout << "\niteration: " << iteration << "\n";
                std::cout << " we found him\n";
            }*/
            
            /*if (iteration == 4151) {
                currentNode->print();
                currentNode->showLB();
                currentNode->showStatus();
                std::cout << "stop";
            }*/
            if (!currentNode->isFathomed()) {
                //split
                currentNode->splitOS(&queue, &U, iteration);
            }
            else {
                if (currentNode->getDepth() >= stat.maxDepth) stat.maxDepth = currentNode->getDepth();
                if (currentNode->getDepth() <= stat.minDepth) stat.minDepth = currentNode->getDepth();
                stat.avgDepth += currentNode->getDepth();
            }

            //currentNode->showLB();
            //currentNode->showStatus();
            //currentNode->getStatistics(&stat);
            delete currentNode;
            iteration++;
            timeLimit.StopTimer();
        }

        //cpuTime.StopTimer();
        //stat.totalTime = cpuTime.ElapsedTime("sec");
        stat.totalTime.StopTimer();
        if (queue.size() >= 1) stat.solved = 0;
        else stat.solved = 1;
        stat.avgDepth /= (stat.nbFathomedDominance + stat.nbFathomedInfeasibility + stat.nbFathomedOptimality);
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

/*! \brief Write the statistic of this run to the stat file (results.txt / .csv)
 */
void BranchAndBound::writeStatistics() {

    std::ofstream file;
    //file.open("C:/Users/au643334/source/repos/LinearRelaxationBasedMultiObjectiveBranchAndBound/Code/instances/expeData/results.txt", std::ios_base::app); // append instead of overwrite
    file.open("C:/Users/au643334/Desktop/expeProject2/results.txt", std::ios_base::app); // append instead of overwrite
    
    // instance and parameters
    file << inst << ","; // instance
    file << lp.get_p() << ","; // p
    file << lp.get_n() << ","; // n
    if (P.LBset == LP_RELAX) file << "LP,"; // LB set
    else if (P.LBset == WARMSTARTED_LP_RELAX) file << "WLP,";
    else file << ",";
    if (P.nodeSelection == DEPTH_FIRST) file << "depth,"; // node sel
    else if (P.nodeSelection == BREADTH_FIRST) file << "breadth,";
    else file << ",";
    if (P.objectiveBranching == NO_OBJECTIVE_BRANCHING) file << "noOB,"; // OB
    else if (P.objectiveBranching == FULL_OBJECTIVE_BRANCHING) file << "fullOB,";
    else if (P.objectiveBranching == CONE_OBJECTIVE_BRANCHING) file << "coneOB,";
    else file << ",";
    if (P.branchingValueSelection == MEDIAN) file << "med,";
    else if (P.branchingValueSelection == MOST_OFTEN_FRACTIONAL_VALUE) file << "mofv,";
    else file << ",";

    // performance and YN
    if (stat.solved == 1) file << "1,"; // solved
    else if (stat.solved == 0) file << "0,";
    else file << ",";
    file << U.getNbNonDominatedPoints() << ","; // YN
    file << stat.totalTime.CumulativeTime("sec") << ","; // cpuTotal
    
    // LB performance profile
    file << stat.lpSolved << ","; // nbLpSolved
    file << stat.timeInitialization.CumulativeTime("sec") << ","; // cpuInitialization
    file << stat.timeFeasibilityCheck.CumulativeTime("sec") << ","; // cpuCplex
    file << stat.timeUpdatePolyhedron.CumulativeTime("sec") << ","; // cpuUpdatePolyhedron

    // Tree structure
    file << stat.nbNodes << ","; // nbNodes
    file << stat.minDepth << ","; // minDepth
    file << stat.maxDepth << ","; // maxDepth
    file << stat.avgDepth << ","; // avgDepth
    file << stat.nbFathomedInfeasibility << ","; // nbFathomedInfeas
    file << stat.nbFathomedOptimality << ","; // nbFathomedOptimality
    file << stat.nbFathomedDominance << ","; // nbFathomedDominance

    // BB performance profile
    file << stat.timeDominanceTest.CumulativeTime("sec") << ",";
    file << stat.timeUpdateUB.CumulativeTime("sec") << ",";
    file << stat.timeLBComputation.CumulativeTime("sec") << ",";
    file << stat.timeComputeOB.CumulativeTime("sec") << ",";
    file << stat.timeVarSel.CumulativeTime("sec") << ",";
    file << stat.timeNodeSel.CumulativeTime("sec") << "\n";
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

    double cputOther = stat.totalTime.CumulativeTime("sec") - stat.timeLBComputation.CumulativeTime("sec") - stat.timeUpdateUB.CumulativeTime("sec") - stat.timeDominanceTest.CumulativeTime("sec");
    std::cout << "\n Elapsed time:   " << stat.totalTime.CumulativeTime("sec") << " sec\n\n";
    std::cout << "    -> LB computation: " << stat.timeLBComputation.CumulativeTime("sec") << " sec ( " << 100 * stat.timeLBComputation.CumulativeTime("sec") / stat.totalTime.CumulativeTime("sec") << " % )\n";
    if (P.LBset == LP_RELAX || P.LBset == WARMSTARTED_LP_RELAX) {
        std::cout << "          -> Initialization: " << stat.timeInitialization.CumulativeTime("sec") << " sec ( " << 100 * stat.timeInitialization.CumulativeTime("sec") / stat.totalTime.CumulativeTime("sec") << " % )\n";
        std::cout << "          -> Feasibility checks: " << stat.timeFeasibilityCheck.CumulativeTime("sec") << " sec ( " << 100 * stat.timeFeasibilityCheck.CumulativeTime("sec") / stat.totalTime.CumulativeTime("sec") << " % )\n";
        //std::cout << "          -> Find furthest feasible point: " << stat.timeFurthestFeasiblePoint << " sec ( " << 100 * stat.timeFurthestFeasiblePoint / stat.totalTime << " % )\n";
        //std::cout << "          -> Find hyperplane (dual benson): " << stat.timeDualBenson << " sec ( " << 100 * stat.timeDualBenson / stat.totalTime << " % )\n";
        std::cout << "          -> Update polyhedron: " << stat.timeUpdatePolyhedron.CumulativeTime("sec") << " sec ( " << 100 * stat.timeUpdatePolyhedron.CumulativeTime("sec") / stat.totalTime.CumulativeTime("sec") << " % )\n";
    }
    std::cout << "    -> Dominance test: " << stat.timeDominanceTest.CumulativeTime("sec") << " sec ( " << 100 * stat.timeDominanceTest.CumulativeTime("sec") / stat.totalTime.CumulativeTime("sec") << " % )\n";
    std::cout << "    -> Update UB: " << stat.timeUpdateUB.CumulativeTime("sec") << " sec ( " << 100 * stat.timeUpdateUB.CumulativeTime("sec") / stat.totalTime.CumulativeTime("sec") << " % )\n";
    std::cout << "    -> others: " << cputOther << " sec ( " << 100 * cputOther / stat.totalTime.CumulativeTime("sec") << " % )\n";

    int nbFathomed = stat.nbFathomedInfeasibility + stat.nbFathomedDominance + stat.nbFathomedOptimality;
    std::cout << "\n Node explored: " << stat.nbNodes << "\n";

    if (P.LBset == LP_RELAX || P.LBset == WARMSTARTED_LP_RELAX) {
        std::cout << "\n LP solved: " << stat.lpSolved << "\n";
        //std::cout << "    -> Feasibility check: " << stat.feasibilityCheckSolved << "\n";
        //std::cout << "    -> Furthest feasible point: " << stat.furthestFeasbilePointSolved << "\n";
        //std::cout << "    -> Dual benson: " << stat.dualBensonSolved << "\n";
    }

    //std::cout << "\n PROFILER: " << stat.profiler << " sec ( " << 100 * stat.profiler / stat.totalTime << " % )\n";
}