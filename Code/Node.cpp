#include "Node.h"

/* ==========================================================
        Constructors
 ========================================================= */

/*! \brief Default constructor of a node.
 */
Node::Node() : P(nullptr), prob(NULL), varFix(NULL), ws(NULL), param(nullptr), LB(), UB(nullptr), score(-1), branchingDecision(), splittingIndex(-1), ndLub(0), lubDomi(0), stat(nullptr), depth(0), timeLBComputation(), timeDominanceTest(), timeUpdateUB(), coverCuts(0), iteration(-1), sLbFatherNode(0), sLbCurrentNode(0), sXiFatherNode(0), sXiCurrentNode(P->get_n(),-1) {} //, status(UNSOLVED)

/* \brief Destructor of the node.
 *
 */
Node::~Node() {
    if (param->LBset == LP_RELAX || param->LBset == WARMSTARTED_LP_RELAX || param->LBset == PRECOMPUTED_WARMSTARTED_LP_RELAX) {
        //dynamic_cast<LinearRelaxation*>(LB)->~LinearRelaxation();
        //delete LB;
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
Node::Node(MathematicalModel* lp, Parameters* par, Statistics* stat, UpperBoundSet* U) : P(lp), param(par), UB(U), score(-1), branchingDecision(lp->get_n(),lp->get_p()), splittingIndex(-1), ndLub(0), lubDomi(0), stat(stat), depth(0), timeLBComputation(), timeDominanceTest(), timeUpdateUB(), coverCuts(0), iteration(-1), sLbFatherNode(0), sLbCurrentNode(0), sXiFatherNode(0), sXiCurrentNode(P->get_n(), -1) { //, status(UNSOLVED)

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
        branchingDecision.slub[k] = 0;
        for (int i = 0; i < lp->get_n(); i++) {
            if (lp->get_objective(k, i) >= 0) {
                branchingDecision.slub[k] += lp->get_objective(k, i) * lp->getUb(i);
            }
        }
        branchingDecision.slub[k] += 2 + lp->get_n();
        std::cout << branchingDecision.slub[k] << "\n";
    }

    // define appropriate lb -> switch lb case ??
    if (par->LBset == LP_RELAX || par->LBset == WARMSTARTED_LP_RELAX || par->LBset == PRECOMPUTED_WARMSTARTED_LP_RELAX) {
        LB = new LinearRelaxation(lp,&branchingDecision,stat,par);
    }
    else {
        //std::cout << "  !!! Undefined lower bound set !!!\n";
        throw std::string("Error: Undefined lower bound set");
    }

    // probing model
    if (param->variableSelection == PROBING || param->variableSelection == VARIABLE_FIXING) {
        prob = new ProbingModel();
        prob->build(*lp, param);
        varFix = NULL;
        if (param->versionProbing == SCORE_WS)
            prob->setUpScoreWS(lp);
    }
    else if (param->variableSelection == PROBING_PRESOLVE_IP || param->variableSelection == PROBING_PRESOLVE_LP) {
        prob = new ProbingModel();
        prob->build(*lp, param);
        varFix = new VariableFixingModel();
        varFix->build(*lp, param);
        if (param->versionProbing == SCORE_WS)
            prob->setUpScoreWS(lp);
    }
    else {
        prob = NULL;
        varFix = NULL;
    }
    U->setCplexModels(dynamic_cast<LinearRelaxation*>(LB)->getFeasibilityCheckModel(), prob);

    if (param->nodeSelection == BEST_BOUND_WS || param->nodeSelection == BEST_BOUND_MAXMIN_GAP) {
        ws = new WeightedSumModel();
        ws->build(*lp);
    }
}

/*! \brief Creates a node given a splitting index and a bound.
 *
 * \param nd Node*. A pointer to the node it is created from.
 * \param index int. The index of the variable to split.
 * \param bound int. States the sign of the new constraint: are we adding an upper or a lower bound value on a variable?
 * \param val int. The value of the bound added on the variable, i.e. the right-hand side of the new branching constraint.
 */
Node::Node(Node& nd, int index, int bound, int val, SLUB& slub) : P(nd.P), prob(nd.prob), varFix(nd.varFix), ws(nd.ws), param(nd.param), UB(nd.UB), score(-1), splittingIndex(index), lubDomi(0), stat(nd.stat), depth(nd.depth + 1), timeLBComputation(), timeDominanceTest(), timeUpdateUB(), coverCuts(0), iteration(-1), sLbFatherNode(nd.sLbCurrentNode), sLbCurrentNode(0), sXiFatherNode(nd.sXiCurrentNode[index]), sXiCurrentNode(P->get_n(), -1) {


    stat->debug.StartTimer();
    ndLub = std::list<int>(nd.ndLub);
    stat->debug.StopTimer();
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
        LB = new LinearRelaxation(nd.P,cst->getWeightedSumModel(),cst->getFeasibilityCheckModel(),cst->getDualBensonModel(),cst->getFurthestFeasiblePointModel(),&branchingDecision,stat,nd.param);
    }
    else if (nd.param->LBset == WARMSTARTED_LP_RELAX || param->LBset == PRECOMPUTED_WARMSTARTED_LP_RELAX) {
        if (depth % CORRECTION_WARMSTART == 0) {
            LinearRelaxation* cst = dynamic_cast<LinearRelaxation*>(nd.LB);
            LB = new LinearRelaxation(nd.P, cst->getWeightedSumModel(), cst->getFeasibilityCheckModel(), cst->getDualBensonModel(), cst->getFurthestFeasiblePointModel(), &branchingDecision, stat,nd.param);
        }
        else {
            stat->timeCopyLB.StartTimer();
            LB = new LinearRelaxation(dynamic_cast<LinearRelaxation*>(nd.LB),&branchingDecision,nd.param);
            stat->timeCopyLB.StopTimer();
        }
    }
    else {
        throw std::string("Error: Lower bound set not supported for node splitting\n");
    }

    if (param->nodeSelection == BEST_BOUND_WS) {
        computeWsScore();
    }
    else if (param->nodeSelection == BEST_BOUND_MAXMIN_GAP) {
        computeMaxMinGapScore();
    }
    else if (param->nodeSelection == MOST_FRACTIONAL) {
        score = dynamic_cast<LinearRelaxation*>(LB)->getPercentageIntegrality();
    }
}


/*! \brief Creates a node given a splitting index and a bound.
 *
 * \param nd Node*. A pointer to the node it is created from.
 * \param index int. The index of the variable to split.
 * \param bound int. States the sign of the new constraint: are we adding an upper or a lower bound value on a variable?
 * \param val int. The value of the bound added on the variable, i.e. the right-hand side of the new branching constraint.
 */
Node::Node(Node& nd, BranchingDecisions* newbd, int index, SLUB& slub) : P(nd.P), prob(nd.prob), varFix(nd.varFix), ws(nd.ws), param(nd.param), UB(nd.UB), score(-1), splittingIndex(index), ndLub(nd.ndLub), stat(nd.stat), depth(nd.depth + 1), timeLBComputation(), timeDominanceTest(), timeUpdateUB(), coverCuts(0), iteration(-1), sLbFatherNode(nd.sLbCurrentNode), sLbCurrentNode(0), sXiFatherNode(nd.sXiCurrentNode[index]), sXiCurrentNode(P->get_n(), -1) {

    branchingDecision = BranchingDecisions(*newbd);
    branchingDecision.depth = depth;
    branchingDecision.lastSplittedIndex = index;

    if (nd.param->LBset == LP_RELAX) {
        LinearRelaxation* cst = dynamic_cast<LinearRelaxation*>(nd.LB);
        LB = new LinearRelaxation(nd.P, cst->getWeightedSumModel(), cst->getFeasibilityCheckModel(), cst->getDualBensonModel(), cst->getFurthestFeasiblePointModel(), &branchingDecision, stat, nd.param);
    }
    else if (nd.param->LBset == WARMSTARTED_LP_RELAX || param->LBset == PRECOMPUTED_WARMSTARTED_LP_RELAX) {
        if (depth % CORRECTION_WARMSTART == 0) {
            LinearRelaxation* cst = dynamic_cast<LinearRelaxation*>(nd.LB);
            LB = new LinearRelaxation(nd.P, cst->getWeightedSumModel(), cst->getFeasibilityCheckModel(), cst->getDualBensonModel(), cst->getFurthestFeasiblePointModel(), &branchingDecision, stat, nd.param);
        }
        else
            LB = new LinearRelaxation(dynamic_cast<LinearRelaxation*>(nd.LB), &branchingDecision, nd.param);
    }
    else {
        throw std::string("Error: Lower bound set not supported for node splitting\n");
    }

    if (ENABLE_EARLY_LB_DOMI) lubDomi = std::list<LocalUpperBound*>(nd.lubDomi);

    // compute LB and update UB
    if (param->LBset == PRECOMPUTED_WARMSTARTED_LP_RELAX) {
        LB->applyBranchingDecisions();
        if (param->generateCuts) generateCuts(branchingDecision);
        stat->timeLBComputation.StartTimer();
        LB->compute(UB, lubDomi);
        stat->timeLBComputation.StopTimer();

        stat->timeUpdateUB.StartTimer();
        LB->gatherIntegerSolutions(*UB);
        stat->timeUpdateUB.StopTimer();

        //dynamic_cast<LinearRelaxation*>(LB)->checkViolatingCut();
    }

    // compute score
    if (param->nodeSelection == BEST_BOUND_WS) {
        computeWsScore();
    }
    else if (param->nodeSelection == BEST_BOUND_MAXMIN_GAP) {
        computeMaxMinGapScore();
    }
}

/* ==========================================================
        Regular Methods
 ========================================================= */

 /*! \brief Process the node
  */
void Node::process(UpperBoundSet& U, int iteration) {

    //P->printObjective(branchingDecision);

    // update cplex models to fit the current node
    
    LB->setIteration(iteration);
    LB->applyBranchingDecisions();

    //double oldScore = score;
    //computeMaxMinGapScore();
    //std::cout << " score difference : " << oldScore << " vs " << score << "\n";

    // compute the LB set

    stat->timeLBComputation.StartTimer();
    if (iteration == 0 || param->LBset != PRECOMPUTED_WARMSTARTED_LP_RELAX) LB->compute(UB, lubDomi);
    stat->timeLBComputation.StopTimer();

    //dynamic_cast<LinearRelaxation*>(LB)->checkViolatingCut();

    // gather statistics

    if (LB->getStatus() == UNSOLVED) {
        throw std::string("Error: lower bound set unsolved");
    }
    else if (LB->getStatus() != INFEASIBLE) {

        // get LB statistics
        if (param->LBset == LP_RELAX || param->LBset == WARMSTARTED_LP_RELAX) {
            int f = dynamic_cast<LinearRelaxation*>(LB)->get_numberFacets();
            int v = dynamic_cast<LinearRelaxation*>(LB)->get_numberVertices();
            if (stat->avgFacets.size() == depth) { // a new depth is reached
                stat->avgFacets.push_back(f);
                stat->minFacets.push_back(f);
                stat->maxFacets.push_back(f);
                stat->avgVertices.push_back(v);
                stat->minVertices.push_back(v);
                stat->maxVertices.push_back(v);
                stat->nbNodesAtStage.push_back(1);
                stat->timeUpdatePoly.push_back(stat->cpuCurrent);
                stat->nbFeasVtx.push_back(stat->nbFeasVtxCurrent);
                stat->nbNewFacets.push_back(stat->nbNewFacetsCurrent);
                stat->vtxOnly.push_back(stat->vtxOnlyCurrent);
                stat->facetsNoRay.push_back(stat->facetsNoRayCurrent);
                stat->newFacetsNoRay.push_back(stat->newFacetsNoRayCurrent);
            }
            else {
                stat->avgFacets[depth] += f;
                if (stat->minFacets[depth] > f) stat->minFacets[depth] = f;
                if (stat->maxFacets[depth] < f) stat->maxFacets[depth] = f;
                stat->avgVertices[depth] += v;
                if (stat->minVertices[depth] > f) stat->minVertices[depth] = f;
                if (stat->maxVertices[depth] < f) stat->maxVertices[depth] = f;
                stat->nbNodesAtStage[depth]++;
                stat->timeUpdatePoly[depth] += stat->cpuCurrent;
                stat->nbFeasVtx[depth] += stat->nbFeasVtxCurrent;
                stat->nbNewFacets[depth] += stat->nbNewFacetsCurrent;
                stat->vtxOnly[depth] += stat->vtxOnlyCurrent;
                stat->facetsNoRay[depth] += stat->facetsNoRayCurrent;
                stat->newFacetsNoRay[depth] += stat->newFacetsNoRayCurrent;
            }
            stat->cpuCurrent = 0;
            stat->nbFeasVtxCurrent = 0;
            stat->nbNewFacetsCurrent = 0;
            stat->vtxOnlyCurrent = 0;
            stat->facetsNoRayCurrent = 0;
            stat->newFacetsNoRayCurrent = 0;
        }

        // Update UB
        stat->timeUpdateUB.StartTimer();
        if (iteration == 0 || param->LBset != PRECOMPUTED_WARMSTARTED_LP_RELAX) LB->gatherIntegerSolutions(U);
        stat->timeUpdateUB.StopTimer();

        if (stat->avgIntegerVtx.size() == depth) { // a new depth is reached
            stat->avgIntegerVtx.push_back(stat->avgIntegerVtxCurrent);
        }
        else {
            stat->avgIntegerVtx[depth] += stat->avgIntegerVtxCurrent;
        }
        stat->avgIntegerVtxCurrent = 0;

        if (LB->getStatus() != OPTIMAL) {
            // domi test
            stat->timeDominanceTest.StartTimer();
            LB->applyDominanceTest(U,param,ndLub,&lubDomi);
            stat->timeDominanceTest.StopTimer();
        }
    }
    else {
        // LB stat
        if (stat->nbInfeasibleNodes.size() <= depth) {
            while (stat->nbInfeasibleNodes.size() < depth) {
                stat->nbInfeasibleNodes.push_back(0);
            }
            stat->nbInfeasibleNodes.push_back(1);
        }
        else {
            stat->nbInfeasibleNodes[depth]++;
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
//void Node::splitOS(std::list<Node*>* Q, UpperBoundSet* U, int iteration) {
void Node::splitOS(TreeManager* T, UpperBoundSet* U, int iteration) {

    if (param->objectiveBranching == NO_OBJECTIVE_BRANCHING) {
        stat->timeComputeOB.StartTimer();
        SLUB slub = SLUB(branchingDecision, param);
        stat->timeComputeOB.StopTimer();

        int nbLub = U->getLubs()->size();
        stat->avgNbLubs += nbLub;
        if (nbLub >= stat->maxNbLubs) stat->maxNbLubs = nbLub;

        splitVS(T, slub, iteration, U);
    }
    else if (param->objectiveBranching == CONE_OBJECTIVE_BRANCHING) {

        //if (iteration == 611) std::cout << "Hi there!\n";

        stat->timeComputeOB.StartTimer();
        SLUB slub = SLUB(P->get_p());
        //slub.print();
        std::list<LocalUpperBound>* NU = U->getLubs();
        std::list<LocalUpperBound>::iterator u;
        std::list<LocalUpperBound*>::iterator u2;
        std::list<int>::iterator nextId = ndLub.begin();
        bool dominated;
        // search for dominated lubs to compute OB. Note that all non-existing lubs have been deleted beforehand
        // in the dominance test, and thus there is no need to take care of this case.
        //for (u = NU->begin(); u != NU->end(); u++) {
        //    dominated = true;
        //    if (nextId != ndLub.end() && *nextId == u->get_id()) {
        //        dominated = false;
        //        nextId++;
        //    }
        //    if (dominated) {
        //        slub.merge(*u);
        //    }
        //    //slub.print();
        //}
        for (u2 = lubDomi.begin(); u2 != lubDomi.end(); u2++) {
            slub.merge(**u2);
            //if (iteration == 1355)
                //(*u2)->print();
        }
        stat->timeComputeOB.StopTimer();

        int nbLub = U->getLubs()->size();
        stat->avgNbLubs += nbLub;
        if (nbLub >= stat->maxNbLubs) stat->maxNbLubs = nbLub;

        /*if (iteration == 611) {
            std::cout << "\n  -> final slub : ";
            slub.print();
        }*/

        splitVS(T, slub, iteration, U);
    }
    else if (param->objectiveBranching == FULL_OBJECTIVE_BRANCHING) {

        stat->timeComputeOB.StartTimer();
        // building the initial set of SLUBs S.
        std::list<SLUB*> S(0);
        std::list<LocalUpperBound>* lubs = U->getLubs();
        std::list<LocalUpperBound>::iterator u;
        std::list<LocalUpperBound*>::iterator u2;
        std::list<int>::iterator nextId = ndLub.begin();
        //for (u = lubs->begin(); u != lubs->end(); u++) { // take only non-dominated ones !!!
        //    if (nextId != ndLub.end() && *nextId == u->get_id()) {
        //        nextId++;
        //    }
        //    else {
        //        S.push_back(new SLUB(*u));
        //    }
        //}
        //stat->cpuSlubsOB.StartTimer();
        for (u2 = lubDomi.begin(); u2 != lubDomi.end(); u2++) {
            //if (iteration == 1170) {
                //std::cout << " new lub considered : ";
                //(*u2)->print();
            //}
            S.push_back(new SLUB(**u2));
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
        //stat->cpuSlubsOB.StopTimer();

        stat->timeComputeOB.StopTimer();
        //if (iteration == 1170) std::cout << " sub-pb created: " << S.size() << std::endl;
        for (s1 = S.begin(); s1 != S.end(); s1++) {
            //if (iteration == 1170) (*s1)->print();
            splitVS(T, **s1, iteration, U);
        }

        //if (S.size() >= 2) std::cout << S.size() << " problems created in the objective space\n";

        int nbLub = U->getLubs()->size();
        stat->avgNbLubs += nbLub;
        if (nbLub >= stat->maxNbLubs) stat->maxNbLubs = nbLub;
        if (S.size() >= 2) {
            stat->avgDepthOB += depth;
            stat->avgSubPbOB += S.size();
            stat->nbOccurencesOB++;
            if (stat->minDepthOB >= depth) stat->minDepthOB = depth;
            if (stat->maxDepthOB <= depth) stat->maxDepthOB = depth;
            if (stat->minSubPbOB >= S.size()) stat->minSubPbOB = S.size();
            if (stat->maxSubPbOB <= S.size()) stat->maxSubPbOB = S.size();

            // stat per depth
            if (stat->nodesOBPerDepth.size() <= depth) {
                while (stat->nodesOBPerDepth.size() < depth) {
                    stat->nodesOBPerDepth.push_back(0);
                    stat->avgSubPbOBPerDepth.push_back(0);
                    stat->minSubPbOBPerDepth.push_back(100000);
                    stat->maxSubPbOBPerDepth.push_back(0);
                }
                stat->nodesOBPerDepth.push_back(1);
                stat->avgSubPbOBPerDepth.push_back(S.size());
                stat->minSubPbOBPerDepth.push_back(S.size());
                stat->maxSubPbOBPerDepth.push_back(S.size());
            }
            else {
                stat->nodesOBPerDepth[depth]++;
                stat->avgSubPbOBPerDepth[depth] += S.size();
                if (stat->minSubPbOBPerDepth[depth] >= S.size()) stat->minSubPbOBPerDepth[depth] = S.size();
                if (stat->maxSubPbOBPerDepth[depth] <= S.size()) stat->maxSubPbOBPerDepth[depth] = S.size();
            }
        }

    }
    else if (param->objectiveBranching == LIMITED_OBJECTIVE_BRANCHING) {

        stat->timeComputeOB.StartTimer();
        // building the initial set of SLUBs S.
        std::list<SLUB*> S(0);
        std::list<LocalUpperBound>* lubs = U->getLubs();
        std::list<LocalUpperBound>::iterator u;
        std::list<LocalUpperBound*>::iterator u2;
        std::list<int>::iterator nextId = ndLub.begin();
        std::list<SLUB*>::iterator s1, s2, sMerge1, sMerge2;

       for (u2 = lubDomi.begin(); u2 != lubDomi.end(); u2++) {
            S.push_back(new SLUB(**u2));
        }

        // merging the SLUBs if necessary
       mergeSlubs(&S);

        // [TODO] add additional merging here
        double dist, smallestDist;
        while (S.size() > param->limitSubPb) { 
            
            smallestDist = DBL_MAX;
            // merge the two closest slubs
            s1 = S.begin();
            while (s1 != S.end()) {
                s2 = s1;
                s2++;
                while (s2 != S.end()) {
                    //std::cout << "dist is : " << (*s1)->distance(*s2) << "\n";
                    if ((*s1)->distance(*s2) <= smallestDist) {
                        sMerge1 = s1;
                        sMerge2 = s2;
                        smallestDist = (*s1)->distance(*s2);
                    }
                    s2++;
                }
                s1++;
            }
            (*sMerge1)->merge(**sMerge2);
            delete* sMerge2;
            S.erase(sMerge2);
            //std::cout << "MERGED!!\n\n";

            // merge subproblems that became redundant
            mergeSlubs(&S);
        }

        stat->timeComputeOB.StopTimer();
        //std::cout << " sub-pb created: " << S.size() << std::endl;
        for (s1 = S.begin(); s1 != S.end(); s1++) {
            //(*s1)->print();
            splitVS(T, **s1, iteration, U);
        }


        //if (S.size() >= 2) std::cout << S.size() << " problems created in the objective space\n";

        int nbLub = U->getLubs()->size();
        stat->avgNbLubs += nbLub;
        if (nbLub >= stat->maxNbLubs) stat->maxNbLubs = nbLub;
        if (S.size() >= 2) {
            stat->avgDepthOB += depth;
            stat->avgSubPbOB += S.size();
            stat->nbOccurencesOB++;
            if (stat->minDepthOB >= depth) stat->minDepthOB = depth;
            if (stat->maxDepthOB <= depth) stat->maxDepthOB = depth;
            if (stat->minSubPbOB >= S.size()) stat->minSubPbOB = S.size();
            if (stat->maxSubPbOB <= S.size()) stat->maxSubPbOB = S.size();

            // stat per depth
            if (stat->nodesOBPerDepth.size() <= depth) {
                while (stat->nodesOBPerDepth.size() < depth) {
                    stat->nodesOBPerDepth.push_back(0);
                    stat->avgSubPbOBPerDepth.push_back(0);
                    stat->minSubPbOBPerDepth.push_back(100000);
                    stat->maxSubPbOBPerDepth.push_back(0);
                }
                stat->nodesOBPerDepth.push_back(1);
                stat->avgSubPbOBPerDepth.push_back(S.size());
                stat->minSubPbOBPerDepth.push_back(S.size());
                stat->maxSubPbOBPerDepth.push_back(S.size());
            }
            else {
                stat->nodesOBPerDepth[depth]++;
                stat->avgSubPbOBPerDepth[depth] += S.size();
                if (stat->minSubPbOBPerDepth[depth] >= S.size()) stat->minSubPbOBPerDepth[depth] = S.size();
                if (stat->maxSubPbOBPerDepth[depth] <= S.size()) stat->maxSubPbOBPerDepth[depth] = S.size();
            }
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
//void Node::splitVS(std::list<Node*>* Q, SLUB& slub, int iteration, UpperBoundSet* U) {
void Node::splitVS(TreeManager* T, SLUB & slub, int iteration, UpperBoundSet * U) {

    //if (iteration > -1) {
    //    std::cout << "\nOB: ";
    //    slub.print();
    //    //std::cout << "\n";
    //}

    stat->timeVarSel.StartTimer();
    //std::vector<int> index(3,0); // 0: index splitted, 1: -ub 1st node, as a negative number, 2: lb 2nd node
    if (param->variableSelection == FIRST_INDEX) {
        if (P->isBinary()) {
            //Q->push_back(new Node(*this, splittingIndex + 1, IS_UB, 0, slub));
            //Q->push_back(new Node(*this, splittingIndex + 1, IS_LB, 1, slub));
            T->pushNode(new Node(*this, splittingIndex + 1, IS_UB, 0, slub));
            T->pushNode(new Node(*this, splittingIndex + 1, IS_LB, 1, slub));
        }
        else {
            throw std::string("Error: integer programs not supported by the current variable selection.");
        }
    }
    else if (param->variableSelection == MOST_OFTEN_FRACTIONAL) {
        if (param->LBset == LP_RELAX || param->LBset == WARMSTARTED_LP_RELAX || param->LBset == PRECOMPUTED_WARMSTARTED_LP_RELAX) {
            int index = dynamic_cast<LinearRelaxation*>(LB)->computeMostOftenFractionalIndex(slub);
            if (P->isBinary()) {
                //stat->debug.StartTimer();
                T->pushNode(new Node(*this, index, IS_UB, 0, slub));
                T->pushNode(new Node(*this, index, IS_LB, 1, slub));
                //stat->debug.StopTimer();
                //if (iteration == 611) std::cout << "two nodes created...\n";
            }
            else {
                int splittingValue = 0;
                if (param->branchingValueSelection == MEDIAN)
                    splittingValue = dynamic_cast<LinearRelaxation*>(LB)->computeMedianSplittingValue(slub, index, branchingDecision.ub[index]);
                else if (param->branchingValueSelection == MOST_OFTEN_FRACTIONAL_VALUE)
                    splittingValue = dynamic_cast<LinearRelaxation*>(LB)->computeMostOftenFractionalSplittingValue(slub, index, branchingDecision.ub[index]);
                else if (param->branchingValueSelection == RANDOM)
                    splittingValue = dynamic_cast<LinearRelaxation*>(LB)->computeRandomSplittingValue(slub, index, branchingDecision.ub[index]);
                else
                    throw std::string("Error: branching value selection not supported.");
                T->pushNode(new Node(*this, index, IS_UB, splittingValue, slub));
                T->pushNode(new Node(*this, index, IS_LB, splittingValue + 1, slub));
                //if (iteration == 1436)
                    //std::cout << " split value is " << splittingValue << " on variable " << splittingIndex << "\n";
                //std::cout << "\niteration " << iteration << " : x_" << index << " <= " << splittingValue << "\n\n";
            }
        }
        else {
            throw std::string("Error: variable selection parameter not supported with this LB set\n");
        }
    }
    else if (param->variableSelection == PROBING || param->variableSelection == PROBING_PRESOLVE_LP || param->variableSelection == PROBING_PRESOLVE_IP) {
        if (P->isBinary()) { // limit to binary pb for now
            BranchingDecisions* newbd0 = new BranchingDecisions(branchingDecision);
            BranchingDecisions* newbd1 = new BranchingDecisions(branchingDecision);
            for (int k = 0; k < P->get_p(); k++) {
                newbd0->slub[k] = slub.get_coordinate(k); // max(slub.get_coordinate(k), 1000000) + 1;
                newbd1->slub[k] = slub.get_coordinate(k); // max(slub.get_coordinate(k), 1000000) + 1;
            }
            stat->cpuProbing.StartTimer();
            //generateCuts(*newbd0);
            bool closeNode = fixBounds(newbd0, newbd1, U, iteration);
            stat->cpuProbing.StopTimer();

            if (!closeNode) {
                if (param->LBset == LP_RELAX || param->LBset == WARMSTARTED_LP_RELAX) {

                    /*int index = -1;
                    if (param->versionProbing == OLD_RULE)
                        index = dynamic_cast<LinearRelaxation*>(LB)->computeMostOftenFractionalIndex(slub, newbd0);
                    else if (param->versionProbing == CLOSEST_DIMENSION)
                        index = selectBranchingIndex1(newbd);
                    else if (param->versionProbing == WORST_BEST_YI)
                        index = selectBranchingIndex2(newbd);
                    else
                        throw std::string("Invalid probind method.\n");*/

                    // create children nodes
                    //newbd->lastSplittedIndex = index;
                    //newbd->lb[index] = 0;
                    //newbd->ub[index] = 0;
                    T->pushNode(new Node(*this, newbd0, newbd0->lastSplittedIndex, slub));
                    //newbd->ub[index] = 1;
                    //newbd->lb[index] = 1;
                    T->pushNode(new Node(*this, newbd1, newbd1->lastSplittedIndex, slub));
                    //if (iteration == 611) std::cout << "two nodes created.\n";
                }
                else {
                    throw std::string("Error: variable selection parameter not supported with this LB set\n");
                }
            }
            //else {
                //if (iteration == 611) std::cout << "we closed the node with variable fixing.\n";
            //}

            delete newbd0;
            delete newbd1;
        }
        else {
            throw std::string("Error: probing not supported with integer variables yet\n");
        }
    }
    else if (param->variableSelection == VARIABLE_FIXING) {
        if (P->isBinary()) {

            // set up the two branching decisions that will be used for branching

            BranchingDecisions* newbd0 = new BranchingDecisions(branchingDecision);
            BranchingDecisions* newbd1 = new BranchingDecisions(branchingDecision);
            for (int k = 0; k < P->get_p(); k++) {
                newbd0->slub[k] = slub.get_coordinate(k);
                newbd1->slub[k] = slub.get_coordinate(k);
            }
            //dynamic_cast<LinearRelaxation*>(LB)->getAverageNormalVector(newbd0);

            // fix variables and select a variable to branching on accordingly

            stat->cpuProbing.StartTimer();
            bool closeNode = fixAndSelectBranchingVariable1(newbd0, newbd1, U);
            stat->cpuProbing.StopTimer();

            // create the children nodes if the node is not fathomed with variable fixing

            if (!closeNode) {
                if (param->LBset == LP_RELAX || param->LBset == WARMSTARTED_LP_RELAX || param->LBset == PRECOMPUTED_WARMSTARTED_LP_RELAX) {
                    if (P->get_objDir(0) == 1) {
                        T->pushNode(new Node(*this, newbd1, newbd1->lastSplittedIndex, slub));
                        T->pushNode(new Node(*this, newbd0, newbd0->lastSplittedIndex, slub));
                    }
                    else {
                        T->pushNode(new Node(*this, newbd0, newbd0->lastSplittedIndex, slub));
                        T->pushNode(new Node(*this, newbd1, newbd1->lastSplittedIndex, slub));
                    }
                }
                else {
                    throw std::string("Error: variable selection parameter not supported with this LB set\n");
                }
            }

        }
        else {
            throw std::string("Error: variable fixing does not support integer variables\n");
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
        if (branchingDecision.lb[i] > P->getLb(i)) { //change for integer var + update lb & ub only if relevant when branching?
            std::cout << " x[" << i << "] >= " << branchingDecision.lb[i] << std::endl;
        }
        if (branchingDecision.ub[i] < P->getUb(i)) { //change for integer var + update lb & ub only if relevant when branching?
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
    std::cout << "lb :";
    for (int i = 0; i < P->get_n(); i++) {
        std::cout << " " << branchingDecision.lb[i];
    }
    std::cout << "\nub :";
    for (int i = 0; i < P->get_n(); i++) {
        std::cout << " " << branchingDecision.ub[i];
    }
    std::cout << "\n";
}

bool Node::isOurCulprit() {

    if (branchingDecision.ub[29] == 0 &&
        branchingDecision.ub[4] == 0 &&
        branchingDecision.ub[9] == 0 &&
        branchingDecision.ub[14] == 0 &&
        branchingDecision.ub[19] == 0 &&
        branchingDecision.ub[24] == 0 &&
        branchingDecision.lb[27] == 1 &&
        branchingDecision.lb[7] == 1 &&
        branchingDecision.ub[5] == 0 &&
        branchingDecision.ub[6] == 0 &&
        branchingDecision.ub[8] == 0 &&
        branchingDecision.ub[18] == 0 &&
        branchingDecision.lb[0] == 1 &&
        branchingDecision.ub[1] == 0 &&
        branchingDecision.ub[2] == 0 &&
        branchingDecision.ub[3] == 0 &&
        branchingDecision.lb[10] == 1 &&
        branchingDecision.ub[11] == 0 &&
        branchingDecision.ub[12] == 0 &&
        branchingDecision.ub[13] == 0 &&
        branchingDecision.ub[15] == 0 &&
        branchingDecision.lb[17] == 1 &&
        branchingDecision.ub[16] == 0 &&
        branchingDecision.ub[21] == 0 &&
        branchingDecision.ub[20] == 0 &&
        branchingDecision.lb[25] == 1 &&
        (branchingDecision.lb[22] == 1 ||
        branchingDecision.ub[23] == 0 ||
        branchingDecision.ub[26] == 0 ||
        branchingDecision.ub[28] == 0)
        
        )
        return true;
    else
        return false;

}

/*! \brief Try to fix integer bounds to each variable by solving LPs and looking at their feasibility, and chose the branching variable (index).
 *
 * The branching decision is stored in newbd, as well as the induces bounds fixed by this decision.
 * \param BranchingDecisions newBd. Describes a potential futur subproblem.
 * \return true if this node is detected as infeasible after checking the LPs.
 */

bool Node::fixBounds(BranchingDecisions* newbd0, BranchingDecisions* newbd1, UpperBoundSet* U, int iteration) {

    std::vector<BranchingDecisions> dec0(P->get_n());
    std::vector<BranchingDecisions> dec1(P->get_n());
    for (int i = 0; i < P->get_n(); i++) {
        dec0[i] = BranchingDecisions(*newbd0);
        dec1[i] = BranchingDecisions(*newbd1);
        // don't fix bounds here, otherwise we will not enter the if ub != lb
    }

    bool closeNode = false;
    bool feasible;

    for (int i = 0; i != P->get_n(); i++) {
        if (!closeNode && newbd0->lb[i] != newbd0->ub[i]) {

            dec1[i].lb[i] = 1;
            feasible = checkSubProblem(&dec1[i], i, 1);
            //std::cout << "feasible: " << feasible << "\n";
            //if (iteration == 1276) std::cout << " x[" << i << "] = 1 is feasible: " << feasible << "\n";

            if (!feasible) { // x_i = 1 is not feasible
                dec0[i].ub[i] = 0;
                feasible = checkSubProblem(&dec0[i], i, 0);
                //std::cout << "feasible: " << feasible << "\n";
                //if (iteration == 1276) std::cout << " x[" << i << "] = 0 is feasible: " << feasible << "\n";
                if (!feasible) {
                    closeNode = true;
                    //print();
                    //std::cout << "We closed the node (x[" << i << "] has no integral value)\n";
                }
                else {
                    newbd0->ub[i] = 0;
                    newbd1->ub[i] = 0;
                    dec0[i].score = 0;
                    dec1[i].score = 0;
                    //std::cout << "We fixed x[" << i << "] to 0\n";

                    // that is true for all branching decisions, so we fix this decision for all
                    for (int j = 0; j < P->get_n(); j++) {
                        dec0[j].ub[i] = 0;
                        dec1[j].ub[i] = 0;
                    }
                }
            }
            else {
                dec0[i].ub[i] = 0;
                //newbd->ub[i] = 0;
                feasible = checkSubProblem(&dec0[i], i, 0);
                //if (iteration == 1276) std::cout << " x[" << i << "] = 0 is feasible: " << feasible << "\n";
                if (!feasible) {
                    //newbd->lb[i] = 1;
                    //newbd->ub[i] = 1;
                    //std::cout << "We fixed x[" << i << "] to 1\n";
                    newbd0->lb[i] = 1;
                    newbd1->lb[i] = 1;
                    dec0[i].score = 0;
                    dec1[i].score = 0;

                    // that is true for all branching decisions, so we fix this decision for all
                    for (int j = 0; j < P->get_n(); j++) {
                        dec0[j].lb[i] = 1;
                        dec1[j].lb[i] = 1;
                    }
                }
                else {
                    //newbd->ub[i] = 1;
                }
            }

        }
    }

    // if we are not fathomed by infeasbility, we test optimality
    if (!closeNode) {

        // we check if all bounds are fixed
        closeNode = true;
        for (int i = 0; i != P->get_n(); i++) {
            if (newbd0->lb[i] != newbd0->ub[i]) { // 0 and 1 are the same ?
                closeNode = false;
                break;
            }
        }

        // if they are, we add the point to the ub set
        if (closeNode) {
            //print();
            //std::cout << "We close the node (all variables has been fixed)\n";
            U->updateUB(newbd0);
        }
    }


    // we select the next branching variable

    if (!closeNode) {
        int index = -1;
        if (param->versionProbing == OLD_RULE) {
            SLUB s = SLUB(newbd0->slub);
            index = dynamic_cast<LinearRelaxation*>(LB)->computeMostOftenFractionalIndex(s, newbd0);
        }
        else if (param->versionProbing == CLOSEST_DIMENSION)
            index = selectBranchingIndex1(newbd0, newbd1, dec0, dec1);
        else if (param->versionProbing == SCORE_WS)
            index = selectBranchingMinScore(newbd0, newbd1, dec0, dec1);
        else
            throw std::string("Invalid probind method.\n");


        // we apply the decision to the branchingDecisions

        for (int i = 0; i < P->get_n(); i++) {
            // lbs
            if (newbd0->lb[i] == 0 && dec0[index].lb[i] == 1)
                newbd0->lb[i] = 1;
            else if (newbd1->lb[i] == 0 && dec1[index].lb[i] == 1)
                newbd1->lb[i] = 1;

            // ubs
            if (newbd0->ub[i] == 1 && dec0[index].ub[i] == 0)
                newbd0->ub[i] = 0;
            else if (newbd1->ub[i] == 1 && dec1[index].ub[i] == 0)
                newbd1->ub[i] = 1;
        }
        newbd0->score = dec0[index].score;
        newbd1->score = dec1[index].score;

        newbd0->lastSplittedIndex = index;
        newbd1->lastSplittedIndex = index;
    }

    return closeNode;
}


bool Node::checkSubProblem(BranchingDecisions* newBd, int i, int v) {

    bool feasible = true;

    // presolve and fix variables
    if (param->variableSelection == PROBING_PRESOLVE_IP || param->variableSelection == PROBING_PRESOLVE_LP) {
        varFix->resetBounds(*newBd);
        feasible = varFix->presolveAndFixVariables(newBd);
    }

    // if relevant, compute ideal point or score
    if (feasible) {
        if (param->versionProbing == CLOSEST_DIMENSION) {
            for (int k = 0; k < P->get_p(); k++) {
                stat->cpuProbingCplex.StartTimer();
                feasible = prob->solve(*newBd, k, P);
                stat->cpuProbingCplex.StopTimer();
                stat->lpProbing++;
                //std::cout << "feas: " << feasible << "\n"; 
                if (feasible) { // why condition 2 ?    && newBd->yIprobing[k] >= prob->retrieveObjectiveValue()
                    newBd->yIprobing[k] = prob->retrieveObjectiveValue();
                    //std::cout << newBd->yIprobing[k] << " \n";
                }
                else
                    break;
            }
        }
        else if (param->versionProbing == SCORE_WS || param->versionProbing == OLD_RULE) {
            stat->cpuProbingCplex.StartTimer();
            feasible = prob->solve(*newBd);
            stat->cpuProbingCplex.StopTimer();
            stat->lpProbing++;
            if (feasible) {
                int v = 0;
                for (int k = 0; k < P->get_p(); k++)
                    v += newBd->slub[k];
                newBd->score = v - prob->retrieveObjectiveValue();
            }
        }
    }

    // reset fixed variables
    //for (std::list<int>::iterator idx = newBd->resetVar.begin(); idx != newBd->resetVar.end(); idx++) {
    //    //std::cout << "we reset.\n";
    //    newBd->lb[*idx] = 0;
    //    newBd->ub[*idx] = 1;
    //}
    //newBd->resetVar.clear();

    return feasible;
}

/* \brief Generate cuts to the node.
 *
 * Generate cover inequalities for minimization problems for OB constraints. Add them to the probing model and to the linear relaxation.
 */
void Node::generateCuts() {

    // clear old cuts from cplex models

    prob->clearCuts(); // if (param->variableSelection == PROBING) 
    dynamic_cast<LinearRelaxation*>(LB)->clearCuts(); // if (param->LBset == LP_RELAX || param->LBset == WARMSTARTED_LP_RELAX) 

    // compute the new cover inequalities

    computeCoverCuts(); // largest to smallest coef

    // add them to the cplex models

    //if (iteration != 316) {
        prob->applyCoverCuts(coverCuts); // if (param->variableSelection == PROBING) 
        dynamic_cast<LinearRelaxation*>(LB)->applyCoverCuts(coverCuts); // if (param->LBset == LP_RELAX || param->LBset == WARMSTARTED_LP_RELAX) 
    //}
}

/* \brief Compute the cover inequalities -> largest to smallest obj coefficient
*/
void Node::computeCoverCuts() {

    int rhs = 0;

    for (int k = 0; k < P->get_p(); k++) { // compute cuts for each objective...

        coverCuts.push_back(std::vector<int>(0));
        if (P->get_objDir(k) == -1 && branchingDecision.slub[k] < INT_MAX - 2 && iteration != 0) { // ... if this objective function is minimized and bounded by an objective branching constraint
            
            // compute the actual rhs by substracting variables fixed to 1 with a positive objective coefficient
            rhs = branchingDecision.slub[k];
            for (int i = 0; i < P->get_n(); i++) {
                if (branchingDecision.lb[i] == 1) {
                    rhs -= max(P->get_objective(k, i), 0);
                }
            }
            coverCuts[k].push_back(rhs);

            // cut #1: largest coefficients + breaking items
            int idx, maxSizeClique = -1, i = 0;
            while (rhs >= 0 && i < P->get_n()) {
                idx = P->get_sortedIndex(k, i);
                if (branchingDecision.lb[idx] != branchingDecision.ub[idx]) {
                    coverCuts[k].push_back(idx);
                    rhs -= P->get_objective(k, idx);
                    maxSizeClique++;
                    //std::cout << "rhs = " << rhs << "\n";
                }
                i++;
            }
            if (i == P->get_n()) {
                maxSizeClique++;
            }
            else { // we try to improve the inequality by adding new variables to the cut
                int j = 0;
                idx = P->get_sortedIndex(k, i);
                int idx2 = P->get_sortedIndex(k, j);
                while (rhs + P->get_objective(k, idx2) - P->get_objective(k, idx) <= 0 && i < P->get_n()) {
                    if (branchingDecision.lb[idx] != branchingDecision.ub[idx]) {
                        coverCuts[k].push_back(idx);
                        rhs = rhs + P->get_objective(k, idx2) - P->get_objective(k, idx);
                        j++;
                        //std::cout << "rhs renewed = " << rhs << "\n";
                    }
                    i++;
                    if (i != P->get_n()) {
                        idx = P->get_sortedIndex(k, i);
                        idx2 = P->get_sortedIndex(k, j);
                    }
                }
            }
            coverCuts[k][0] = maxSizeClique;
            
            // print
            /*if (coverCuts[k].size() >= 2) {
                std::cout << "\n -> a new cut was generated on objective " << k + 1 << " : ";
                for (int l = 1; l < coverCuts[k].size() - 1; l++) {
                    std::cout << "x[" << coverCuts[k][l] << "] + ";
                }
                std::cout << "x[" << coverCuts[k][coverCuts[k].size() - 1] << "] <= " << maxSizeClique << "\n";
            }
            else {
                std::cout << "\n --> No cut generated !\n";
            }*/
        }

    }

}

/* \brief Generate cuts to the node.
 *
 * Generate cover inequalities for minimization problems for OB constraints. Add them to the probing model and to the linear relaxation.
 */
void Node::generateCuts(BranchingDecisions& bd) {

    /*std::cout << " Generate new cover cuts on : ";
    for (int k = 0; k < P->get_p(); k++) {
        std::cout << bd.slub[k] << " ";
    }
    std::cout << "\n";*/

    // clear old cuts from cplex models

    prob->clearCuts(); // if (param->variableSelection == PROBING) 
    dynamic_cast<LinearRelaxation*>(LB)->clearCuts(); // if (param->LBset == LP_RELAX || param->LBset == WARMSTARTED_LP_RELAX) 

    // compute the new cover inequalities

    //std::cout << "hai\n";
    computeCoverCuts(bd);
    //computeMofCoverCuts(bd);

    // add them to the cplex models

    //if (iteration != 316) {
    prob->applyCoverCuts(coverCuts); // if (param->variableSelection == PROBING) 
    dynamic_cast<LinearRelaxation*>(LB)->applyCoverCuts(coverCuts); // if (param->LBset == LP_RELAX || param->LBset == WARMSTARTED_LP_RELAX) 
    //}
}

/* \brief Compute the cover inequalities
*/
void Node::computeCoverCuts(BranchingDecisions& bd) {

    int rhs = 0;

    for (int k = 0; k < P->get_p(); k++) { // compute cuts for each objective...

        if (coverCuts.size() <= k) coverCuts.push_back(std::vector<int>(0));
        else coverCuts[k] = std::vector<int>(0);

        if (P->get_objDir(k) == -1 && bd.slub[k] < INT_MAX - 2 && iteration != 0) { // ... if this objective function is minimized and bounded by an objective branching constraint

            // compute the actual rhs by substracting variables fixed to 1 with a positive objective coefficient
            rhs = bd.slub[k];
            for (int i = 0; i < P->get_n(); i++) {
                if (bd.lb[i] == 1) {
                    rhs -= max(P->get_objective(k, i), 0);
                }
            }
            coverCuts[k].push_back(rhs);

            // cut #1: largest coefficients + breaking items
            int idx, maxSizeClique = -1, i = 0;
            while (rhs >= 0 && i < P->get_n()) {
                idx = P->get_sortedIndex(k, i);
                if (bd.lb[idx] != bd.ub[idx]) {
                    coverCuts[k].push_back(idx);
                    rhs -= P->get_objective(k, idx);
                    maxSizeClique++;
                    //std::cout << "rhs = " << rhs << "\n";
                }
                i++;
            }
            if (i == P->get_n()) {
                maxSizeClique++;
            }
            else { // we try to improve the inequality by adding new variables to the cut
                int j = 0;
                idx = P->get_sortedIndex(k, i);
                int idx2 = P->get_sortedIndex(k, j);
                while (rhs + P->get_objective(k, idx2) - P->get_objective(k, idx) <= 0 && i < P->get_n()) {
                    if (bd.lb[idx] != bd.ub[idx]) {
                        coverCuts[k].push_back(idx);
                        rhs = rhs + P->get_objective(k, idx2) - P->get_objective(k, idx);
                        j++;
                        //std::cout << "rhs renewed = " << rhs << "\n";
                    }
                    i++;
                    if (i != P->get_n()) {
                        idx = P->get_sortedIndex(k, i);
                        idx2 = P->get_sortedIndex(k, j);
                    }
                }
            }
            coverCuts[k][0] = maxSizeClique;

            // print
            if (coverCuts[k].size() >= 2) {
                std::cout << "\n -> a new cut was generated on objective " << k << " : ";
                for (int l = 1; l < coverCuts[k].size() - 1; l++) {
                    std::cout << "x[" << coverCuts[k][l] << "] + ";
                }
                std::cout << "x[" << coverCuts[k][coverCuts[k].size() - 1] << "] <= " << maxSizeClique << "\n";
            }
            else {
                std::cout << "\n --> No cut generated !\n";
            }
        }

    }

}

/* \brief Compute the cover inequalities
 *
 * Compute cover inequalities for minimization problems, and add them to coverCuts, their description. The cover cut is computed by considering variables from the most to the least often fractional
 * in the extreme points of the lower bound set.
 * Note: assume WARMSTARTING of LP relax
 */
void Node::computeMofCoverCuts(BranchingDecisions& bd) {

    std::vector<int>* I = new std::vector<int>(P->get_n());
    SLUB s = SLUB(branchingDecision.slub);
    dynamic_cast<LinearRelaxation*>(LB)->computeFractionalProportion(I, s);

    int rhs = 0;
    for (int k = 0; k < P->get_p(); k++) {

        if (coverCuts.size() <= k) coverCuts.push_back(std::vector<int>(0));
        else coverCuts[k] = std::vector<int>(0);

        if (P->get_objDir(k) == -1 && bd.slub[k] < INT_MAX - 2 && iteration != 0) { // ... if this objective function is minimized and bounded by an objective branching constraint

            // compute the actual rhs by substracting variables fixed to 1 with a positive objective coefficient
            
            rhs = bd.slub[k];
            for (int i = 0; i < P->get_n(); i++) {
                if (bd.lb[i] == 1) {
                    rhs -= max(P->get_objective(k, i), 0);
                }
            }
            coverCuts[k].push_back(rhs);

            // compute the cut in function of order given by I, the proportion of fractional value of each variable

            std::vector<int> h = std::vector<int>(0);
            make_heap(h.begin(), h.end());
            int idx, maxSizeClique = -1, i = 0;
            
            while (rhs >= 0 && i < P->get_n()) {
                idx = (*I)[i];
                if (bd.lb[idx] != bd.ub[idx]) {
                    coverCuts[k].push_back(idx);
                    rhs -= P->get_objective(k, idx);
                    maxSizeClique++;
                    h.push_back(P->get_objective(k, idx));
                    push_heap(h.begin(), h.end());
                }
                i++;
            }

            if (i == P->get_n()) {
                maxSizeClique++;
            }
            else { // we try to improve the inequality by adding new variables to the cut
                int max = -1;
                while (i > 0 && i != P->get_n()) {
                    idx = (*I)[i];
                    if (bd.lb[idx] != bd.ub[idx]) {
                        max = h.front();
                        if (rhs + max - P->get_objective(k, idx) < 0) {
                            coverCuts[k].push_back(idx);
                            rhs = rhs + max - P->get_objective(k, idx);
                            pop_heap(h.begin(), h.end());
                            h.pop_back();
                            h.push_back(P->get_objective(k, idx));
                            push_heap(h.begin(), h.end());
                        }
                    }
                    i++;
                }
            }

            coverCuts[k][0] = maxSizeClique;

            // print
            if (coverCuts[k].size() >= 2) {
                std::cout << "\n -> a new cut was generated on objective " << k + 1 << " : ";
                for (int l = 1; l < coverCuts[k].size() - 1; l++) {
                    std::cout << "x[" << coverCuts[k][l] << "] + ";
                }
                std::cout << "x[" << coverCuts[k][coverCuts[k].size() - 1] << "] <= " << maxSizeClique << "\n";
            }
            else {
                std::cout << "\n --> No cut generated !\n";
            }
        }
    }
    delete I;
}

/* \brief Record at which iteration this node is explored.
 */
void Node::setIteration(int it) {
    iteration = it;
}

/*! \brief Try to fix integer bounds to each variable by solving LPs and looking at their feasibility
 *
 * \param BranchingDecisions newBd0. Describes a potential futur subproblem where a variable is set to 0.
 * \param BranchingDecisions newBd1. Describes a potential futur subproblem where a variable is set to 1.
 * \return true if this node is detected as infeasible after checking the LPs.
 */
bool Node::fixAndSelectBranchingVariable(BranchingDecisions* newBd0, BranchingDecisions* newBd1, UpperBoundSet* U) {

    bool feasible;
    bool nodeClosed = false;
    bool modified = false;
    bool fixToBound = false;
    int i = 0;

    while (!nodeClosed && i < P->get_n()) {

        if (newBd0->ub[i] != newBd0->lb[i]) {

            //std::cout << "\n -> we test x[" << i << "]\n";

            // try update bounds manually

            nodeClosed = manualInspection(newBd0, newBd1, i, &modified);

            //std::cout << "     . lb = " << newBd0->lb[i] << "\n";
            //std::cout << "     . ub = " << newBd0->ub[i] << "\n";

            if (!nodeClosed) { // go to next step if node not closed

                // test x = 0, if relevant


                if (newBd0->lb[i] != 1 && !dynamic_cast<LinearRelaxation*>(LB)->takeValue(newBd0, i, 0)) { // not fixed to 1 with preprocessing & value not taken by any of the extreme points

                    //std::cout << " testing x = 0...\n";
   
                    fixToBound = newBd0->ub[i] == 0;

                    newBd0->ub[i] = 0;
                    feasible = prob->solve(*newBd0);
                    if (!fixToBound) newBd0->ub[i] = 1;

                    if (!feasible) { // x = 0 not feasible => fix to 1

                        newBd0->lb[i] = 1;
                        newBd1->lb[i] = 1;
                        //std::cout << " x[" << i << "] fixed to 1 with LP\n";

                        // test x = 1, if relevant


                        if (newBd0->ub[i] == 0) { // x cannot be 1 -> node closed by infeasibility.
                            nodeClosed = true;
                            //std::cout << "    -> no integral value for x[" << i << "]\n";
                        }
                        //else if (!dynamic_cast<LinearRelaxation*>(LB)->takeValue(newBd0, i, 1)) { // x can be 1 -> x fixed to 1 => solve lp if not the case

                        //    std::cout << " testing x = 1...\n";
    
                        //    feasible = prob->solve(*newBd0);

                        //    if (!feasible) {
                        //        nodeClosed = true;
                        //        std::cout << "    -> no integral value for x[" << i << "]\n";
                        //    }
                        //}

                    }

                }

                // test x = 1, if x cannot be 0, and if relevant

                if (!nodeClosed && newBd0->ub[i] != 0 && !dynamic_cast<LinearRelaxation*>(LB)->takeValue(newBd0, i, 1)) { // else 

                    //std::cout << " testing x = 1...\n";

                    fixToBound = newBd0->lb[i] == 1;

                    newBd0->lb[i] = 1;
                    feasible = prob->solve(*newBd0);
                    if (!fixToBound) newBd0->lb[i] = 0;

                    if (!feasible) { // x = 1 not feasible => fix to 0

                        newBd0->ub[i] = 0;
                        newBd1->ub[i] = 0;
                        //std::cout << " x[" << i << "] fixed to 0 with LP\n";

                        if (newBd0->lb[i] == 1) { // x cannot be 0 -> node closed by infeasibility
                            nodeClosed = true;
                        }

                    }
                }
            }
        }

        i++;
    }


    // test fathoming by optimality

    bool allFixed = true;
    int id = 0;
    while (!nodeClosed && allFixed && id != P->get_n()) {
        if (newBd0->lb[id] != newBd0->ub[id]) {
            allFixed = false;
        }
        id++;
    }
    if (!nodeClosed && allFixed) {
        nodeClosed = true; // fathomed by optimality
        U->updateUB(newBd0);
        //std::cout << "\n     FATHOMED BY OPTIMALITY\n";
    }


    // select the branching variable

    if (!nodeClosed) {
        SLUB s = SLUB(newBd0->slub);
        int index = dynamic_cast<LinearRelaxation*>(LB)->computeMostOftenFractionalIndex(s, newBd0);

        newBd0->ub[index] = 0;
        newBd1->lb[index] = 1;
        newBd0->lastSplittedIndex = index;
        newBd1->lastSplittedIndex = index;
    }

    return nodeClosed;
}

/*! \brief Try to fix integer bounds to each variable by solving LPs and looking at their feasibility
 *
 * \param BranchingDecisions newBd0. Describes a potential futur subproblem where a variable is set to 0.
 * \param BranchingDecisions newBd1. Describes a potential futur subproblem where a variable is set to 1.
 * \return true if this node is detected as infeasible after checking the LPs.
 */
bool Node::fixAndSelectBranchingVariable1(BranchingDecisions* newBd0, BranchingDecisions* newBd1, UpperBoundSet* U) {

    std::vector<bool> free(P->get_n(), true);
    bool nodeClosed = false;
    bool modified = true;
    bool atLeastOneLpWasSolved = false;
    bool allInfinite = true;
    int id;

    //std::vector<std::vector<double>> SBscores(P->get_n(), std::vector<double>(2, 0));

    std::vector<double> sol0(P->get_n(), -1);
    std::vector<double> sol1(P->get_n(), -1);
    std::vector<double> reducedCosts(P->get_n(), -1);

    // apply the average normal vector to the WS of the probing model
    //std::vector<double> l = dynamic_cast<LinearRelaxation*>(LB)->getAverageNormalVector(newBd0); // everything is forced to 1 currently !! see function to change that
    std::vector<double> l(P->get_p(), 1);
    for (int k = 0; k < P->get_p(); k++) {
        if (newBd0->slub[k] >= 10000000) {
            l[k] = 0;
        }
        else {
            l[k] = 1;
            allInfinite = false;
        }
    }
    if (true || allInfinite) l = std::vector<double>(param->searchDir); // P->get_p(), 1
    /*std::cout << "\n l = ";
    for (int k = 0; k < P->get_p(); k++) {
        l[k] = round(10 * l[k]);
        std::cout << l[k] << " ";
    }
    std::cout << "\n";*/
    if (param->domiVarFix == 2) prob->updateWSCoef(P, l);

    // get free variables
    for (int i = 0; i < P->get_n(); i++) {
        if (newBd0->lb[i] == newBd0->ub[i]) {
            free[i] = false;
        }
    }

    // test whether OB constraints were generated
    bool OBgenerated = false;
    for (int k = 0; k < P->get_p(); k++) {
        if (newBd0->slub[k] != branchingDecision.slub[k]) {
            OBgenerated = true;
        }
    }

    // test preprocessing
    
    stat->cpuProbingManual.StartTimer();
    while (!nodeClosed && modified) {

        id = 0;
        modified = false;
        while (!nodeClosed && id < P->get_n()) {
            //std::cout << "\nWe try to fix x_" << id << "\n";
            if (free[id]) {
                nodeClosed = manualInspection(newBd0, newBd1, id, &modified);
            }
            else {
                //std::cout << " is already fixed lol.\n";
            }
            //std::cout << " node closed : " << nodeClosed << "\n";
            id++;
        }

    }
    stat->cpuProbingManual.StopTimer();

    // test LPs

    bool feasible, dominated{}, dominanceOccured = false;
    bool wasFixed = false;
    double objSol0 = -1, objSol1 = -1;
    double maxWsVal = UB->getLargestWeightedSumValue(newBd0, l);
    bool feasibleForSubpb;
    int worstVar = -1000000, minWsVal; // min ws val for a variable for each bound, and min ws val considering variable branched on.
    bool bothLpSolved, lp0Solved, lp1solved;
    id = 0;
    //std::cout << " ws value: " << maxWsVal << "\n";

    while (OBgenerated && !nodeClosed && id < P->get_n()) {

        wasFixed = false;

        lp0Solved = lp1solved = false;
        if (free[id]) {

            // test x_i = 0 if not proven impossible beforehand

            feasibleForSubpb = sol0[id] != 0 && sol1[id] != 0;
            if (newBd0->lb[id] != 1 && feasibleForSubpb) { // toTest[id][0] // dynamic_cast<LinearRelaxation*>(LB)->takeValue(newBd0, id, 0)
                if (newBd0->ub[id] == 0) {
                    wasFixed = true;
                    //for (int i = 0; i < P->get_n(); i++) sol1[i] = -1; // we won't test x = 1 -> it's not up to date anymore
                }
                else {
                    newBd0->ub[id] = 0;
                    newBd1->ub[id] = 0;
                }
                stat->cpuProbingCplex.StartTimer();
                feasible = prob->solve(*newBd0);
                stat->cpuProbingCplex.StopTimer();
                stat->lpProbing++;
                if (!feasible) { // x_i = 0 is not possible

                    //std::cout << "   -> x[" << id << "] fixed to 1 with LP\n";
                    stat->nbFixingLp++;
                    if (wasFixed) { // x_i = 1 proven impossible before -> we close the node by infeasibility
                        nodeClosed = true;
                        //if (iteration > 3590) std::cout << "\n     FATHOMED BY INFEASIBILITY\n";
                    }
                    else { // otherwise we fix the variable to 1
                        newBd0->ub[id] = 1;
                        newBd1->ub[id] = 1;
                        newBd0->lb[id] = 1;
                        newBd1->lb[id] = 1;
                        nodeClosed = manualInspection(newBd0, newBd1, id, &modified);
                    }
                }
                else {
                    prob->getSolution(&sol0);
                    //prob->getReducedCosts(&reducedCosts);
                    objSol0 = prob->retrieveObjectiveValue();
                    atLeastOneLpWasSolved = true;
                    lp0Solved = true;
                    stat->timeUpdateUB.StartTimer();
                    U->updateUB(sol0);
                    stat->timeUpdateUB.StopTimer();
                    //std::cout << " obj val 0 : " << objSol0 << "\n";

                    // proceed to dominance check

                    //dominated = ACTIVATE_DOMINANCE_VARIABLE_FIXING && objSol0 > maxWsVal; // UB->testWeightedSumValue(&branchingDecision, l, objSol0);
                    dominated = param->domiVarFix == 2 && objSol0 > maxWsVal; // UB->testWeightedSumValue(&branchingDecision, l, objSol0);
                    if (dominated) {
                        //std::cout << "   -> x[" << id << "] fixed to 1 by dominance\n";
                        if (wasFixed) { // x_i = 1 proven impossible before -> we close the node by infeasibility
                            nodeClosed = true;
                            //std::cout << "\n     FATHOMED BY DOMINANCE\n";
                        }
                        else {
                            newBd0->ub[id] = 1;
                            newBd1->ub[id] = 1;
                            newBd0->lb[id] = 1;
                            newBd1->lb[id] = 1;
                            nodeClosed = manualInspection(newBd0, newBd1, id, &modified);
                        }
                    }
                    

                    if (dominated && wasFixed) {
                        nodeClosed = true;
                    }
                    else if (!dominated && !wasFixed) { // we reset the bounds in case x_i = 1 is a possible value
                        newBd0->ub[id] = 1;
                        newBd1->ub[id] = 1;
                    }
                }
            }

            // test x_i = 1 if not proven impossible beforehand


            feasibleForSubpb = sol0[id] != 1 && sol1[id] != 1;
            if (newBd0->ub[id] != 0 && feasibleForSubpb) { // dynamic_cast<LinearRelaxation*>(LB)->takeValue(newBd0, id, 1)
                if (newBd0->lb[id] == 1) {
                    wasFixed = true;
                    //for (int i = 0; i < P->get_n(); i++) sol0[i] = -1; // we won't test x = 0 -> it's not up to date anymore
                }
                else {
                    newBd0->lb[id] = 1;
                    newBd1->lb[id] = 1;
                }
                stat->cpuProbingCplex.StartTimer();
                feasible = prob->solve(*newBd0);
                stat->cpuProbingCplex.StopTimer();
                stat->lpProbing++;
                if (!feasible) { // x_i = 1 is not possible

                    //std::cout << "   -> x[" << id << "] fixed to 0 with LP\n";
                    stat->nbFixingLp++;
                    if (wasFixed) { // x_i = 0 proven impossible before -> we close the node by infeasibility
                        nodeClosed = true;
                        //if (iteration > 3590) std::cout << "\n     FATHOMED BY INFEASIBILITY\n";
                    }
                    else { // otherwise we fix the variable to 0
                        newBd0->ub[id] = 0;
                        newBd1->ub[id] = 0;
                        newBd0->lb[id] = 0;
                        newBd1->lb[id] = 0;
                        nodeClosed = manualInspection(newBd0, newBd1, id, &modified);
                    }
                }
                else {
                    prob->getSolution(&sol1);
                    //prob->getReducedCosts(&reducedCosts);
                    objSol1 = prob->retrieveObjectiveValue();
                    atLeastOneLpWasSolved = true;
                    stat->timeUpdateUB.StartTimer();
                    U->updateUB(sol1);
                    stat->timeUpdateUB.StopTimer();
                    lp1solved = true;
                    //std::cout << " obj val 1 : " << objSol1 << "\n";

                    // proceed to dominance check

                    //dominated = ACTIVATE_DOMINANCE_VARIABLE_FIXING && objSol1 > maxWsVal; // UB->testWeightedSumValue(&branchingDecision, l, objSol1);
                    dominated = param->domiVarFix == 2 && objSol1 > maxWsVal; // UB->testWeightedSumValue(&branchingDecision, l, objSol1);
                    if (dominated) {
                        //std::cout << "   -> x[" << id << "] fixed to 0 by dominance\n";
                        if (wasFixed) { // x_i = 0 proven impossible before -> we close the node by dominance
                            nodeClosed = true;
                            //std::cout << "\n     FATHOMED BY DOMINANCE\n";
                        }
                        else { // otherwise we fix the variable to 0
                            newBd0->ub[id] = 0;
                            newBd1->ub[id] = 0;
                            newBd0->lb[id] = 0;
                            newBd1->lb[id] = 0;
                            nodeClosed = manualInspection(newBd0, newBd1, id, &modified);
                        }
                    }

                    if (dominated && wasFixed) {
                        nodeClosed = true;
                    }
                    if (!dominated && !wasFixed) { // we reset the bounds in case x_i = 1 is a possible value
                        newBd0->lb[id] = 0;
                        newBd1->lb[id] = 0;
                    }
                }
            }

            //if (lp0Solved && lp1solved) {
                //std::cout << "subpb: " << objSol0 << " and " << objSol1 << "\n";
            //}
        }
        id++;
    }

    // test fathoming by optimality

    bool allFixed = true;
    id = 0;
    while (!nodeClosed && allFixed && id != P->get_n()) {
        if (newBd0->lb[id] != newBd0->ub[id]) {
            allFixed = false;
        }
        id++;
    }
    if (!nodeClosed && allFixed) {
        nodeClosed = true; // fathomed by optimality
        U->updateUB(newBd0);
        //if (iteration > 3590) std::cout << "\n     FATHOMED BY OPTIMALITY\n";
    }

    // if the node is not fathomed, choose a free variable to branch on (rule: MOST_OFTEN_FRACTIONAL)

    if (!nodeClosed) {
        SLUB s = SLUB(newBd0->slub);
        int index = 0;
        index = dynamic_cast<LinearRelaxation*>(LB)->computeMostOftenFractionalIndex(s, newBd0);
        //index = getFullStrongBranchingWSIndex(newBd0);
        //index = dynamic_cast<LinearRelaxation*>(LB)->computeMostViolatingIndex(newBd0);

        //if (OBgenerated) index = P->getLargestFreeCoef(newBd0);
        //else index = dynamic_cast<LinearRelaxation*>(LB)->computeMostOftenFractionalIndex(s, newBd0);

        //if (atLeastOneLpWasSolved) index = getLargestFractionalReducedCost(newBd0, reducedCosts, sol0, sol1);
        //else index = dynamic_cast<LinearRelaxation*>(LB)->computeMostOftenFractionalIndex(s, newBd0);

        newBd0->ub[index] = 0;
        newBd1->lb[index] = 1;
        newBd0->lastSplittedIndex = index;
        newBd1->lastSplittedIndex = index;
    }

    return nodeClosed;
}

/*! \brief Try to fix integer bounds to each variable by solving LPs and looking at their feasibility
 *
 * \param BranchingDecisions newBd0. Describes a potential futur subproblem where a variable is set to 0.
 * \param BranchingDecisions newBd1. Describes a potential futur subproblem where a variable is set to 1.
 * \return true if this node is detected as infeasible after checking the LPs.
 */
bool Node::fixAndSelectBranchingVariable2(BranchingDecisions* newBd0, BranchingDecisions* newBd1, UpperBoundSet* U) {

    bool modificationOccured = true;
    bool varFixed;
    bool nodeClosed = false;
    bool feasible0 = true, feasible1 = true;
    bool wasFixed = false;

    std::vector<std::vector<bool>> toTest(P->get_n());
    std::vector<bool> inList(P->get_n());
    std::list<int> L(0); // contains variables to explore next

    // initialise lists

    for (int i = 0; i < P->get_n(); i++) {
        toTest[i] = std::vector<bool>(2, true);
        inList[i] = false;
        if (newBd0->lb[i] == newBd0->ub[i]) { // if the variable is fixed, we don't check it
            toTest[i][0] = false;
            toTest[i][1] = false;
            //std::cout << " x[" << i << "] is fixed with branching\n";
        }
    }
    if (param->LBset == LP_RELAX || param->LBset == WARMSTARTED_LP_RELAX) dynamic_cast<LinearRelaxation*>(LB)->checkFeasibleValues(newBd0, &toTest); // check values of extreme points in the LB set
    else throw std::string("Error: lower bound set configuration not supported for VARIABLE_FIXING");

    // try to fix variables -> change that to list exploration

    int id = 0;
    while (id < P->get_n()) { // L.size() == 0 &&  // find first index to explore
        if (toTest[id][0] || toTest[id][1]) {
            L.push_back(id);
            inList[id] = true;
            std::cout << " x[" << id << "] added to the exploration list\n";
        }
        id++;
    }

    while (!nodeClosed && L.size() != 0) {

        varFixed = false;

        id = L.front();
        L.pop_front();
        inList[id] = false; // let to true to explore each variable only once.
        std::cout << "\nWe try to fix x_" << id << "\n";

        //std::cout << "before manInspec: x9 in [" << newBd0->lb[9] << " , " << newBd0->ub[9] << "]\n";
        nodeClosed = manualInspection(newBd0, newBd1, &toTest, id, &inList, &L);
        //std::cout << "after manInspec: x9 in [" << newBd0->lb[9] << " , " << newBd0->ub[9] << "]\n";

        // LB set test if variable fixed
        if (newBd0->lb[id] == 1) dynamic_cast<LinearRelaxation*>(LB)->updateFeasibleValues(&toTest, newBd0, id, 1);
        else if (newBd0->ub[id] == 0) dynamic_cast<LinearRelaxation*>(LB)->updateFeasibleValues(&toTest, newBd0, id, 0);

        // LP test with x_id = 0
        if (toTest[id][0]) {
            if (newBd0->ub[id] == 0) wasFixed = true;
            newBd0->ub[id] = 0;
            newBd1->ub[id] = 0;
            feasible0 = prob->solve(*newBd0);
            if (!wasFixed) {
                newBd0->ub[id] = 1;
                newBd1->ub[id] = 1;
            }

            std::cout << "An LP is solved for x = 0\n";

            if (!feasible0) {
                varFixed = true;
                newBd0->lb[id] = 1;
                newBd1->lb[id] = 1;
                std::cout << "   -> x[" << id << "] fixed to 1 with LP\n";
                toTest[id][0] = false;
                dynamic_cast<LinearRelaxation*>(LB)->updateFeasibleValues(&toTest, newBd0, id, 1);
                updateListVarToFix(&toTest, &inList, &L, id);
            }
        }
        else {
            if (newBd0->lb[id] == 1) feasible0 = false;
            else feasible0 = true; // just to make sure that the next loop is correct when fixing variables
        }

        // LP test with x_id = 1
        if (toTest[id][1]) {
            //std::cout << " ub = " << newBd0->ub[id] << "\n";
            if (newBd0->ub[id] == 0) {
                std::cout << " AH!\n";
            }
            newBd0->lb[id] = 1;
            newBd1->lb[id] = 1;
            feasible1 = prob->solve(*newBd0);
            std::cout << "An LP is solved for x = 1\n";
            if (feasible0) { // if x_id = 0 was feasible, we reset the lower bounds to 0
                newBd0->lb[id] = 0;
                newBd1->lb[id] = 0;
            }
            if (!feasible1) {
                varFixed = true;
                newBd0->ub[id] = 0;
                newBd1->ub[id] = 0;
                std::cout << "   -> var[" << id << "] fixed to 0 with LP\n";
                toTest[id][1] = false;
                dynamic_cast<LinearRelaxation*>(LB)->updateFeasibleValues(&toTest, newBd0, id, 0);
                updateListVarToFix(&toTest, &inList, &L, id);
            }
        }

        // test for infeasibility
        if (newBd0->lb[id] == 1 && newBd0->ub[id] == 0) {
            nodeClosed = true; // fathomed by infeasibility
            //std::cout << "\n     FATHOMED BY INFEASIBILITY\n";
        }
    }

    // test for optimality
    bool allFixed = true;
    id = 0;
    while (allFixed && id != P->get_n()) {
        if (newBd0->lb[id] != newBd0->ub[id]) {
            allFixed = false;
        }
        id++;
    }
    if (allFixed) {
        nodeClosed = true; // fathomed by optimality
        U->updateUB(newBd0);
        //std::cout << "\n     FATHOMED BY OPTIMALITY\n";
    }

    // choose a free variable to branch on (rule: MOST_OFTEN_FRACTIONAL)

    if (!nodeClosed) {
        SLUB s = SLUB(newBd0->slub);
        int index = dynamic_cast<LinearRelaxation*>(LB)->computeMostOftenFractionalIndex(s, newBd0);

        newBd0->ub[index] = 0;
        newBd1->lb[index] = 1;
        newBd0->lastSplittedIndex = index;
        newBd1->lastSplittedIndex = index;
    }


    return nodeClosed;
}

/* \brief Tries to fix variables manually by looking at the constraints
 *
 * \param BranchingDecisions newBd0. Describes a potential futur subproblem where a variable is set to 0.
 * \param BranchingDecisions newBd1. Describes a potential futur subproblem where a variable is set to 1.
 */
bool Node::manualInspection(BranchingDecisions* newBd0, BranchingDecisions* newBd1, std::vector<std::vector<bool>>* toTest) {

    bool modificationOccured = false, internalModificationOccured = true;
    int lhs, rhs;
    std::vector<bool> free(P->get_n(), true);

    // find free variables

    for (int i = 0; i < P->get_n(); i++) {
        if (newBd0->lb[i] == newBd0->ub[i]) 
            free[i] = false;
    }

    // search in the regular constraints

    while (internalModificationOccured) {

        internalModificationOccured = false;

        for (int j = 0; j < P->get_m(); j++) {

            // test for <= constraints (also includes == constraints)
            if (P->get_signCte(j) == 1 || P->get_signCte(j) == 2) {

                // compute the right-hand side of the constraint, considering already fixed variables
                rhs = P->get_rhs(j);
                for (int i = 0; i < P->get_n(); i++) {
                    if (!free[i] && newBd0->lb[i] == 1) // variable fixed to 1
                        rhs -= P->get_constraint(j, i);
                    else if (free[i] && P->get_constraint(j, i) < 0) // free variable with negative coefficient -> considered as fixed to 1 here
                        rhs -= P->get_constraint(j, i);
                }

                // test variables
                for (int i = 0; i < P->get_n(); i++) {
                    if (free[i] && P->get_constraint(j, i) > rhs) { // fix variable to 0: rhs violated
                        //printCteFixed(free, j);
                        newBd0->ub[i] = 0;
                        newBd1->ub[i] = 0;
                        free[i] = false;
                        modificationOccured = true;
                        internalModificationOccured = true;
                        //std::cout << "x[" << i << "] fixed to 0 in <= or == cte\n";
                    }
                }

                // compute sum of negative coefficients & actual rhs
                rhs = P->get_rhs(j);
                lhs = 0;
                for (int i = 0; i < P->get_n(); i++) {
                    if (!free[i] && newBd0->lb[i] == 1) rhs -= P->get_constraint(j, i);
                    if (free[i] && P->get_constraint(j, i) < 0) lhs += P->get_constraint(j, i);
                }
                if (lhs == rhs) {
                    for (int i = 0; i < P->get_n(); i++) {
                        if (free[i] && P->get_constraint(j, i) < 0) {
                            //printCteFixed(free, j);
                            newBd0->lb[i] = 1;
                            newBd1->lb[i] = 1;
                            free[i] = false;
                            modificationOccured = true;
                            internalModificationOccured = true;
                            //std::cout << "x[" << i << "] fixed to 1 in <= or == cte\n";
                        }
                    }
                }
            }

            // test for >= constraints (also includes == constraints)
            if (P->get_signCte(j) == 0 || P->get_signCte(j) == 2) {

                // compute the right-hand side of the constraint, considering already fixed variables
                rhs = P->get_rhs(j);
                for (int i = 0; i < P->get_n(); i++) {
                    if (!free[i] && newBd0->lb[i] == 1) { // variable fixed to 1
                        rhs -= P->get_constraint(j, i);
                    }
                }

                // add all positive coefficient and check with rhs
                lhs = 0;
                for (int i = 0; i < P->get_n(); i++) {
                    if (free[i] && P->get_constraint(j, i) > 0) {
                        lhs += P->get_constraint(j, i);
                    }
                }
                if (lhs == rhs) {
                    for (int i = 0; i < P->get_n(); i++) {
                        if (free[i] && P->get_constraint(j, i) > 0) { // fix to 1
                            //printCteFixed(free, j);
                            newBd0->lb[i] = 1;
                            newBd1->lb[i] = 1;
                            free[i] = false;
                            modificationOccured = true;
                            internalModificationOccured = true;
                            //std::cout << "x[" << i << "] fixed to 1 in >= or == cte\n";
                        }
                        else if (free[i] && P->get_constraint(j, i) < 0) { // fix to 0
                            //printCteFixed(free, j);
                            newBd0->ub[i] = 0;
                            newBd0->lb[i] = 0;
                            free[i] = false;
                            internalModificationOccured = true;
                            //std::cout << "x[" << i << "] fixed to 0 in >= or == cte\n";
                        }
                    }
                }
                else if (lhs < rhs) {
                    std::cout << "we need a feasbility check!\n";
                }

            }

            // check feasibility final fix
        }

    }

    // update lists

    for (int i = 0; i < P->get_n(); i++) {
        if (newBd0->lb[i] == newBd0->ub[i]) { // if the variable is fixed, we don't check it
            (*toTest)[i][0] = false;
            (*toTest)[i][1] = false;
        }
    }

    return modificationOccured;
}

/* \brief Tries to fix variable i manually by looking at the constraints
 *
 * \param BranchingDecisions newBd0. Describes a potential futur subproblem where a variable is set to 0.
 * \param BranchingDecisions newBd1. Describes a potential futur subproblem where a variable is set to 1.
 * \param int i. Index of the variable to fix.
 */
bool Node::manualInspection(BranchingDecisions* newBd0, BranchingDecisions* newBd1, std::vector<std::vector<bool>>* toTest, int i, std::vector<bool>* inList, std::list<int>* L) {

    bool closeNode = false;
    bool variableFixed = false;
    int lhs, rhs;
    int oldLb = newBd0->lb[i];
    int oldUb = newBd0->ub[i];
    std::vector<bool> free(P->get_n(), true);

    // find free variables

    for (int i = 0; i < P->get_n(); i++) {
        if (newBd0->lb[i] == newBd0->ub[i])
            free[i] = false;
    }

    // fix variable x_i

    for (int j = 0; j < P->get_m(); j++) {

        if (P->get_constraint(j, i) != 0) { // if x_i is in constraint j

            //std::cout << " ... checking x_" << i << " in constraint " << j + 1 << "\n";
            

            // test for <= constraints (also includes == constraints)
            if (P->get_signCte(j) == 1 || P->get_signCte(j) == 2) {

                rhs = P->get_rhs(j);
                for (int l = 0; l < P->get_n(); l++) {
                    if (!free[l] && newBd0->lb[l] == 1) rhs -= P->get_constraint(j, l);
                    else if (free[l] && P->get_constraint(j, l) < 0) rhs -= P->get_constraint(j, l);
                }

                if (P->get_constraint(j, i) > 0) {
                    //std::cout << " ub = " << rhs << " / " << P->get_constraint(j, i) << " = " << rhs / P->get_constraint(j, i) << "\n";
                    newBd0->ub[i] = min(max(0, min(1, trunc(rhs / P->get_constraint(j, i)))), newBd0->ub[i]); // remove min(1, .) for general integer?
                    newBd1->ub[i] = min(max(0, min(1, trunc(rhs / P->get_constraint(j, i)))), newBd0->ub[i]); // remove min(1, .) for general integer?
                }
                else {
                    //std::cout << " lb = " << rhs << " / " << P->get_constraint(j, i) << " = " << rhs / P->get_constraint(j, i) << "\n";
                    newBd0->lb[i] = max(min(1, max(0, ceil(rhs / P->get_constraint(j, i)))), newBd0->lb[i]); // remove min(1, .) for general integer?
                    newBd1->lb[i] = max(min(1, max(0, ceil(rhs / P->get_constraint(j, i)))), newBd0->lb[i]); // remove min(1, .) for general integer?
                }
            }

            // test for >= constraints (also includes == constraints)
            if (P->get_signCte(j) == 0 || P->get_signCte(j) == 2) {

                rhs = P->get_rhs(j);
                for (int l = 0; l < P->get_n(); l++) {
                    if (!free[l] && newBd0->lb[l] == 1) rhs -= P->get_constraint(j, l);
                    else if (free[l] && P->get_constraint(j, l) > 0) rhs -= P->get_constraint(j, l);
                }

                if (P->get_constraint(j, i) > 0) {
                    //std::cout << " lb = " << rhs << " / " << P->get_constraint(j, i) << " = " << rhs / P->get_constraint(j, i) << "\n";
                    newBd0->lb[i] = max(min(1, max(0, ceil(rhs / P->get_constraint(j, i)))), newBd0->lb[i]);
                    newBd1->lb[i] = max(min(1, max(0, ceil(rhs / P->get_constraint(j, i)))), newBd0->lb[i]);
                }
                else {
                    //std::cout << " ub = " << rhs << " / " << P->get_constraint(j, i) << " = " << rhs / P->get_constraint(j, i) << "\n";
                    newBd0->ub[i] = min(max(0, min(1, trunc(rhs / P->get_constraint(j, i)))), newBd0->ub[i]);
                    newBd1->ub[i] = min(max(0, min(1, trunc(rhs / P->get_constraint(j, i)))), newBd0->ub[i]);
                }

            }
        }
    }

    // feasibility test

    //std::cout << " ub old vs new : " << oldUb << " vs " << newBd0->ub[i] << "\n";

    if (newBd0->lb[i] > newBd0->ub[i]) {
        closeNode = true; // fathomed by infeasibility
        (*toTest)[i][0] = false;
        (*toTest)[i][1] = false;
    }
    else if (newBd0->lb[i] == 1 && oldLb == 0) {
        (*toTest)[i][0] = false;
        variableFixed = true;
        std::cout << "    -> var[" << i << "] fixed to 1 with preprocessing\n";
    }
    else if (newBd0->ub[i] == 0 && oldUb == 1) {
        (*toTest)[i][1] = false;
        variableFixed = true;
        std::cout << "    -> var[" << i << "] fixed to 0 with preprocessing\n";
    }

    // if x_i was fixed, we have to explore variables in the same constraints as x_i

    if (variableFixed) {    
        updateListVarToFix(toTest, inList, L, i);
    }

    return closeNode;
}

/* \brief Tries to fix variable i manually by looking at the constraints
 *
 * \param BranchingDecisions newBd0. Describes a potential futur subproblem where a variable is set to 0.
 * \param BranchingDecisions newBd1. Describes a potential futur subproblem where a variable is set to 1.
 * \param int i. Index of the variable to fix.
 */
bool Node::manualInspection(BranchingDecisions* newBd0, BranchingDecisions* newBd1, int i, bool* modified) {

    bool closeNode = false;
    bool variableFixed = false;
    int lhs, rhs;
    int oldLb = newBd0->lb[i];
    int oldUb = newBd0->ub[i];
    std::vector<bool> free(P->get_n(), true);


    // find free variables

    for (int i = 0; i < P->get_n(); i++) {
        if (newBd0->lb[i] == newBd0->ub[i])
            free[i] = false;
    }


    // fix variable x_i in initial constraints

    for (int j = 0; j < P->get_m(); j++) {

        if (P->get_constraint(j, i) != 0) { // if x_i is in constraint j

            //std::cout << " ... checking x_" << i << " in constraint " << j + 1 << "\n";

            // test for <= constraints (also includes == constraints)
            if (P->get_signCte(j) == 1 || P->get_signCte(j) == 2) {

                //std::cout << " test 1 -> ";

                rhs = P->get_rhs(j);
                for (int l = 0; l < P->get_n(); l++) {
                    if (l != i && !free[l] && newBd0->lb[l] == 1) rhs -= P->get_constraint(j, l);
                    else if (l != i && free[l] && P->get_constraint(j, l) < 0) rhs -= P->get_constraint(j, l);
                }

                if (P->get_constraint(j, i) > 0) {
                    //std::cout << " ub = " << rhs << " / " << P->get_constraint(j, i) << " = " << rhs / P->get_constraint(j, i) << "\n";
                    newBd0->ub[i] = min(max(0, min(1, trunc(rhs / P->get_constraint(j, i)))), newBd0->ub[i]); // remove min(1, .) for general integer?
                    newBd1->ub[i] = min(max(0, min(1, trunc(rhs / P->get_constraint(j, i)))), newBd0->ub[i]); // remove min(1, .) for general integer?
                }
                else {
                    //std::cout << " lb = " << rhs << " / " << P->get_constraint(j, i) << " = " << rhs / P->get_constraint(j, i) << "\n";
                    newBd0->lb[i] = max(min(1, max(0, ceil(rhs / P->get_constraint(j, i)))), newBd0->lb[i]); // remove min(1, .) for general integer?
                    newBd1->lb[i] = max(min(1, max(0, ceil(rhs / P->get_constraint(j, i)))), newBd0->lb[i]); // remove min(1, .) for general integer?
                }
            }

            // test for >= constraints (also includes == constraints)
            if (P->get_signCte(j) == 0 || P->get_signCte(j) == 2) {

                //std::cout << " test 2 -> ";

                rhs = P->get_rhs(j);
                for (int l = 0; l < P->get_n(); l++) {
                    if (l != i && !free[l] && newBd0->lb[l] == 1) rhs -= P->get_constraint(j, l);
                    else if (l != i && free[l] && P->get_constraint(j, l) > 0) rhs -= P->get_constraint(j, l);
                }

                if (P->get_constraint(j, i) > 0) {
                    //std::cout << " lb = " << rhs << " / " << P->get_constraint(j, i) << " = " << rhs / P->get_constraint(j, i) << "\n";
                    newBd0->lb[i] = max(min(1, max(0, ceil(rhs / P->get_constraint(j, i)))), newBd0->lb[i]);
                    newBd1->lb[i] = max(min(1, max(0, ceil(rhs / P->get_constraint(j, i)))), newBd0->lb[i]);
                }
                else {
                    //std::cout << " ub = " << rhs << " / " << P->get_constraint(j, i) << " = " << rhs / P->get_constraint(j, i) << "\n";
                    newBd0->ub[i] = min(max(0, min(1, trunc(rhs / P->get_constraint(j, i)))), newBd0->ub[i]);
                    newBd1->ub[i] = min(max(0, min(1, trunc(rhs / P->get_constraint(j, i)))), newBd0->ub[i]);
                }

            }
        }
    }


    // fix variable x_i in objective branching constraints

    for (int k = 0; k < P->get_p(); k++) { // assume all coefficients of the same sign

        if (P->get_objDir(k) == -1) { // minimization objective
            
            rhs = 0;
            for (int l = 0; l < P->get_n(); l++) {
                if (l != i && newBd0->lb[l] == 1) rhs += P->get_objective(k, l);
            }

            //std::cout << " coef vs rhs: " << P->get_objective(k, i) << " vs " << rhs << "\n";
            if (rhs >= newBd0->slub[k]) {
                closeNode = true;
                //std::cout << " OB infeasibility !!\n";
            }
            else if (rhs + P->get_objective(k, i) >= newBd0->slub[k]) {
                newBd0->ub[i] = 0;
                newBd1->ub[i] = 0;
            }
        }
        else {

            rhs = 0;
            //std::cout << "\nrhs = ";
            for (int l = 0; l < P->get_n(); l++) {
                if (l != i && newBd0->ub[l] != 0) {
                    rhs += P->get_objective(k, l);
                    //std::cout << P->get_objective(k, l) << " + ";
                }
            }
            //std::cout << "\n";

            //std::cout << "coef = " << P->get_objective(k, i) << "\n";
            //std::cout << "OB = " << newBd0->slub[k] << "\n";

            if (rhs + P->get_objective(k, i) > newBd0->slub[k]) { // assume negative coefficients only??
                closeNode = true;
                //std::cout << " OB infeasibility !!\n";
            }
            else if (rhs > newBd0->slub[k]) {
                newBd0->lb[i] = 1;
                newBd1->lb[i] = 1;
            }
        }
    }


    // feasibility test

    //std::cout << " ub old vs new : " << oldUb << " vs " << newBd0->ub[i] << "\n";
    //std::cout << " lb vs ub : " << newBd0->lb[i]  << " vs " << newBd0->ub[i] << "\n";

    if (newBd0->lb[i] > newBd0->ub[i]) {
        closeNode = true; // fathomed by infeasibility
        //std::cout << "bruh\n";
        //if (iteration > 3590) std::cout << "\n     FATHOMED BY INFEASIBILITY\n";
    }
    else if (newBd0->lb[i] == 1 && oldLb == 0) {
        variableFixed = true;
        *modified = true;
        //if (iteration > -1) std::cout << "    -> var[" << i << "] fixed to 1 with preprocessing\n";
    }
    else if (newBd0->ub[i] == 0 && oldUb == 1) {
        variableFixed = true;
        *modified = true;
        //if (iteration > -1) std::cout << "    -> var[" << i << "] fixed to 0 with preprocessing\n";
    }

    return closeNode;
}

/* \brief Update the list of variable to potentially fix (L) after variable x_i was fixed in the model.
 *
 * \param vector of vector of bool, toTest. Used to know which variables should be tested.
 * \param vector of bool, inList. Vector that states whether each variable is already in the list of variables to explore or not. Used to avoid redundancies.
 * \param list of int. List of variables to explore for fixing.
 * \param int i. Index of the variable fixed.
 */
void Node::updateListVarToFix(std::vector<std::vector<bool>>* toTest, std::vector<bool>* inList, std::list<int>* L, int i) {

    for (int j = 0; j < P->get_m(); j++) {
        if (P->get_constraint(j, i) != 0) {
            for (int i2 = 0; i2 != P->get_n(); i2++) {
                if (i2 != i && P->get_constraint(j, i2) != 0 && !(*inList)[i2] && ((*toTest)[i2][0] || (*toTest)[i2][1])) {
                    L->push_back(i2);
                    (*inList)[i2] = true;
                }
            }
        }
    }
}

/* \brief Compute and returns the pseudo-cost of the last branching variable in the current node.
 * Note: returns ceta in Achteberg.
 */
double Node::getBranchingVariableScore() {

    double f = abs(branchingDecision.lb[branchingDecision.lastSplittedIndex] - sXiFatherNode);
    double delta = sLbCurrentNode - sLbFatherNode; // Comptue sLbCurrentNode !!

    return delta / f;
}

/* \brief returns the last branching variable.
 */
int Node::getLastBranchingVariable() {
    return branchingDecision.lastSplittedIndex;
}

/* \brief returns the last branching value.
     * Note: valid for binary problems only
     */
int getLastBranchingValue() {
    return -2;
}

/* \brief Returns the agregated value of variable i in the current node.
 */
double Node::getVariableScoreValue(int i) {
    return sXiCurrentNode[i];
}

/* Get the percentage of integer extreme point in the lower bound set computed at this node.
 *
 * \return the percentage as a double between 0 and 1.
 */
double Node::getPercentageIntegralityLB() {

    double pct = -2;

    if (iteration == 0) {
        pct = 0;
    }
    else {
        if (param->LBset == LP_RELAX || param->LBset == WARMSTARTED_LP_RELAX || param->LBset == PRECOMPUTED_WARMSTARTED_LP_RELAX) {
            pct = dynamic_cast<LinearRelaxation*>(LB)->getPercentageIntegrality();
        }
        else {
            throw std::string("Error: impossible to retreive the percentage of integrality of the linear relaxation. Invalid lower bound set.");
        }
    }  

    return pct;
}

/* Get the percentage of extreme point in the lower bound set computed at this node that satisfy the branching decisions of this node.
 *
 * \return the percentage as a double between 0 and 1.
 */
double Node::getPercentageFeasibilityLB() {

    double pct = -2;

    if (param->LBset == LP_RELAX || param->LBset == WARMSTARTED_LP_RELAX) {
        pct = dynamic_cast<LinearRelaxation*>(LB)->getPercentageFeasibility();
    }
    else {
        throw std::string("Error: impossible to retreive the percentage of integrality of the linear relaxation. Invalid lower bound set.");
    }

    return pct;
}

/*! \brief Compute the rhs of the weighted sum scalarization when the node is created.
 *
 * Computed at the creation of a node.
 */
void Node::computeWsScore() {

    stat->timeNodeSel.StartTimer();

    // compute (1,...,1) weights

    std::vector<double> l;

    //if (ADJUST_SEARCH_DIR_TO_COEF) l = std::vector<double>(param->searchDir);
    if (param->adjustBBWSbounds) l = std::vector<double>(param->searchDir);
    else l = std::vector<double>(P->get_p(), 1.0);

    if (!ENABLE_LB_SEARCH_IN_BEST_BOUND || param->LBset == WARMSTARTED_LP_RELAX || param->LBset == LP_RELAX) { // solve the lp to get ws value
        ws->adjustBounds(branchingDecision);
        bool solved = ws->solve(*P, l);

        if (solved) {
            score = ws->retrieveObjectiveValue(); // *P, l
            //std::cout << "we set score to " << score << "\n";
        }
        else {
            score = 0;
            //std::cout << "we created an infreasible node\n";
        }
    }
    else if (param->LBset == PRECOMPUTED_WARMSTARTED_LP_RELAX) { // the lp relax is already computed => search for extr pts with the smallest w.s. value
        score = dynamic_cast<LinearRelaxation*>(LB)->getSmallestWeightedSumValue(l);
    }

    stat->timeNodeSel.StopTimer();
}

/*! \brief Compute the the minimum difference between the rhs of an hyperplane and the w.s. of an lub using the normal vector of the hyperplane,
 * for each Lub. Then, return the maximal value obtained.
 *
 * Note: valid with PRECOMPUTED_WARMSTARTED_LP_RELAX only.
 */
void Node::computeMaxMinGapScore() {

    if (param->LBset == PRECOMPUTED_WARMSTARTED_LP_RELAX) {

        double maxGapAll = -1000000, minGapLub;
        std::list<LocalUpperBound>* lubs = UB->getLubs();
        std::list<LocalUpperBound>::iterator u;
        SLUB s = SLUB(branchingDecision.slub);

        for (u = lubs->begin(); u != lubs->end(); u++) {
            if (s.dominated(*u)) {
                minGapLub = dynamic_cast<LinearRelaxation*>(LB)->getMinFacetGap(*u);
                if (minGapLub > maxGapAll) {
                    maxGapAll = minGapLub;
                }
            }
        }

        score = maxGapAll;
    }
    else {
        throw std::string("Error: invalid lower bound set for this node selection strategy.");
    }
}

/* Compute idx with a full strong branching approach on a w.s. in function of the branching decisions bd.
 *
 * \param branchingDecisions* bd. The branching decisions made so far.
 */
int Node::getFullStrongBranchingWSIndex(BranchingDecisions* bd) {

    std::vector<bool> isCandidate(P->get_n());
    //dynamic_cast<LinearRelaxation*>(LB)->getFractionalCandidateSet(isCandidate, bd); //mof
    for (int i = 0; i < P->get_n(); i++) {
        isCandidate[i] = bd->lb[i] != bd->ub[i];
    }

    std::vector<std::vector<double>> scorel(P->get_n(), std::vector<double>(2, 0));
    std::vector<double> l(P->get_p(), 1);
    //std::vector<double> l(param->searchDir);
    bool solved;

    for (int i = 0; i < P->get_n(); i++) {
        if (isCandidate[i]) { // check free variables // bd->ub[i] != bd->lb[i]

            // branch to 0
            bd->ub[i] = 0;
            ws->adjustBounds(*bd);
            solved = ws->solve(*P, l);
            //std::cout << ws->retrieveObjectiveValue() << " vs " << score << "\n";
            if (solved) {
                if (score >= 0) scorel[i][0] = ws->retrieveObjectiveValue(); //  / score
                // max(0, ws->retrieveObjectiveValue() - EPS_INT - score); // *P, l
                else scorel[i][0] = ws->retrieveObjectiveValue(); //  score /
            }
            else {
                scorel[i][0] = 0; // 1000000
            }

            // branch to 1
            bd->ub[i] = 1;
            bd->lb[i] = 1;
            ws->adjustBounds(*bd);
            solved = ws->solve(*P, l);
            //std::cout << ws->retrieveObjectiveValue() << " vs " << score << "\n";
            if (solved) {
                if (score >= 0) scorel[i][1] = ws->retrieveObjectiveValue(); // / score
                // max(0, ws->retrieveObjectiveValue() - EPS_INT - score); // *P, l
                else scorel[i][1] = ws->retrieveObjectiveValue(); // score / 
            }
            else {
                scorel[i][1] = 0; // 1000000
            }

            // reset
            bd->lb[i] = 0;
        }
    }

    // get min score

    double s, mins = 1000000, maxs = -1000000;
    int idx = -1;

    for (int i = 0; i < P->get_n(); i++) {
        if (isCandidate[i]) { // bd->ub[i] != bd->lb[i]
            //s = (1 - MU) * min(score[i][0], score[i][1]) + MU * max(score[i][0], score[i][1]);
            s = scorel[i][0] * scorel[i][1];
            //std::cout << "score of " << i << " : " << s << "\n";
            std::cout << "x[" << i << "] : " << scorel[i][0] << " and " << scorel[i][1] << "\n"; //  << " => " << s 
            if (s > maxs) { // s < mins
                maxs = s; // mins
                idx = i;
            }
        }
    }

    SLUB sss(bd->slub);
    idx = dynamic_cast<LinearRelaxation*>(LB)->computeMostOftenFractionalIndex(sss, bd);
    std::cout << "  -> " << idx << " is chosen.\n";

    return idx;
}






bool Node::isMorePromising(Node* nd) {
    return score <= nd->getScore();
}

int Node::getLargestFractionalReducedCost(BranchingDecisions* bd, std::vector<double>& rc, std::vector<double>& sol0, std::vector<double>& sol1) {

    int idx = -1;
    int val = -1;
    double eps = 0.000001;
    std::vector<bool> candidate(P->get_n(), false);
    bool atLeastOneCandidate = false;

    // search for fractional non-fixed variables

    for (int i = 0; i < P->get_n(); i++) {
        if (bd->lb[i] != bd->ub[i]) {
            if (abs(sol0[i] - 1) <= eps || abs(sol0[i]) <= eps) {
                candidate[i] = true;
                atLeastOneCandidate = true;
            }
        }
    }

    // search for non-fixed variables if no fractional value occur

    if (!atLeastOneCandidate) {
        for (int i = 0; i < P->get_n(); i++) {
            candidate[i] = bd->lb[i] != bd->ub[i];
        }
    }

    // search for the variable with largest reduced cost

    for (int i = 0; i < P->get_n(); i++) {
        if (candidate[i]) {
            if (rc[i] > val) {
                val = rc[i];
                idx = i;
            }
        }
    }

    return idx;
}

/* Choose the variable that has a dimension the closest to the slub
*/
int Node::selectBranchingIndex1(BranchingDecisions* newbd0, BranchingDecisions* newbd1, std::vector<BranchingDecisions>& dec0, std::vector<BranchingDecisions>& dec1) {

    int index = -1;
    double closest = 10000000;

    //std::cout << "\n region: ";
    for (int i = 0; i < P->get_n(); i++) {
        for (int k = 0; k < P->get_p(); k++) {
            dec0[i].slub[k] = min(dec0[i].slub[k], P->getUbObj(k));
            dec1[i].slub[k] = min(dec1[i].slub[k], P->getUbObj(k));
            //std::cout << newbd->slub[k] << " ";
        }

    }
    for (int k = 0; k < P->get_p(); k++) {
        newbd0->slub[k] = min(newbd0->slub[k], P->getUbObj(k));
        //std::cout << newbd->slub[k] << " ";
    }
    //std::cout << "\n";

    // search in the x_i = 0 sub-pb
    for (int i = 0; i < P->get_n(); i++) { // look at each variable
        if (newbd0->lb[i] != newbd0->ub[i]) {
            for (int k = 0; k < P->get_p(); k++) { // compute the distance for each objective
                //std::cout << " val = " << dec0[i].slub[k] - dec0[i].yIprobing[k] << " vs " << closest << "\n";
                dec0[i].score = dec0[i].slub[k] - dec0[i].yIprobing[k];
                if (dec0[i].score < closest) {
                    closest = dec0[i].slub[k] - dec0[i].yIprobing[k];
                    index = i;
                    //std::cout << "okok\n";
                }
            }
        }
    }

    // search in the x_i = 1 sub-pb
    for (int i = 0; i < P->get_n(); i++) { // look at each variable
        if (newbd1->lb[i] != newbd1->ub[i]) {
            for (int k = 0; k < P->get_p(); k++) { // compute the distance for each objective
                dec1[i].score = dec1[i].slub[k] - dec1[i].yIprobing[k];
                if (dec1[i].score < closest) {
                    closest = dec1[i].slub[k] - dec1[i].yIprobing[k];
                    index = i;
                    //std::cout << "okok\n";
                }
            }
        }
    }
    //std::cout << "closest is : " << closest << "\n";

    //std::cout << "index = " << index << "\n";
    
    if (index == -1)
        throw std::string("Not valid index splitting\n");

    return index;
}

/* Branch on the variable with the smallest volume of the cube defined by (yI,slub)
 */
//int Node::selectBranchingIndex2(BranchingDecisions* newbd0, BranchingDecisions* newbd1) {
//
//    int index = -1;
//    int smallestVolume = INT_MAX;
//    int M = 100000000;
//    int minForEach;
//    int maxForAll;
//
//    //std::cout << "\n region: ";
//    for (int k = 0; k < P->get_p(); k++) {
//        newbd->slub[k] = min(newbd->slub[k], M);
//        //std::cout << newbd->slub[k] << " ";
//    }
//    //std::cout << "\n";
//
//    maxForAll = 200000000;
//    for (int i = 0; i < P->get_n(); i++) {
//        if (newbd->lb[i] != newbd->ub[i]) {
//            minForEach = 0;
//            for (int v = 0; v <= 1; v++) {
//                //std::cout << " x[" << i << "]_" << v << " -> ";
//                for (int k = 0; k < P->get_p(); k++) {
//                    //std::cout << newbd->yIprobing[i][v][k] << " ";
//                    if (newbd->slub[k] - newbd->yIprobing[k] > minForEach) {
//                        minForEach = newbd->slub[k] - newbd->yIprobing[k];
//                    }
//                }
//                if (maxForAll > minForEach) {
//                    maxForAll = minForEach;
//                    index = i;
//                }
//            }
//
//        }
//    }
//
//    if (index == -1)
//        throw std::string("Not valid index splitting\n");
//
//    return index;
//}

/* Branch on the variable that has the worst weighted sum
 */
int Node::selectBranchingMinScore(BranchingDecisions* newBd0, BranchingDecisions* newBd1, std::vector<BranchingDecisions>& dec0, std::vector<BranchingDecisions>& dec1) {

    int index(-1);
    double minScore = 1000000000000000;

    for (int i = 0; i < P->get_n(); i++) {
        if (newBd0->lb[i] != newBd0->ub[i]) {
            //std::cout << " Score is " << dec0[i].score << "\n";
            if (dec0[i].score < minScore) {
                minScore = dec0[i].score;
                index = i;
            }
            if (dec1[i].score < minScore) {
                minScore = dec0[i].score;
                index = i;
            }
        }
    }

    return index;
}

void Node::setUpScoreWS() {
    prob->setUpScoreWS(P);
}


double Node::getScore() {
    return score; // branchingDecision.
}










void Node::mergeSlubs(std::list<SLUB*>* S) {

    std::vector<int> intersectionPoint(P->get_p(), 0);
    std::list<SLUB*>::iterator s1 = S->begin(), s2, sTempo;
    bool aMergingOperationOccured = true;
    bool jobFinished = false;
    while (aMergingOperationOccured) { // it can be improved, e.g. by taking care of a list of merged slubs and checking those only
        s1 = S->begin();
        aMergingOperationOccured = false;
        while (s1 != S->end()) { // !jobFinished
            s2 = s1;
            s2++;
            while (s2 != S->end()) {
                for (int k = 0; k < P->get_p(); k++) {
                    intersectionPoint[k] = min((*s1)->get_coordinate(k), (*s2)->get_coordinate(k));
                }
                if (LB->dominates(intersectionPoint)) {
                    (*s1)->merge(**s2);
                    delete* s2;
                    S->erase(s2);
                    s2 = s1;
                    aMergingOperationOccured = true;
                }
                s2++;
            }
            s1++;
        }
    }

}











void Node::printCteFixed(std::vector<bool>& free, int j) {

    P->printConstraints(j);
    std::cout << "----------------------------\n";
    for (int i = 0; i < P->get_n(); i++) {
        std::cout << free[i] << " ";
    }
    std::cout << "\n";

}

bool Node::operator>(const Node& nd) const {
    return score > nd.score;
}
