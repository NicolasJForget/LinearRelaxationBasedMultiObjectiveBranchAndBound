#include "Tree.h"

TreeManager::TreeManager() : P(NULL), threshold(-2), L(0), nbLayers(0), ruleLayer(0), currentNodeRule(BREADTH_FIRST) {
}

/* Constructor of the Tree, initialized with the root node.
 *
 * \param Node* nd. A pointer to the root node.
 * \param Parameters* P. Parameters of the algorithm, that in particular tells how the tree should be explored.
 */
TreeManager::TreeManager(Node* nd, Parameters* param) : P(param), L(0), heapL(), nbLayers(1), ruleLayer(0), currentNodeRule(BREADTH_FIRST) {

    // set the threshold parameter
    if (P->nodeSelection == BREADTH_FIRST) threshold = -1;
    else if (P->nodeSelection == DEPTH_FIRST) threshold = 2;
    else threshold = P->threshold;

    // create a new layer and initialize it with the root node.
    L.push_back(std::list<Node*>(0));
    L.front().push_back(nd);
    heapL.push(nd);
    if (P->nodeSelection == BREADTH_FIRST) ruleLayer.push_back(BREADTH_FIRST);
    else if (P->nodeSelection == DEPTH_FIRST) ruleLayer.push_back(DEPTH_FIRST);
    else if (P->nodeSelection == HYBRID) ruleLayer.push_back(BREADTH_FIRST);
    else if (P->nodeSelection == MOST_FRACTIONAL) ruleLayer.push_back(MOST_FRACTIONAL);
    else if (P->nodeSelection == BEST_BOUND_WS) ruleLayer.push_back(BEST_BOUND_WS);
    else if (P->nodeSelection == BEST_BOUND_MAXMIN_GAP) ruleLayer.push_back(BEST_BOUND_MAXMIN_GAP);

    //if (P->nodeSelection == BEST_BOUND_WS) make_heap(L.back().begin(), L.back().end(), &nodeComparator::operator()); // min heap
    //if (P->nodeSelection == BEST_BOUND_WS) ; // min heap

}

/* Extract the next node to explore. Set currentNode to the extracted node.
 *
 * \return a pointer to the next node explored.
 */
Node* TreeManager::extractNode() {

    Node* nd = NULL;

    if (ruleLayer.back() == DEPTH_FIRST) {
        nd = L.back().back();
        L.back().pop_back();
    }
    else if (ruleLayer.back() == BREADTH_FIRST) {
        nd = L.back().front();
        L.back().pop_front();
    }
    else if (ruleLayer.back() == MOST_FRACTIONAL) {
        std::list<Node*>::iterator ndd, nd2;
        double minProp = 2, nodeProp;
        for (ndd = L.back().begin(); ndd != L.back().end(); ndd++) {
            nodeProp = (*ndd)->getScore(); // PercentageIntegralityLB
            if (nodeProp <= minProp) {
                minProp = nodeProp;
                nd2 = ndd;
            }
        }
        nd = *nd2;
        L.back().erase(nd2);
    }
    else if (ruleLayer.back() == BEST_BOUND_WS) { //  || ruleLayer.back() == BEST_BOUND_MAXMIN_GAP

        nd = heapL.top();
        heapL.pop();

        //std::cout << "\n score = " << nd->getScore() << "\n";
        //std::cout << "      -> size : " << heapL.size() + 1 << "";
    }
    else if (ruleLayer.back() == BEST_BOUND_MAXMIN_GAP) {

        std::list<Node*>::iterator n1, n2;
        double maxScore = -1000000, currentScore;

        // recompute scores
        //std::cout << "\n updating node costs: \n";
        for (n1 = L.back().begin(); n1 != L.back().end(); n1++) {
            //std::cout << (*n1)->getScore() << " --> ";
            (*n1)->computeMaxMinGapScore();
            //std::cout << (*n1)->getScore() << "\n";
        }

        // get node with minimal score
        for (n1 = L.back().begin(); n1 != L.back().end(); n1++) {
            currentScore = (*n1)->getScore();
            //std::cout << currentScore << " is tested... \n";
            if (currentScore >= maxScore) {
                maxScore = currentScore;
                n2 = n1;
            }
        }

        nd = *n2;
        L.back().erase(n2);
        //std::cout << "\n score = " << nd->getScore() << "\n";
    }
    else {
        throw std::string("Error: undefined rule for the current layer of nodes.");
    }
    
    // compute the new rule, if relevant
    if (P->nodeSelection == HYBRID) {
        double pctInt = nd->getPercentageIntegralityLB();
        //double pctFeas = nd->getPercentageFeasibilityLB();
        if (pctInt != -1) {
            if (pctInt > threshold) {
                currentNodeRule = BREADTH_FIRST;
            }
            else {
                currentNodeRule = DEPTH_FIRST;
            }
        }
    }
    else {
        currentNodeRule = P->nodeSelection;
    }

    return nd;
}

/* Add a layer to the tree, i.e. create a new sub-tree with a new exploration rule.
 */
void TreeManager::addLayer(int rule) {
    L.push_back(std::list<Node*>(0));
    ruleLayer.push_back(rule);
    nbLayers++;
    //std::cout << nbLayers << " layers\n";
    //std::cout << "switch to " << rule << "\n";
}

/* Add a new node to the tree.
 *
 * \param Node* nd. A pointer to the new node added.
 */
void TreeManager::pushNode(Node* nd) {

    // if we observe a change in the rule, we add a new layer
    if (currentNodeRule != ruleLayer.back()) {
        addLayer(currentNodeRule);
    }

    // add the node to the last layer.
    if (P->nodeSelection == BEST_BOUND_WS) { //  || P->nodeSelection == BEST_BOUND_MAXMIN_GAP
        heapL.push(nd);
    }
    else {
        L.back().push_back(nd);
    }
    //if (ruleLayer.back() == BEST_BOUND_WS) push_heap(L.back().begin(), L.back().end());
}

/* Check whether the tree is empty, i.e. there is no node to explore remaining. Remove the last empty layers.
 *
 * \return true if L is empty.
 */
bool TreeManager::isEmpty() {

    while (L.size() != 0 && L.back().size() == 0) {
        L.pop_back();
        nbLayers--;
        //std::cout << nbLayers << " layers\n";
    }

    return (L.size() == 0) || ((P->nodeSelection == BEST_BOUND_WS || P->nodeSelection == BEST_BOUND_MAXMIN_GAP) && heapL.empty());
}




bool CompNodePtrs::operator()(const Node* n1, const Node* n2) const {
    return n1->operator>(*n2);
}