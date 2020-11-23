#include "LB.h"

/* ==========================================================
        Constructors
 ========================================================= */

/*! \brief Default constructor of a lower bound set.
 *
 * \param lp MathematicalModel*. A pointer to the problem this lower bound set refers to.
 */
LowerBoundSet::LowerBoundSet(MathematicalModel* lp, BranchingDecisions* branch) : lp(lp), status(UNSOLVED), branchDec(branch) {}

/* ==========================================================
        Regular Methods
 ========================================================= */

/*! \brief A virtual function, to call the correct computation.
 *
 * This function will be implemented for each lower bound set in the subclasses. If the program enters in this compute(),
 * it throws a message.
 */
void LowerBoundSet::compute() {
	throw std::string("You tried to compute a non-defined lower bound set");
    //std::cout << "\n\n    !!!! YOU ARE NOT SUPPOSED TO GET HERE !!!!\n   -> you entered the compute() function from a lower bound set object, instead of an inherited one\n\n";
}

/*! \brief A virtual function that prints the lower bound set.
 *
 * This function will be implemented for each lower bound set in the subclasses. If the program enters in this print(),
 * it throws a message.
 */
void LowerBoundSet::print() {
	//std::cout << "\n Nothing to print, you are not in a specific lower bound set !!";
	throw std::string("You tried to print a non-defined lower bound set");
}

/*! \brief A virtual function that updates the upper bound set with integer solutions from the lower bound set.
 *
 * \param U UpperBoundSet. The upper bound set updated.
 */
void LowerBoundSet::gatherIntegerSolutions(UpperBoundSet& U) {
	throw std::string("You tried to gather integer solutions in a non-defined lower bound set");
	//std::cout << "\n\n !!!! You attempt to gather integer solutions in a non-defined lower bound set !!!!\n\n";
}

/*! \brief A virtual function that adjust the bounds of the variable given the bounds in the node nd.
 *
 * \param branchingDec BranchingDecisions. The data structure that describes the branching decisions to apply.
 */
void LowerBoundSet::applyBranchingDecisions() {
	throw std::string("You tried to apply branching decisions to a non-defined lower bound set");
}

/*! \brief A virtual function that checks whether this lower bound set is dominated by the upper bound set U.
 *
 * \param U UpperBoundSet. The upper bound set used to do the dominance test.
 * \param param Parameters. Used to check whether all the local upper bounds have to be tested.
 */
void LowerBoundSet::applyDominanceTest(UpperBoundSet& U, Parameters* param) {
	throw std::string("You tried to do a dominance test with a non-defined lower bound set");
}

/* ==========================================================
		Getters
 ========================================================= */

 /*! \brief Return the current status of the lower bound set.
  *
  * \return the identifier of the status, as an int.
  */
int LowerBoundSet::getStatus() {
	return status;
}


// ===============================================================================================================================
//							LinearRelaxation
// ===============================================================================================================================

/* ==========================================================
        Constructors & Destructor
 ========================================================= */

 /* \brief Destructor of a LinearRelaxation
  *
  */
LinearRelaxation::~LinearRelaxation() {

	std::list<Point*>::iterator pts;
	std::list<Hyperplane*>::iterator H;

	for (pts = extrPoints.begin(); pts != extrPoints.end(); pts++) {
		delete *pts;
	}
	for (H = facets.begin(); H != facets.end(); H++) {
		delete* H;
	}
	for (int k = 0; k < lp->get_p(); k++) {
		delete boundingBox[k];
	}
	//std::cout << "\nLB set destroyed\n";
}

LinearRelaxation::LinearRelaxation(MathematicalModel* lp, WeightedSumModel* ws, FeasibilityCheckModel* feas, DualBensonModel* db, FurthestFeasiblePointModel* bvp, BranchingDecisions* branch) : LowerBoundSet(lp,branch), weightedSum(ws), feasibilityCheck(feas), dualBenson(db), furthestFeasiblePoint(bvp), antiIdealPoint(lp->get_p()), interiorPoint(lp->get_p()), boundingBox(lp->get_p()), warmstarted(false), facets(0), extrPoints(0), nbIterations(0) {
    //initialize(*lp);
}

/*! \brief Default constructor of a linear relaxation, with models built internally.
 *
 * \param lp MathematicalModel*. A pointer to the problem this lower bound set refers to.
 */
LinearRelaxation::LinearRelaxation(MathematicalModel* lp, BranchingDecisions* branch) : LowerBoundSet(lp,branch), antiIdealPoint(lp->get_p()), interiorPoint(lp->get_p()), boundingBox(lp->get_p()), warmstarted(false), facets(0), extrPoints(0), nbIterations(0) {

	weightedSum = new WeightedSumModel();
	weightedSum->build(*lp);
	feasibilityCheck = new FeasibilityCheckModel();
	feasibilityCheck->build(*lp);
	dualBenson = new DualBensonModel();
	dualBenson->build(*lp);
	furthestFeasiblePoint = new FurthestFeasiblePointModel();
	furthestFeasiblePoint->build(*lp);

	//initialize(*lp);
}

/*! \brief The copy constructor for a linear relaxation
 *
 * \param LPrelax LinearRelaxation. The LinearRelaxation this object is a copy of.
 */
LinearRelaxation::LinearRelaxation(LinearRelaxation* LPrelax, BranchingDecisions* branch) : LowerBoundSet(LPrelax->lp,branch), weightedSum(LPrelax->weightedSum), feasibilityCheck(LPrelax->feasibilityCheck), dualBenson(LPrelax->dualBenson), furthestFeasiblePoint(LPrelax->furthestFeasiblePoint), antiIdealPoint(LPrelax->antiIdealPoint), interiorPoint(LPrelax->interiorPoint), warmstarted(true), nbIterations(0) { // CHANGE WARMSTARTED TO TRUE LATER

	// Bouding box
	Hyperplane* H;
	boundingBox = std::vector<Hyperplane*>(LPrelax->lp->get_p());
	for (int k = 0; k < LPrelax->lp->get_p(); k++) {
		H = new Hyperplane(*LPrelax->boundingBox[k]);
		LPrelax->boundingBox[k]->setCopy(H); // remember the new address in the old object
		boundingBox[k] = H; // allocate a new address for the copy
	}

	// Copy hyperplane representation
	std::list<Hyperplane*>::iterator f;
	facets = std::list<Hyperplane*>(0);
	for (f = LPrelax->facets.begin(); f != LPrelax->facets.end(); f++) {
		H = new Hyperplane(**f);
		(*f)->setCopy(H);
		facets.push_back(H);
	}

	// Copy extreme points
	std::list<Point*>::iterator vertex;
	Point* P;
	extrPoints = std::list<Point*>(0);
	for (vertex = LPrelax->extrPoints.begin(); vertex != LPrelax->extrPoints.end(); vertex++) {
		P = new Point(**vertex);
		(*vertex)->setCopy(P);
		extrPoints.push_back(P);
	}

	// Link adjecent extreme points and hyperplanes with new addresses
	std::list<Point*>::iterator adjacentVertex;
	std::list<Point*>* adjacencyList;
	std::list<Hyperplane*>* activeHyperplanes;
	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); vertex++) {
		adjacencyList = (*vertex)->get_adjList();
		for (adjacentVertex = adjacencyList->begin(); adjacentVertex != adjacencyList->end(); adjacentVertex++) {
			*adjacentVertex = (*adjacentVertex)->get_copy();
		}
		
		activeHyperplanes = (*vertex)->get_activeHyperplanes();
		for (f = activeHyperplanes->begin(); f != activeHyperplanes->end(); f++) {
			*f = (*f)->get_copy();
		}
	}
}

/* ==========================================================
        Regular Methods
 ========================================================= */

 /*! \brief Initialise the linear relaxation by building a simplex that contains it.
  *
  * \param linprog MathematicalModel*. A pointer to the problem this lower bound set refers to.
  */
void LinearRelaxation::initialize() { //MathematicalModel& lp

    // init anti-ideal & interior point
    //antiIdealPoint = std::vector<double>(lp.get_p());
    //interiorPoint = std::vector<double>(lp.get_p());
    for (int k = 0; k < lp->get_p(); k++) {
        for (int i = 0; i < lp->get_n(); i++) {
            if (lp->get_objective(k, i) >= 0) {
                antiIdealPoint[k] += lp->get_objective(k, i);
            }
        }
        antiIdealPoint[k] += 2;
        interiorPoint[k] = antiIdealPoint[k] - 1; // arbitrary, to be modified ? + do we actually need interior pts ?
    }

	// solve z1(x) + ... + zp(z)

	std::vector<double> normalVector(lp->get_p(), 1);
	double ws = weightedSum->retrieveObjectiveValue(*lp, normalVector); // ws value -> does it defines rhs of this hyperplane ?

	// add the new hyperplan to the LB set

	Hyperplane* h = new Hyperplane(normalVector, ws);
	facets.push_back(h);

	// update the set of points
	Point* newPts = new Point(antiIdealPoint);
	extrPoints.push_back(newPts);
	for (int k = 0; k < lp->get_p(); k++) {
		extrPoints.push_back(new Point(antiIdealPoint));
		extrPoints.back()->setObjVector(k, extrPoints.back()->findMissingCoordinate(*h, k));
		extrPoints.back()->addActiveHyperplane(h);
		h->addVertex();
	}

	// create adjacency lists for extreme points
	std::list<Point*>::iterator it1;
	int k = 0;
	std::list<Point*>::iterator it2;
	int l;
	for (it1 = extrPoints.begin(); it1 != extrPoints.end(); ++it1) {
		l = 0;
		for (it2 = extrPoints.begin(); it2 != extrPoints.end(); ++it2) {
			if (k != l) {
				(*it1)->addAdjacentPoint((*it2)->get_adress());
			}
			l++;
		}
		k++;
	}

	// add p artificial dominated plans (search cone)
	std::vector<double> nv(lp->get_p());
	boundingBox = std::vector<Hyperplane*>(lp->get_p());
	for (int k = 1; k <= lp->get_p(); k++) {
		for (int l = 0; l < lp->get_p(); l++) {
			nv[l] = 0;
		}
		nv[k - 1] = 1;
		boundingBox[k - 1] = new Hyperplane(nv, antiIdealPoint[k - 1]);
	}

	int obj = 0;
	for (it1 = extrPoints.begin(); it1 != extrPoints.end(); ++it1) {
		for (int k = 1; k <= lp->get_p(); k++) {
			if (obj != k) {
				(*it1)->addActiveHyperplane(boundingBox[k - 1]);
			}
		}
		obj++;
	}
}

/*! \brief This function computes the linear relaxation
 */
void LinearRelaxation::compute() {

	try {
		if (!isFeasible()) {
			status = INFEASIBLE;
		}
		else {

			if (!warmstarted) {
				initialize();
			}

			bool feasible(true);
			std::vector<double> y;
			Hyperplane* H;

			// For exploring and updating list of extreme points at the same time
			std::list<Point*>::iterator checkpoint;
			std::list<Point*>::iterator currentPoint = extrPoints.begin();
			bool firstCheckpointReached = false;

			do
			{
				nbIterations++;

				if (!warmstarted) {
					//std::cout << "solving feasibility check...\n";
					feasible = feasibilityCheck->solve(*(*currentPoint)->get_objVector());
				}
				else {
					// to do (case where LB is warmstarted) -> use point.feasible flag instead ???
				}

				if (feasible) { // if this extreme point is a feasible objective vector, store preimage
					if (!warmstarted) { // preimage already computed if feasible when warmstarted
						(*currentPoint)->isNowFeasible(feasibilityCheck);
					}
					checkpoint = currentPoint;
					firstCheckpointReached = true;
					currentPoint++;
				}
				else // else, compute the cutting hyperplane
				{

					// first, get a point on the boundary
					//std::cout << "solving further feasible point...\n";
					y = furthestFeasiblePoint->extractPoint(*(*currentPoint)->get_objVector(), interiorPoint);

					// afterwards, get the hyperplane of the facet y is located on
					//std::cout << "solving dual benson...\n";
					dualBenson->solve(y);
					std::vector<double> normalVector = dualBenson->extractNormalVector();
					double rhs = dualBenson->extractConstant(*lp,branchDec);
					H = new Hyperplane(normalVector, rhs);

					// then, update the polyhedron that defines the lienar relaxation
					updatePolyhedron(H);

					// delete the non-feasible point -> some list management
					if (firstCheckpointReached) {
						currentPoint = checkpoint;
						currentPoint++;
					}
					else
					{
						currentPoint = extrPoints.begin();
					}
				}

				// check if on bounding box

				// filter out-of-polytope hyperplanes using point ids

			} while (currentPoint != extrPoints.end());
			// --> include the last point to check completely ???? seems to be done actually
			// --> skip the first point ??? (anti-ideal)
			filterExtremePoints();
			status = SOLVED;
		}
	}
	catch (IloException& ie)
	{
		std::cerr << "Error in the constructor of the ModelClass : " << ie.getMessage() << ". Terminating!\n";
		exit(1);
	}
}

/*! \brief This function update the representations of the linear relaxation by adding a new hyperplane.
 *
 * \param H Hyperplane*. A pointer to the new hyperplane added to the representation.
 */
void LinearRelaxation::updatePolyhedron(Hyperplane* H) {

	std::list<Point*>::iterator vertex;
	std::list<Point*>::iterator vertex2;
	std::list <Hyperplane*>::iterator f;
	std::list <Hyperplane*>::iterator fCurrent;

	facets.push_back(H);

	//std::cout << " \n\n ---- new it ----\n";
	// search for infeasible and degenerated points
	vertex = extrPoints.begin();
	while (vertex != extrPoints.end()) {
		//(*vertex)->print();
		//std::cout << " at address" << (*vertex) << "\n";
		(*vertex)->becomesNonDegenerate();
		if ((*vertex)->locatedOn(*H)) {
			(*vertex)->becomesDegenerate();
		}
		else if ((*vertex)->below(*H)) {
			(*vertex)->becomesDiscarded();
		}
		++vertex;
	}

	// search for new extreme points (comparison with edges, by adjacency lists)
	std::list<Point*>* adjacentVertex;
	Point* newPts;
	std::vector<Point*> allNewVertices(0);

	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); ++vertex) { // for each extreme point

		if ((*vertex)->isDiscarded()) { // if it is discarded

			adjacentVertex = (*vertex)->get_adjList(); // we look at its adjacency list

			for (vertex2 = adjacentVertex->begin(); vertex2 != adjacentVertex->end(); ++vertex2) {

				if (!(*vertex2)->isDiscarded() && !(*vertex2)->isDegenerate()) { // for each non-discarded adjacent vertex, we compute the new point

					//std::cout << "\n New pts: ";
					//newPts->print();
					//std::cout << "at intersection of " << *vertex << " and " << *vertex2;
					//(*vertex)->print();
					//(*vertex2)->print();
					newPts = new Point((*vertex2)->edgeIntersection(**vertex, *H)); // new pts as intersection of H and edge
					(*vertex2)->replaceAdjVertex((*vertex)->get_adress(), newPts); // update adj vertex for feasible one
					newPts->addAdjacentPoint((*vertex2)->get_adress()); // init adj vertex for new one
					newPts->updateActiveHyperplanes(**vertex, **vertex2, H); // add a new active hyperplanes [ICI] !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					extrPoints.push_back(newPts); //new pts added to list of extreme points
					allNewVertices.push_back(newPts); // remember the new vertex				
				}
				else if ((*vertex2)->isDegenerate()) {
					(*vertex2)->removeAdjacentPoint(*vertex);
				}
			}
		}
		else if ((*vertex)->isDegenerate()) {
			(*vertex)->addActiveHyperplane(H);
			H->addVertex(); //new point to hyperplane H
			allNewVertices.push_back(*vertex); // remember the degenerate vertex
		}
	}

	// Update the adjacency lists for the new vertices + counters hyperplanes
	int s = allNewVertices.size();
	std::list<Hyperplane*>* hpp;

	for (int i = 0; i < s; i++) {
		for (int j = i + 1; j < s; j++) {
			if (allNewVertices[i]->sizeIntersectionActiveHyperplanes(*allNewVertices[j]) == H->get_dim()) { // is this true for p >= 3 ?
				allNewVertices[i]->addAdjacentPoint(allNewVertices[j]);
				allNewVertices[j]->addAdjacentPoint(allNewVertices[i]);
			}
		}
		if (!allNewVertices[i]->isDegenerate()) {
			hpp = allNewVertices[i]->get_activeHyperplanes();
			for (f = hpp->begin(); f != hpp->end(); f++) {
				(*f)->addVertex();
			}
		}
	}

	// delete non-feasible vertices & update hyperplanes
	std::list <Hyperplane*>* listHpp;
	vertex = extrPoints.begin();

	while (vertex != extrPoints.end()) { // relies on the fact that the first extreme point is the antiIdeal and is thus always feasible. This simplifies the list management.
		vertex2 = vertex;
		++vertex;
		if ((*vertex2)->isDiscarded()) {
			listHpp = (*vertex2)->get_activeHyperplanes();
			for (f = listHpp->begin(); f != listHpp->end(); f++) {
				(*f)->removeVertex();
			}
			//(*vertex2)->~Point(); // destroy this point
			delete* vertex2;
			extrPoints.erase(vertex2);
		}
	}

	// delete redundants hyperplanes
	f = facets.begin();
	while (f != facets.end()) {
		fCurrent = f;
		++f;
		if ((*fCurrent)->isRedundant()) {
			//(*fCurrent)->~Hyperplane();
			delete* fCurrent;
			facets.erase(fCurrent);
		}
	}
}

/*! \brief Identify the actual extreme point of the lower bound set
 */
void LinearRelaxation::filterExtremePoints() {
	for (std::list<Point*>::iterator pts = extrPoints.begin(); pts != extrPoints.end(); ++pts) {
		int k = 0;
		while (!(*pts)->isOnBoundingBox() && k < lp->get_p()) {
			if ((*pts)->get_objVector(k) == antiIdealPoint[k]) {
				(*pts)->becomesOnBoundingBox();
			}
			++k;
		}
	}
}

/*! \brief Check whether the problem solved is feasible.
 *
 * \return true is it is feasible, false otherwise.
 */
bool LinearRelaxation::isFeasible() {
	return feasibilityCheck->solve(antiIdealPoint);
}

/*! \brief Prints the lower bound set.
 */
void LinearRelaxation::print() {

	std::cout << "\n============= Lower Bound Set ==============\n\n";

	std::cout << "Anti-ideal point: ( ";
	for (int k = 0; k < lp->get_p() - 1; k++) {
		std::cout << antiIdealPoint[k] << " , ";
	}
	std::cout << antiIdealPoint[lp->get_p() - 1] << " )\n";

	// Hyperplane representation
	std::cout << "\nHyperplane representation(" << facets.size() << " facets):\n";
	std::list<Hyperplane*>::iterator it;
	for (it = facets.begin(); it != facets.end(); ++it) {
		std::cout << "( ";
		for (int k = 0; k < lp->get_p() - 1; k++) {
			std::cout << (*it)->get_normalVector(k) << " , ";
		}
		std::cout << (*it)->get_normalVector(lp->get_p() - 1) << " ) " << " = " << (*it)->get_rhs() << std::endl; // << "  -> defined by " << (*it)->get_nbDefiningPts() << " points" << std::endl; // , at address " << *it << std::endl;
	}

	// Extreme points representation

	std::list<Point*>::iterator it2;
	int s = 0;
	for (it2 = extrPoints.begin(); it2 != extrPoints.end(); ++it2) {
		if (!(*it2)->isOnBoundingBox()) {
			++s;
		}
	}

	std::cout << "\n Extreme points representation (" << s << " extreme points + " << extrPoints.size() - s << " artificial points):\n";
	for (it2 = extrPoints.begin(); it2 != extrPoints.end(); ++it2) {
		if (!(*it2)->isOnBoundingBox()) {
			std::cout << " ( ";
			for (int k = 0; k < lp->get_p() - 1; k++) {
				std::cout << (*it2)->get_objVector(k) << " , ";
			}
			std::cout << (*it2)->get_objVector(lp->get_p() - 1) << " )\n"; // ", at address " << *it2 << "\n";
		}
	}
	std::cout << "\n";
}

/*! \brief Updates the upper bound set with integer solutions found in the linear relaxation
 *
 * \param U UpperBoundSet. The upper bound set updated.
 */
void LinearRelaxation::gatherIntegerSolutions(UpperBoundSet& U) {

	std::list<Point*>::iterator vertex;
	int nbPts = 0;
	bool existsIntegerPoints = false;

	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); vertex++) {
		if (!(*vertex)->isOnBoundingBox()) { // if the point is not on the bouding box and is integer     && (*vertex)->isInteger()
			nbPts++;
			if (!(*vertex)->isInUB() && (*vertex)->isInteger()) {
				//std::cout << "\n new integer extreme point: ";
				//(*vertex)->print();
				(*vertex)->setAsIntegratedInUB(); // marche pas !!!!
				U.updateUB(**vertex);
				existsIntegerPoints = true;
			}
		}
	}

	if (nbPts == 1 && existsIntegerPoints) {
		status = OPTIMAL;
	}

}

/*! \brief A virtual function that adjust the bounds of the variable given the bounds in the node nd.
 *
 * \param branchingDec BranchingDecisions. The data structure that describes the branching decisions to apply.
 */
void LinearRelaxation::applyBranchingDecisions() {
	dualBenson->adjustBounds(*branchDec);
	furthestFeasiblePoint->adjustBounds(*branchDec);
	feasibilityCheck->adjustBounds(*branchDec);
	weightedSum->adjustBounds(*branchDec);
}

/*! \brief A virtual function that checks whether this lower bound set is dominated by the upper bound set U.
 *
 * \param U UpperBoundSet. The upper bound set used to do the dominance test.
 * \param param Parameters. Used to check whether all the local upper bounds have to be tested.
 */
void LinearRelaxation::applyDominanceTest(UpperBoundSet& U, Parameters* param) {

	std::list<Hyperplane*>::iterator H;
	std::list<LocalUpperBound>* NU = U.getLubs();
	std::list<LocalUpperBound>::iterator u;
	bool lbDominated = true;
	
	if (param->objectiveBranching == NO_OBJECTIVE_BRANCHING) {
		bool lubDominated = true;

		u = NU->begin();
		while (lbDominated && u != NU->end()) {
			//(*u).print();
			H = facets.begin();
			lubDominated = true;
			while (lubDominated && H != facets.end()) {
				if (!(u->above(*H,param))) { // if lub below an hyperplane H, not dominated
					lubDominated = false;
				}
				H++;
			}
			if (lubDominated) {
				lbDominated = false;
			}
			u++;
		}

		if (lbDominated) {
			status = DOMINATED;
		}
	}
	else {
		throw std::string("Error: this objective branching parameter is not supported yet for dominance test.");
	}
}

/* ==========================================================
		Getters
 ========================================================= */

/*! \brief Return the pointer to the weighted sum model.
 *
 * \return a pointer to a WeightedSumModel.
 */
WeightedSumModel* LinearRelaxation::getWeightedSumModel() {
	return weightedSum;
}

/*! \brief Return the pointer to the feasibility check model.
 *
 * \return a pointer to a FeasiblitityCheckModel.
 */
FeasibilityCheckModel* LinearRelaxation::getFeasibilityCheckModel() {
	return feasibilityCheck;
}

/*! \brief Return the pointer to the dual benson model.
 *
 * \return a pointer to a DualBensonModel.
 */
DualBensonModel* LinearRelaxation::getDualBensonModel() {
	return dualBenson;
}

/*! \brief Return the pointer to the furthest feasible point model.
 *
 * \return a pointer to a FurthestFeasiblePointModel.
 */
FurthestFeasiblePointModel* LinearRelaxation::getFurthestFeasiblePointModel() {
	return furthestFeasiblePoint;
}