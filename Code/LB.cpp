#include "LB.h"

/* ==========================================================
        Constructors
 ========================================================= */

/*! \brief Default constructor of a lower bound set.
 *
 * \param lp MathematicalModel*. A pointer to the problem this lower bound set refers to.
 */
LowerBoundSet::LowerBoundSet(MathematicalModel* lp, BranchingDecisions* branch) : lp(lp), status(UNSOLVED), branchDec(branch), iteration(0) {}

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
void LowerBoundSet::applyDominanceTest(UpperBoundSet& U, Parameters* param, std::list<int>& ndLub) {
	throw std::string("You tried to do a dominance test with a non-defined lower bound set");
}

/*! \brief A virtual function that checks whether this lower bound set the super local upper bound a point y.
 *
 * \param U UpperBoundSet. The upper bound set used to do the dominance test.
 * \return true if y is dominated by the lower bound set.
 */
bool LowerBoundSet::dominates(std::vector<int>& y) {
	throw std::string("You tried to do a dominance test with a non-defined lower bound set");

	return false;
}

/*! \brief Set the number of the iteration.
 *
 * \param it int. The number of the iteration.
 * \return true if y is dominated by the lower bound set.
 */
void LowerBoundSet::setIteration(int it) {
	iteration = it;
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
	for (int k = 0; k < boundingBox.size(); k++) {
		delete boundingBox[k];
	}
	//std::cout << "\nLB set destroyed\n";
}

LinearRelaxation::LinearRelaxation(MathematicalModel* lp, WeightedSumModel* ws, FeasibilityCheckModel* feas, DualBensonModel* db, FurthestFeasiblePointModel* bvp, BranchingDecisions* branch, Statistics* S) : LowerBoundSet(lp,branch), weightedSum(ws), feasibilityCheck(feas), dualBenson(db), furthestFeasiblePoint(bvp), antiIdealPoint(lp->get_p()), interiorPoint(lp->get_p()), boundingBox(0), warmstarted(false), facets(0), extrPoints(0), nbIterations(0), call(0), firstGeneration(true), S(S), checkPointDestroyed(false) { //stat()
    //initialize(*lp);
	for (int k = 0; k < lp->get_p(); k++) {
		for (int i = 0; i < lp->get_n(); i++) {
			//std::cout << "C[" << k << "," << i << "] = " << lp->get_objective(k, i) << "\n";
			if (lp->get_objective(k, i) >= 0) {
				antiIdealPoint[k] += lp->get_objective(k, i);
			}
			else {
				interiorPoint[k] += lp->get_objective(k, i);
			}
		}
		antiIdealPoint[k] += 2;
		//antiIdealPoint[k] *= 100000;
		interiorPoint[k] -= 1; // arbitrary, to be modified ? + do we actually need interior pts ?
	}
}

/*! \brief Default constructor of a linear relaxation, with models built internally.
 *
 * \param lp MathematicalModel*. A pointer to the problem this lower bound set refers to.
 */
LinearRelaxation::LinearRelaxation(MathematicalModel* lp, BranchingDecisions* branch, Statistics* S) : LowerBoundSet(lp,branch), antiIdealPoint(lp->get_p()), interiorPoint(lp->get_p()), boundingBox(0), warmstarted(false), facets(0), extrPoints(0), nbIterations(0), call(0), firstGeneration(true), S(S), checkPointDestroyed(false) { // stat()

	weightedSum = new WeightedSumModel();
	weightedSum->build(*lp);
	feasibilityCheck = new FeasibilityCheckModel();
	feasibilityCheck->build(*lp);
	dualBenson = new DualBensonModel();
	dualBenson->build(*lp);
	furthestFeasiblePoint = new FurthestFeasiblePointModel();
	furthestFeasiblePoint->build(*lp);

	//initialize(*lp);
	for (int k = 0; k < lp->get_p(); k++) {
		for (int i = 0; i < lp->get_n(); i++) {
			//std::cout << "C[" << k << "," << i << "] = " << lp->get_objective(k, i) << "\n";
			if (lp->get_objective(k, i) >= 0) {
				antiIdealPoint[k] += lp->get_objective(k, i);
			}
			else {
				interiorPoint[k] += lp->get_objective(k, i);
			}
		}
		antiIdealPoint[k] += 2;
		//antiIdealPoint[k] *= 100000;
		interiorPoint[k] -= 1;
		//interiorPoint[k] = antiIdealPoint[k] - 1; // arbitrary, to be modified ? + do we actually need interior pts ?
	}
}

/*! \brief The copy constructor for a linear relaxation
 *
 * \param LPrelax LinearRelaxation. The LinearRelaxation this object is a copy of.
 */
LinearRelaxation::LinearRelaxation(LinearRelaxation* LPrelax, BranchingDecisions* branch) : LowerBoundSet(LPrelax->lp,branch), weightedSum(LPrelax->weightedSum), feasibilityCheck(LPrelax->feasibilityCheck), dualBenson(LPrelax->dualBenson), furthestFeasiblePoint(LPrelax->furthestFeasiblePoint), antiIdealPoint(LPrelax->antiIdealPoint), interiorPoint(LPrelax->interiorPoint), warmstarted(true), nbIterations(0), call(LPrelax->call), firstGeneration(LPrelax->firstGeneration), S(LPrelax->S) { // CHANGE WARMSTARTED TO TRUE LATER stat()

	// Bouding box
	Hyperplane* H;
	if(LPrelax->boundingBox.size() == 1) {
		boundingBox = std::vector<Hyperplane*>(1);
		boundingBox[0] = new Hyperplane(*LPrelax->boundingBox[0]);
		LPrelax->boundingBox[0]->setCopy(boundingBox[0]);
	}
	else if (LPrelax->boundingBox.size() == LPrelax->lp->get_p()) {
		boundingBox = std::vector<Hyperplane*>(LPrelax->lp->get_p());
		for (int k = 0; k < boundingBox.size(); k++) {
			H = new Hyperplane(*LPrelax->boundingBox[k]);
			LPrelax->boundingBox[k]->setCopy(H); // remember the new address in the old object
			boundingBox[k] = H; // allocate a new address for the copy
		}
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
		P->becomesNonCheckpoint();
	}

	// Link adjecent extreme points and hyperplanes with new addresses
	std::list<Point*>::iterator adjacentVertex;
	std::list<Point*>* adjacencyList;
	std::list<Hyperplane*>* activeHyperplanes;
	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); vertex++) {
		adjacencyList = (*vertex)->get_adjList();
		for (adjacentVertex = adjacencyList->begin(); adjacentVertex != adjacencyList->end(); adjacentVertex++) {
			//if (*adjacentVertex != NULL) {
				*adjacentVertex = (*adjacentVertex)->get_copy();
			//}
		}
		
		activeHyperplanes = (*vertex)->get_activeHyperplanes();
		for (f = activeHyperplanes->begin(); f != activeHyperplanes->end(); f++) {
			*f = (*f)->get_copy();
		}
	}

	// Link Hyperplanes and defining points with new addresses
	//std::cout << "\n And a new generation comes in...\n";
	for (f = facets.begin(); f != facets.end(); f++) { // LPrelax->
		//std::cout << "    -> new facet : ";
		//(*f)->print();
		adjacencyList = (*f)->get_defPts();
		for (adjacentVertex = adjacencyList->begin(); adjacentVertex != adjacencyList->end(); adjacentVertex++) {
			//(*adjacentVertex)->print();
			//std::cout << "\n";
			*adjacentVertex = (*adjacentVertex)->get_copy();
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

	std::vector<double> refPt(lp->get_p());
	std::vector<double> y(lp->get_p());
	std::vector<Hyperplane*> H(lp->get_p());
	Point* newPts;
	Point* pt;
	std::vector<Point*> pts(lp->get_p());
	Hyperplane* hpp;
	double rhs;
	bool dummy;

	// compute first face (sum of objectives, scalarization)

	std::vector<double> normalVector(lp->get_p(), 1);
	double ws = weightedSum->retrieveObjectiveValue(*lp, normalVector);
	Hyperplane* h = new Hyperplane(normalVector, ws);

	// compute the set of points

	pt = new Point(antiIdealPoint, true);
	for (int k = 0; k < lp->get_p(); k++) {
		pts[k] = new Point(antiIdealPoint);
		pts[k]->setObjVector(k, pts[k]->findMissingCoordinate(*h, k));
		pts[k]->addActiveHyperplane(h);
		h->addVertex(pts[k]);
	}

	// connect the vertex together (adj lists)

	for (int k = 0; k < lp->get_p(); k++) {
		for (int l = 0; l < lp->get_p(); l++) {
			if (k != l)
				pts[k]->addAdjacentPoint(pts[l]);
		}
		pts[k]->addAdjacentPoint(pt);
		pt->addAdjacentPoint(pts[k]);
	}

	// create the bounding box

	boundingBox = std::vector<Hyperplane*>(lp->get_p());
	for (int k = 0; k < lp->get_p(); k++) {
		for (int l = 0; l < lp->get_p(); l++) {
			normalVector[l] = 0;
		}
		normalVector[k] = 1;
		boundingBox[k] = new Hyperplane(normalVector, antiIdealPoint[k]);
	}

	// connect the rays and the bounding box

	for (int k = 0; k < lp->get_p(); k++) {
		for (int l = 0; l < lp->get_p(); l++) {
			if (k != l) {
				pts[k]->addActiveHyperplane(boundingBox[l]);
				boundingBox[l]->addVertex(pts[k]);
			}
		}
		pt->addActiveHyperplane(boundingBox[k]);
		boundingBox[k]->addVertex(pt);
	}

	// add everything to the lb set

	facets.push_back(h);
	pt->becomesNonDegenerate();
	extrPoints.push_back(pt);
	for (int k = 0; k < lp->get_p(); k++) {
		pts[k]->becomesNonDegenerate();
		extrPoints.push_back(pts[k]);
	}
}

/*! \brief This function computes the linear relaxation
 */
void LinearRelaxation::compute() {

	call++;
	try {
		if (!isFeasible()) {
			status = INFEASIBLE;
		}
		else {

			if (!warmstarted) {
				S->timeInitialization.StartTimer();
				initialize();
				S->timeInitialization.StopTimer();
			}

			bool feasible(true);
			std::vector<double> y;
			Hyperplane* H;

			// For exploring and updating list of extreme points at the same time
			std::list<Point*>::iterator currentPoint;
			std::list<Point*>::iterator checkpoint;
			std::list<Point*>::iterator cleaner;
			std::list<Point*>* debuglol;
			std::list<Hyperplane*>::iterator f;
			Timer tps;

			bool firstCheckpointReached = false;
			int it = 0;
			bool allFeasible = false;
			bool boundaryUpdated = false;
			call = 0;

			currentPoint = extrPoints.begin();
			checkpoint = currentPoint;
			(*currentPoint)->becomesCheckpoint();
			do
			{
				tps.StartTimer();
				if ((*currentPoint) == NULL || (*currentPoint)->get_nbObj() == 0) { // just a security check
					std::cout << "call to cleaner !!\n";
					cleaner = currentPoint;
					currentPoint++;
					extrPoints.erase(cleaner);
				}
				else if (shouldBeChecked(*currentPoint)) {

					//if (iteration == 1) { // 17239
					//	std::cout << "\n\n ררררררררררררררררררררררר\n point checked : ";
					//	(*currentPoint)->print();
					//	std::cout << " (at " << *currentPoint << ") ";
					//}

					if ((*currentPoint)->isNew() || (*currentPoint)->get_nbVar() == 0) {
						S->timeFeasibilityCheck.StartTimer();
						feasible = feasibilityCheck->solve(*(*currentPoint)->get_objVector());
						S->timeFeasibilityCheck.StopTimer();
						S->lpSolved++;
					}
					else {
						if ((*currentPoint)->satisfyBranchingDecisions(branchDec)) {
							feasible = true;
						}
						else { // check for alternative pre-images
							S->timeFeasibilityCheck.StartTimer();
							feasible = feasibilityCheck->solve(*(*currentPoint)->get_objVector());
							S->timeFeasibilityCheck.StopTimer();
							S->lpSolved++;
							if (feasible) {
								(*currentPoint)->becomesNew();
							}
						}
					}

					if (feasible) { // if this extreme point is a feasible objective vector, store preimage
						//if (iteration == 1) std::cout << " is feasible\n";
						if ((*currentPoint)->isNew()) { // preimage already computed if feasible when warmstarted.. except newly computed pts !!!
							(*currentPoint)->isNowFeasible(feasibilityCheck);
						}
						(*currentPoint)->becomesCplexChecked();
						(*checkpoint)->becomesNonCheckpoint();
						checkpoint = currentPoint;
						(*currentPoint)->becomesCheckpoint();
						currentPoint++;
					}
					else // else, compute the cutting hyperplane
					{
						//if (iteration == 1) std::cout << "is not feasible\n";
						
						S->timeFeasibilityCheck.StartTimer();
						std::vector<double> normalVector = feasibilityCheck->extractNormalVector();
						double rhs = feasibilityCheck->extractConstant(*lp, branchDec);
						S->timeFeasibilityCheck.StopTimer();

						H = new Hyperplane(normalVector, rhs);

						S->timeUpdatePolyhedron.StartTimer();
						updatePolyhedron2(H,*currentPoint);
						S->timeUpdatePolyhedron.StopTimer();

						// if the previous checkpoint has been destroyed, search for a new one

						if (checkPointDestroyed) { 
							checkpoint = extrPoints.begin();
							while (checkpoint != extrPoints.end() && !shouldBeChecked(*checkpoint)) { //(*checkpoint)->isCheckedByCplex
								checkpoint++;
							}
							checkPointDestroyed = false;
							(*checkpoint)->becomesCheckpoint();
						}					
						currentPoint = checkpoint;
					}
				}
				else { // if the point is weakly dominated, just look at the next one
					//(*currentPoint)->becomesCplexChecked();
					currentPoint++;
				}
				it++;
				tps.StopTimer();
				//std::cout << "timer : " << tps.CumulativeTime("sec") << "\n";

			} while (currentPoint != extrPoints.end()); // && tps.CumulativeTime("sec") <= 1000

			filterExtremePoints();
			status = SOLVED;
			//std::cout << " \n LP solved : " << S->lpSolved << "\n";
			/*if (iteration == 1) {
				print();
				std::cout << "oof\n";
			}*/
			print();
			firstGeneration = false;
			/*if (iteration == 93) {
				throw std::string("Debug\n");
			}*/

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
 * USING THE MERGING METHOD
 * \param H Hyperplane*. A pointer to the new hyperplane added to the representation.
 * \param ptsDiscardedByH Point*. A pointer to the point that is cut out of the lower bound set by the new hyperplane H.
 */
void LinearRelaxation::updatePolyhedron2(Hyperplane* H, Point* ptsDiscardedByH) {

	/*if (iteration == 0) {
		std::cout << "\n\n------------------- \n\n";
		H->print();
		std::cout << "   -> cuts out : ";
		ptsDiscardedByH->print();
		std::cout << " at " << ptsDiscardedByH << "\n\n";
	}*/

	// compute the new points

	double lambda;
	std::list<Point*> N(0), D(0), E(0), M(0); // new pts, degenerate pts, points to explore, modified feasible points
	std::list<Point*>::iterator v, w, v1, v2;
	std::list<Point*>* adjacentVertex;
	Point* u, * newPts;
	E.push_back(ptsDiscardedByH);
	ptsDiscardedByH->becomesDiscarded();
	ptsDiscardedByH->becomesVisited();

	while (E.size() != 0) {

		u = E.front();
		E.pop_front();

		adjacentVertex = u->get_adjList();
		w = adjacentVertex->begin();

		while (w != adjacentVertex->end()) {

			v = w;
			w++;

			if ((*v)->isDiscarded(H)) {
				if (!(*v)->isVisited()) {
					(*v)->becomesVisited();
					E.push_back(*v);
				}
			}

			else if (!(*v)->isDegenerate()) {

				if (!(*v)->isModified())
					M.push_back(*v);
				(*v)->becomesModified();

				lambda = H->edgeIntersection(u, *v);
				if (lambda >= 1 - EPS_LAMBDA) { // 0.001
					u->becomesDegenerate(&N, H);
					D.push_back(u);
					M.push_back(*v);
					(*v)->becomesModified();
					break; // the unfeasible point is actually degenerate, we can leave the adjacent-vertex loop since we will only generate duplicates
				}
				else if (lambda <= 0 + EPS_LAMBDA) { // 0.001
					(*v)->becomesDegenerate(&N, H);
					D.push_back(*v);
				}
				else {
					newPts = new Point(u, *v, lambda, H);
					N.push_back(newPts);
				}
			}
		}
	}

	// debug
	std::list<Hyperplane*>* hpp;
	std::list<Hyperplane*>::iterator f;
	bool flg;
	for (v1 = extrPoints.begin(); v1 != extrPoints.end(); v1++) {
		hpp = (*v1)->get_activeHyperplanes();
		for (f = hpp->begin(); f != hpp->end(); f++) {
			flg = false;
			adjacentVertex = (*f)->get_defPts();
			for (v2 = adjacentVertex->begin(); v2 != adjacentVertex->end(); v2++) {
				if (*v1 == *v2)
					flg = true;
			}
			if (flg = false)
				std::cout << " inadequation hpp & vtx middle !!\n";
		}
	}
	

	// clear the discarded points

	clearDiscardedPoints();

	// clear the non-facet-defining hyperplane
	
	clearRedundantFacets();

	// connect the new vertices and the degenerate vertices (update adjacency lists)

	updateAjacencyLists(N, D);

	// clear status before the concerned points before the next iteration

	clearStatus(N, M);

	// add the new facet and the new points & rays to the polyhedron description

	facets.push_back(H);
	for (v = N.begin(); v != N.end(); v++)
		extrPoints.push_back(*v);

	// debug
	for (v1 = extrPoints.begin(); v1 != extrPoints.end(); v1++) {
		hpp = (*v1)->get_activeHyperplanes();
		for (f = hpp->begin(); f != hpp->end(); f++) {
			flg = false;
			adjacentVertex = (*f)->get_defPts();
			for (v2 = adjacentVertex->begin(); v2 != adjacentVertex->end(); v2++) {
				if (*v1 == *v2)
					flg = true;
			}
			if (flg = false)
				std::cout << " inadequation hpp & vtx !!\n";
		}
	}

	//print();
	call++; // go to next iteration
}

/*! \brief Return true if the point should be checked.
 *
 * This function checks whether (1) the point has already been checked by Cplex and (2) if not, if it is weakly dominated by
 * another existing point. If none of these conditions are met, the point should be checked for feasibility and the function
 * returns true.
 * \param y Point*. A pointer to the point being tested.
 * \return true if the point should be checked.
 */
bool LinearRelaxation::shouldBeChecked(Point* y) {

	// if the point y has already been checked by cplex or if it is weakly dominated by another point, there is no need to check it

	if (y->isCheckedByCplex() || isWeeklyDominated(y))
		return false;
	else
		return true;
}

/*! \brief Clear the redundant facets.
 */
void LinearRelaxation::clearRedundantFacets() {

	std::list<Hyperplane*>::iterator f = facets.begin(), fCurrent;

	while (f != facets.end()) {

		// take the next facet and apply the redundancy test
		fCurrent = f;
		f++;
		(*fCurrent)->checkRedundancy();

		// if the facet is redundant, notify its defining points and discard it
		if ((*fCurrent)->isRedundant()) {
			(*fCurrent)->print();
			(*fCurrent)->notifyDeletion();
			delete* fCurrent;
			facets.erase(fCurrent);
		}
	}
}

/*! \brief Clear the dsicarded points.
 */
void LinearRelaxation::clearDiscardedPoints() {

	std::list<Point*>::iterator v1 = extrPoints.begin(), v2;

	while (v1 != extrPoints.end()) {

		v2 = v1;
		v1++;

		if ((*v2)->isDiscarded()) {

			(*v2)->notifyDeletion();
			if ((*v2)->isCheckpoint())
				checkPointDestroyed = true;
			delete* v2;
			extrPoints.erase(v2);
		}
	}
}

/*! \brief Updates the adjacency lists of the points contained in lists N and D
 *
 * Done using the double description method (vertex enumeration).
 * \param N, list of Point*. The list of the polyhedron's new vertices.
 * \param D, list of Point*. The list of the polyhedron's degenerate vertices.
 */
void LinearRelaxation::updateAjacencyLists(std::list<Point*>& N, std::list<Point*>& D) {

	std::list<Point*>::iterator v1, v2;
	bool it1;

	// testing new vertices with ...
	for (v1 = N.begin(); v1 != N.end(); v1++) {

		// ... other new vertices
		it1 = true;
		for (v2 = v1; v2 != N.end(); v2++) {
			if (it1)
				it1 = false;
			else
				(*v1)->checkAdjacency(*v2, N, D);
		}

		// ... degenerate vertices
		for (v2 = D.begin(); v2 != D.end(); v2++)
			(*v1)->checkAdjacency(*v2, N, D);
	}

	// testing degenerate vertices with other degenerate vertices
	for (v1 = D.begin(); v1 != D.end(); v1++) {
		it1 = true;
		for (v2 = v1; v2 != D.end(); v2++) {
			if (it1)
				it1 = false;
			else
				(*v1)->checkAdjacency(*v2, N, D);
		}
	}
}

/*! \brief Clear the status of all the concerned points before the next iteration.
 *
 * Note: D (degenerate points) is included in M (modified points).
 * \param N, list of Point*. The list of the polyhedron's new vertices.
 * \param M, list of Point*. The list of the polyhedron's degenerate vertices.
 */
void LinearRelaxation::clearStatus(std::list<Point*>& N, std::list<Point*>& M) {

	std::list<Point*>::iterator u;

	// clear new points (N)

	for (u = N.begin(); u != N.end(); u++) {
		(*u)->becomesNonDegenerate();
		(*u)->becomesNonVisited();
		(*u)->clearModifications();
	}

	// clear modified points (M)

	for (u = M.begin(); u != M.end(); u++) {
		(*u)->becomesNonDegenerate();
		(*u)->becomesNonVisited();
		(*u)->clearModifications();
	}
}

void LinearRelaxation::computeStatusPoints2(Hyperplane* H, Point* u) {

	std::list<Point*>* adjVtx;
	std::list<Point*>::iterator vertex, vertex2;
	double lambda0, lambda1, val, d;

	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); vertex++) {
		(*vertex)->becomesNonDegenerate();
		d = (*vertex)->distToHyperplane(H,antiIdealPoint);
		if (*vertex == u) {
			(*vertex)->becomesDiscarded();
		}
		else if (abs(d) <= 0.001) {
			(*vertex)->becomesDegenerate();
			(*vertex)->addActiveHyperplane(H);
			H->addVertex(*vertex);
			//allDegenerateVertices.push_back(*vertex);
		}
		else if (d <= 0) {
			(*vertex)->becomesDiscarded();
		}
	}
}

void LinearRelaxation::computeStatusPoints(Hyperplane* H, Point* u) {

	std::list<Point*>* adjVtx;
	std::list<Point*>::iterator vertex, vertex2;
	double lambda0, lambda1, val;

	// reset degenerancy
	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); vertex++)
		(*vertex)->becomesNonDegenerate();

	// compute degenerancy
	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); vertex++) {

		//if (!(*vertex)->isDegenerate()) {
			adjVtx = (*vertex)->get_adjList();
			lambda0 = 100000;
			lambda1 = 100000;
			for (vertex2 = adjVtx->begin(); vertex2 != adjVtx->end(); vertex2++) {

				//if (!(*vertex2)->isDegenerate()) {
					val = (*vertex)->getEdgeCrossingValue(*vertex2, H);
					// test equality to 0 we put all everything vertex
					if (abs(val) <= 0.001) {
						//std::cout << " val = " << val << "\n";
						(*vertex)->becomesDegenerate();
					} // test equality to 1, we put everything on vertex2
					else if (abs(val - 1) <= 0.001) {
						//std::cout << " val = " << val << "\n";
						(*vertex2)->becomesDegenerate();
					}
				//}
			}
		//}
	}

	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); vertex++) {
		if (!(*vertex)->isDegenerate() && (*vertex)->below2(*H))
			(*vertex)->becomesDiscarded();
	}
}

void LinearRelaxation::debugCrashError(Point* vtx) {

	std::list<Point*>* adjl, adjl2;
	std::list<Point*>::iterator v1, v2, v3, v4;
	for (v1 = extrPoints.begin(); v1 != extrPoints.end(); v1++) {
		adjl = vtx->get_adjList();
		for (v2 = adjl->begin(); v2 != adjl->end(); v2++) {
			if (*v2 == vtx) {
				std::cout << " ok ca va crash\n";
			}
		}
	}
}


/*! \brief Correct the set of facet-defining hyperplanes by removing only face-defining ones.
 */
void LinearRelaxation::correctHyperplanes() {

	std::list<Hyperplane*>::iterator f, fCurrent;
	std::list<Point*>::iterator vertex, vertex2, vSave;
	std::list<Hyperplane*>* hpp;
	std::list<Point*>* adjList;

	// delete redundants hyperplanes
	f = facets.begin();
	while (f != facets.end()) {
		(*f)->checkRedundancy();
		f++;
	}

	// go through extreme points and delete redundant hyperplanes
	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); vertex++) {
		hpp = (*vertex)->get_activeHyperplanes();
		f = hpp->begin();
		while (f != hpp->end()) {
			fCurrent = f;
			f++;
			if ((*fCurrent)->isRedundant()) {
				hpp->erase(fCurrent);
			}
		}
	}

	// destroy the redundant hyperplanes
	f = facets.begin();
	while (f != facets.end()) {
		fCurrent = f;
		f++;
		if ((*fCurrent)->isRedundant()) {
			delete* fCurrent;
			facets.erase(fCurrent);
		}
	}

	// update the adjacency lists
	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); vertex++) {
		adjList = (*vertex)->get_adjList();
		vSave = adjList->begin();
		while (vSave != adjList->end()) {
			vertex2 = vSave;
			++vSave;
			if ((*vertex)->sizeIntersectionActiveHyperplanes(**vertex2) <= lp->get_p() - 2) {
				adjList->erase(vertex2);
			}
		}
	}

}

/*! \brief Identify the actual extreme point of the lower bound set
 */
void LinearRelaxation::filterExtremePoints() {

	for (std::list<Point*>::iterator pts = extrPoints.begin(); pts != extrPoints.end(); ++pts) {

		int k = 0;
		while (!(*pts)->isOnBoundingBox() && k < lp->get_p()) {
			if (abs((*pts)->get_objVector(k) - antiIdealPoint[k]) <= 0.00000001) { // approx here for equality. 10^-8
				(*pts)->becomesOnBoundingBox();
			}
			++k;
		}
	}
}

/*! \brief Check whether a point y is weakly dominated by an existing extreme point
 *
 * \param y Point. A pointer to the point to look at.
 * \return true if there exists an extreme point that weakly dominates y.
 */
bool LinearRelaxation::isWeeklyDominated(Point* y) {

	std::list<Point*>* adjList = y->get_adjList();
	std::list<Point*>::iterator ePts = adjList->begin();

	while (!y->is_ray() && ePts != adjList->end()) { // check in adjacency list whether a point weakly dominates y
		if ((*ePts) != y && (*ePts)->dominates(y)) {
			y->becomes_ray(); // field ray is true if the point is weakly dominated
			//y->becomesOnBoundingBox();
		}
		ePts++;
	}

	return y->is_ray();
}

/*! \brief Check whether the problem solved is feasible.
 *
 * \return true is it is feasible, false otherwise.
 */
bool LinearRelaxation::isFeasible() {

	S->lpSolved++;
	S->timeFeasibilityCheck.StartTimer();
	feasibilityCheck->solve(antiIdealPoint); // b
	S->timeFeasibilityCheck.StopTimer();

	return feasibilityCheck->getStatus();
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
		if (!(*it2)->isOnBoundingBox()) { //   (*it2)->is_ray()
			++s;
		}
	}

	std::list<Point*>* adjList;
	std::list<Point*>::iterator v;
	std::cout << "\n Extreme points representation (" << s << " extreme points + " << extrPoints.size() - s << " artificial points):\n";
	for (it2 = extrPoints.begin(); it2 != extrPoints.end(); ++it2) {
		if (!(*it2)->isOnBoundingBox()) { // is_ray     || true
			std::cout << " ( ";
			for (int k = 0; k < lp->get_p() - 1; k++) {
				std::cout << (*it2)->get_objVector(k) << " , ";
			}
			std::cout << (*it2)->get_objVector(lp->get_p() - 1) << " )\n"; // ", at address " << *it2 << "\n";

			/*adjList = (*it2)->get_adjList();
			for (v = adjList->begin(); v != adjList->end(); v++) {
				std::cout << "     -> ";
				(*v)->print();
				std::cout << "\n";
			}*/
			/*if ((*it2)->get_nbVar() != 0) {
				std::cout << " ( ";
				for (int k = 0; k < lp->get_n() - 1; k++) {
					std::cout << (*it2)->get_preImage(k) << " , ";
				}
				std::cout << (*it2)->get_preImage(lp->get_n() - 1) << " )\n";
			}*/
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
			//(*vertex)->print();
			if ((*vertex)->isInteger()) {
				//std::cout << " is integer";
				//(*vertex)->setAsIntegratedInUB(); // marche pas !!!!
				if ((*vertex)->isNew()) {
					U.updateUB(**vertex);
				}
				existsIntegerPoints = true;
			}
			//std::cout << "\n";
		}
		(*vertex)->becomesOld();
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
void LinearRelaxation::applyDominanceTest(UpperBoundSet& U, Parameters* param, std::list<int>& ndLub) {

	std::list<Hyperplane*>::iterator H;
	std::list<LocalUpperBound>* NU = U.getLubs();
	std::list<LocalUpperBound>::iterator u;
	std::list<int>::iterator locInsert;
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
	else if (param->objectiveBranching == FULL_OBJECTIVE_BRANCHING || param->objectiveBranching == CONE_OBJECTIVE_BRANCHING) {

		std::list<int>::iterator nextId = ndLub.begin();
		std::list<int>::iterator jump;
		bool lubDominated = true;
		bool weCheckedAllKnownNdLubs = false;
		for (u = NU->begin(); u != NU->end(); u++) {

			lubDominated = true; // by default, the lub is dominated
			// detect if next lub is known as dominated
			if (!weCheckedAllKnownNdLubs) {
				// we go through the list until an existing lub is found.
				while (nextId != ndLub.end() && u->get_id() > * nextId) { // while we hit a lub that does not exists, we delete it
					jump = nextId;
					jump++;
					ndLub.erase(nextId);
					nextId = jump;
				}
				// the actual test
				if (nextId == ndLub.end()) { // we went trough all the lus with a known status (new ones excluded)
					weCheckedAllKnownNdLubs = true; // we want to skip this test now because we checked all the lubs
				}
				else if (u->get_id() == *nextId) { // lub is non-dominated
					nextId++;
					lubDominated = false; // the next test is skipped
					//std::cout << " jump " << u->get_id() << "\n";
				}
			}
			
			// we proceed to the dominance test if we don't know the status of u
			if (lubDominated) {
				H = facets.begin();
				while (lubDominated && H != facets.end()) {
					if (!(u->above(*H, param))) { // if lub below an hyperplane H, not dominated
						lubDominated = false;
					}
					H++;
				}
				if (lubDominated) {
					lbDominated = false;
				}
				else { // insert sorted
					locInsert = ndLub.begin();
					while (locInsert != ndLub.end() && *locInsert <= u->get_id()) {
						locInsert++;
					}
					if (locInsert == ndLub.end()) {
						ndLub.push_back(u->get_id());
					}
					else {
						ndLub.insert(locInsert,u->get_id());
					}
					/*std::cout << " lub nd : " << u->get_id() << std::endl;
					for (locInsert = ndLub.begin(); locInsert != ndLub.end(); locInsert++) {
						std::cout << " " << *locInsert;
					}
					std::cout << "\n";*/
				}
			}
		}
		if (lbDominated) {
			status = DOMINATED;
		}

	}
	else {
		throw std::string("Error: this objective branching parameter is not supported yet for dominance test.");
	}
}

/*! \brief A virtual function that checks whether this lower bound set the super local upper bound a point y.
 *
 * \param U UpperBoundSet. The upper bound set used to do the dominance test.
 * \return true if y is dominated by the lower bound set.
 */
bool LinearRelaxation::dominates(std::vector<int>& y) {

	bool domi = true;
	std::list<Hyperplane*>::iterator f = facets.begin();

	while (domi && f != facets.end()) {
		if (!(*f)->dominates(y)) {
			domi = false;
		}
		f++;
	}

	return domi;
}

/*! \brief Compute the extreme rays derived from hyperplane H and attached to vertex y.
*
* \param vtx Point*. A pointer to the vertex.
* \param H Hyperplane*. A pointer to the hyperplane potentially generating extreme rays.
*/
void LinearRelaxation::computeExtremeRays(Point* vtx, Hyperplane* H, std::vector<Point*>* allNewRays) {

	std::vector<double> y(lp->get_p());
	Point* newRay;
	bool isValid;
	std::vector<bool> alreadyExists(lp->get_p());
	std::list<Point*>* adjVtx = vtx->get_adjList();
	std::list<Point*>::iterator vertex;

	// check if a ray already exists in each dimension in case the point is degenerate

	if (vtx->isDegenerate()) {
		for (vertex = adjVtx->begin(); vertex != adjVtx->end(); vertex++) {
			if ((*vertex)->is_ray()) {
				for (int k = 0; k < lp->get_p(); k++) {
					if ((*vertex)->get_objVector(k) == antiIdealPoint[k]) {
						alreadyExists[k] = true;
					}
				}
			}
		}
	}
	else {
		for (int k = 0; k < lp->get_p(); k++) {
			alreadyExists[k] = false;
		}
	}
	
	// compute the rays
	for (int k = 0; k < lp->get_p(); k++) {

		if (H->get_normalVector(k) == 0) { // unbounded dimension => potential extreme ray
			if (!alreadyExists[k]) {
				// compute the coordinates of the new ray and create the ray
				for (int l = 0; l < lp->get_p(); l++) {
					if (k == l) { // the unbounded dimension being considered
						y[l] = antiIdealPoint[l];
					}
					else {
						y[l] = vtx->get_objVector(l);
					}
				}

				// check the validity of the ray
				newRay = new Point(y, true);
				//newRay->setObjVector(k, newRay->findMissingCoordinate(*H, k));
				//newRay->print();
				//std::cout << " is created (ray)\n";
				isValid = newRay->updateActiveHyperplanes(*vtx, H, k);

				// connect the ray if it is valid
				if (isValid) {
					vtx->receive_ray();
					vtx->addAdjacentPoint(newRay);
					newRay->addAdjacentPoint(vtx);
					allNewRays->push_back(newRay);
				}
				else {
					delete newRay;
				}
			}
		}
	}

}

/*! \brief Computes the index of the most often fractional variable among the extreme points.
 *
 * \param slub SLUB. The slub that defines the part of the objective space to search in.
 * \return the index of the most often fractional variable, as an int.
 */
int LinearRelaxation::computeMostOftenFractionalIndex(SLUB& slub) {

	std::list<Point*>::iterator pt;
	std::vector<int> counter(lp->get_n(), 0);
	std::vector<double> avgValue(lp->get_n(), 0);
	double val;
	int nbIncludedPoints = 0;
	bool allInteger = true;
	bool empty = true;
	//std::cout << "\n------- splitting -------\n";
	for (pt = extrPoints.begin(); pt != extrPoints.end(); pt++) {
		//(*pt)->print();
		//std::cout << "\n";
		if (!(*pt)->isOnBoundingBox() && slub.dominated(*pt) && (*pt)->get_nbVar() != 0) { // last condition happen if timer threshold is reached during LB computation
			//(*pt)->printPreImage();
			nbIncludedPoints++;
			empty = false;
			for (int i = 0; i < lp->get_n(); i++) {
				val = (*pt)->get_preImage(i);
				if (val - trunc(val + 0.00000000001) >= 0.00000000001 ) { // for numerical instabilities
					counter[i]++;
					allInteger = false;
				}
				avgValue[i] += val;
			}
		}
	}

	/*std::cout << "counter:";
	for (int i = 0; i < lp->get_n(); i++) {
		std::cout << " " << counter[i];
	}
	std::cout << "\navg val:";
	for (int i = 0; i < lp->get_n(); i++) {
		std::cout << " " << avgValue[i] / nbIncludedPoints;
	}
	std::cout << "\n";*/

	int index = -1;
	if (!allInteger) { // if there exists at least one fractional variable in one of the extreme points
		//std::cout << "here\n";
		int max = 0;
		for (int i = 0; i < lp->get_n(); i++) {
			if (counter[i] > max) {
				max = counter[i];
				index = i;
			}
		}
	}
	else if (nbIncludedPoints >= 2) { // if there exists at least two extreme point (integer only)  // !empty
		double minDiff = 1;
		for (int i = 0; i < lp->get_n(); i++) {
			avgValue[i] = avgValue[i] / nbIncludedPoints;
			avgValue[i] = avgValue[i] - trunc(avgValue[i]);
			if (abs(avgValue[i] - 0.5) < minDiff) {
				index = i;
				minDiff = abs(avgValue[i] - 0.5);
			}
		}
	}
	else { // if there is no extreme point at all, or a unique integer extreme point
		//std::cout << "empty sub-pb !\n";
		int i = 0;
		while (i < lp->get_n()) { // index == -1 && 
			if (branchDec->ub[i] != branchDec->lb[i]) {
				index = i;
			}
			i++;
		}
	}

	//if (lp->isBinary() && (branchDec->ub[index] != 1 || branchDec->lb[index] != 0)) {

	//	std::cout << "\n At iteration: " << iteration << std::endl;
	//	print();
	//	std::cout << "\n OB region: ";
	//	slub.print();
	//	std::cout << "\n nb included pts: " << nbIncludedPoints << ":\n";
	//	for (pt = extrPoints.begin(); pt != extrPoints.end(); pt++) {
	//		if (!(*pt)->isOnBoundingBox() && slub.dominated(*pt)) {
	//			(*pt)->print();
	//			//(*pt)->printPreImage();
	//			std::cout << "\n";
	//		}
	//	}
	//	std::cout << "  -> index splitted: " << index << std::endl;

	//	std::cout << "counter:";
	//	for (int i = 0; i < lp->get_n(); i++) {
	//		std::cout << " " << counter[i];
	//	}
	//	std::cout << "\navg val:";
	//	for (int i = 0; i < lp->get_n(); i++) {
	//		std::cout << " " << avgValue[i] / nbIncludedPoints;
	//	}
	//	std::cout << "\n";

		/*for (int i = 0; i < lp->get_n(); i++) {
			std::cout << "----\n i = " << i << "u-l = " << branchDec->lb[i] << "-" << branchDec->ub[i] << "\n count = " << counter[i] << "\n val = " << avgValue[i];
		}*/
	/*	throw std::string("\n\nDebug: pb found in research of splitting index for a binary problem.");
	}*/

	if (index == -1) {
		//print();
		std::cout << "\n -> at iteration " << iteration << "\n";
		throw std::string("Error: invalid splitting index in ComputeMostOftenFractionalIndex");
	}
	
	return index;
}

/*! \brief Computes the median value of variable $x_i$ among the extreme points of the linear relaxation.
 *
 * \param slub SLUB. The slub that defines the part of the objective space to search in.
 * \param i int. Index of the splitting variable.
 * \return the index of the most often fractional variable, as an int.
 */
int LinearRelaxation::computeMedianSplittingValue(SLUB& slub, int i) {

	std::vector<double> valDec(0); // decimal values
	std::vector<double> valInt(0); // integer values

	std::list<Point*>::iterator pt;
	double val;
	for (pt = extrPoints.begin(); pt != extrPoints.end(); pt++) {
		if (!(*pt)->isOnBoundingBox() && slub.dominated(*pt) && (*pt)->get_nbVar() != 0) {
			val = (*pt)->get_preImage(i);
			if (val - trunc(val + 0.00000000001) >= 0.00000000001) { // for numerical instabilities -> if integer
				valDec.push_back(val);
			}
			else {
				valInt.push_back(val);
			}
		}
	}

	int splittingValue = 0;
	if (valInt.size() == 0 && valDec.size() == 0) { // we have no point in the cone
		splittingValue = branchDec->lb[i];
	}
	else if (valDec.size() == 0) { // no decimal values
		std::sort(valInt.begin(),valInt.end());
		splittingValue = valInt[static_cast<int>(floor(valInt.size()/2))];
	}
	else {
		std::sort(valDec.begin(), valDec.end());
		splittingValue = valDec[static_cast<int>(floor(valDec.size() / 2))];
	}
	if (splittingValue == lp->getUb(i)) { // see if that does not creates an indesired behaviour, e.g. infinite loop in branching
		splittingValue -= 1;
	}

	return splittingValue;
}

/* Generate the 2^p points of the box defined by antiIdealPoint and yI.
 *
 * Done in a recursive way.
 * \param yI. The ideal point, used for the definition of the box.
 * \param int index. Used for the recursion.
 */
void LinearRelaxation::generateBox(std::vector<double> yI, std::vector<double> coord, int index, std::vector<Hyperplane*> H) {

	if (index == lp->get_n()) {

		Point* newPt = new Point(coord);
		for (int k = 0; k < lp->get_p(); k++) {
			if (coord[k] == antiIdealPoint[k]) {
				newPt->addActiveHyperplane(boundingBox[k]);
				boundingBox[k]->addVertex(newPt);
			}
			else {
				newPt->addActiveHyperplane(H[k]);
				H[k]->addVertex(newPt);
			}
		}
		extrPoints.push_back(newPt);
	}
	else {

		coord[index] = yI[index];
		generateBox(yI, coord, index + 1, H);

		coord[index] = antiIdealPoint[index];
		generateBox(yI, coord, index + 1, H);
	}
}

/*! \brief Write into stat the statistics about this LP-relax
 *
 * This function update the lpStat field of stat.
 * \param stat Statistics. The data structure that store statistics for the whole Branch and Bound
 */
//void LinearRelaxation::getStatistics(Statistics* statBB) {
//
//	statBB->lpSolved += stat.lpSolved;
//	statBB->feasibilityCheckSolved += stat.feasibilityCheckSolved;
//	statBB->dualBensonSolved += stat.dualBensonSolved;
//	statBB->furthestFeasbilePointSolved += stat.furthestFeasbilePointSolved;
//
//	statBB->timeDualBenson += stat.timeDualBenson.CumulativeTime("sec");
//	statBB->timeFeasibilityCheck += stat.timeFeasibilityCheck.CumulativeTime("sec");
//	statBB->timeInitialization += stat.timeInitialization.CumulativeTime("sec");
//	statBB->timeFurthestFeasiblePoint += stat.timeFurthestFeasiblePoint.CumulativeTime("sec");
//	statBB->timeUpdatePolyhedron += stat.timeUpdatePolyhedron.CumulativeTime("sec");
//
//	//statBB->profiler += stat.profiler.CumulativeTime("sec");
//	//std::cout << statBB->profiler << "\n";
//}

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