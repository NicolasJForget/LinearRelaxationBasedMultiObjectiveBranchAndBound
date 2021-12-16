#include "LB2.h"

/* ==========================================================
        Constructors
 ========================================================= */

/*! \brief Default constructor of a lower bound set.
 *
 * \param lp MathematicalModel*. A pointer to the problem this lower bound set refers to.
 */
LowerBoundSet::LowerBoundSet(MathematicalModel* lp, BranchingDecisions* branch) : lp(lp), status(UNSOLVED), branchDec(branch), iteration(0) {}

LowerBoundSet::~LowerBoundSet() { }

/* ==========================================================
        Regular Methods
 ========================================================= */

/*! \brief A virtual function, to call the correct computation.
 *
 * This function will be implemented for each lower bound set in the subclasses. If the program enters in this compute(),
 * it throws a message.
 */
void LowerBoundSet::compute(UpperBoundSet* U, std::list<LocalUpperBound*> lubDomi) {
	throw std::string("You tried to compute a non-defined lower bound set");
    //std::cout << "\n\n    !!!! YOU ARE NOT SUPPOSED TO GET HERE !!!!\n   -> you entered the compute() function from a lower bound set object, instead of an inherited one\n\n";
}

/*! \brief A virtual function, to call the correct computation.
 *
 * This function will be implemented for each lower bound set in the subclasses. If the program enters in this compute(),
 * it throws a message.
 */
void LowerBoundSet::computeFull() {
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
void LowerBoundSet::applyDominanceTest(UpperBoundSet& U, Parameters* param, std::list<int>& ndLub, std::list<LocalUpperBound*>* domiLub) {
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

	for (pts = extrPoints.begin(); pts != extrPoints.end(); pts++)
		delete *pts;
	
	for (pts = extrRays.begin(); pts != extrRays.end(); pts++)
		delete *pts;

	for (H = facets.begin(); H != facets.end(); H++)
		delete* H;

	for (int k = 0; k < boundingBox.size(); k++)
		delete boundingBox[k];
}

LinearRelaxation::LinearRelaxation(MathematicalModel* lp, WeightedSumModel* ws, FeasibilityCheckModel* feas, DualBensonModel* db, FurthestFeasiblePointModel* bvp, BranchingDecisions* branch, Statistics* S, Parameters* param) : LowerBoundSet(lp,branch), weightedSum(ws), feasibilityCheck(feas), dualBenson(db), furthestFeasiblePoint(bvp), antiIdealPoint(lp->get_p()), interiorPoint(lp->get_p()), boundingBox(0), warmstarted(false), firstGeneration(true), checkPointDestroyed(false), attemptedDeleteNewHpp(false), facets(0), extrPoints(0), extrRays(0), nbIterations(0), call(0), S(S) { //stat()
    //initialize(*lp);
	for (int k = 0; k < lp->get_p(); k++) {
		for (int i = 0; i < lp->get_n(); i++) {
			//std::cout << "C[" << k << "," << i << "] = " << lp->get_objective(k, i) << "\n";
			if (lp->get_objective(k, i) >= 0) {
				antiIdealPoint[k] += lp->get_objective(k, i) * lp->getUb(i);
			}
			else {
				interiorPoint[k] += lp->get_objective(k, i) * lp->getUb(i);
			}
		}
		antiIdealPoint[k] += 2;
		//antiIdealPoint[k] *= 100000;
		interiorPoint[k] -= 1; // arbitrary, to be modified ? + do we actually need interior pts ?
	}

	/*if (param->variableSelection == PROBING_CUTS && lp->isBinary()) {
		if (branch->lb[branch->lastSplittedIndex] == 1)
			feasibilityCheck->addSumVarCut(lp, branch->cut1[branch->lastSplittedIndex]);
		else if (branchDec->lb[branchDec->lastSplittedIndex] == 0)
			feasibilityCheck->addSumVarCut(lp, branch->cut0[branch->lastSplittedIndex]);
		else
			std::cout << "eh ?!\n";
	}*/
}

/*! \brief Default constructor of a linear relaxation, with models built internally.
 *
 * \param lp MathematicalModel*. A pointer to the problem this lower bound set refers to.
 */
LinearRelaxation::LinearRelaxation(MathematicalModel* lp, BranchingDecisions* branch, Statistics* S, Parameters* param) : LowerBoundSet(lp,branch), antiIdealPoint(lp->get_p()), interiorPoint(lp->get_p()), boundingBox(0), warmstarted(false), firstGeneration(true), checkPointDestroyed(false), attemptedDeleteNewHpp(false), facets(0), extrPoints(0), extrRays(0), nbIterations(0), call(0), S(S) { // stat()

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
				antiIdealPoint[k] += lp->get_objective(k, i) * lp->getUb(i);
			}
			else {
				interiorPoint[k] += lp->get_objective(k, i) * lp->getUb(i);
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
LinearRelaxation::LinearRelaxation(LinearRelaxation* LPrelax, BranchingDecisions* branch, Parameters* param) : LowerBoundSet(LPrelax->lp,branch), weightedSum(LPrelax->weightedSum), feasibilityCheck(LPrelax->feasibilityCheck), dualBenson(LPrelax->dualBenson), furthestFeasiblePoint(LPrelax->furthestFeasiblePoint), antiIdealPoint(LPrelax->antiIdealPoint), interiorPoint(LPrelax->interiorPoint), warmstarted(true), firstGeneration(LPrelax->firstGeneration), checkPointDestroyed(LPrelax->checkPointDestroyed), attemptedDeleteNewHpp(LPrelax->attemptedDeleteNewHpp), nbIterations(0), call(LPrelax->call), S(LPrelax->S) { // CHANGE WARMSTARTED TO TRUE LATER stat()

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
		H->becomes_old();
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

	// Copy extreme rays
	extrRays = std::list<Point*>(0);
	for (vertex = LPrelax->extrRays.begin(); vertex != LPrelax->extrRays.end(); vertex++) {
		P = new Point(**vertex);
		(*vertex)->setCopy(P);
		extrRays.push_back(P);
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

	for (vertex = extrRays.begin(); vertex != extrRays.end(); vertex++) {
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

	/*if (param->variableSelection == PROBING_CUTS && lp->isBinary()) {
		if (branch->lb[branch->lastSplittedIndex] == 1)
			feasibilityCheck->addSumVarCut(lp, branch->cut1[branch->lastSplittedIndex]);
		else if (branchDec->lb[branchDec->lastSplittedIndex] == 0)
			feasibilityCheck->addSumVarCut(lp, branch->cut0[branch->lastSplittedIndex]);
		else
			std::cout << "eh ?!\n";
	}*/
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
	std::vector<double> normalVector;
	std::vector<double> y(lp->get_p());
	std::vector<Hyperplane*> H(lp->get_p());
	Hyperplane* hpp;
	double rhs;
	bool dummy;

	// compute the p first hyperplanes

	for (int k = 0; k < lp->get_p(); k++) {

		//std::cout << " new pts for hpp computation : ";
		// adjust the point used for computation
		for (int l = 0; l < lp->get_p(); l++) {
			if (k == l) {
				refPt[l] = interiorPoint[l];
			}
			else {
				refPt[l] = antiIdealPoint[l];
			}
			//std::cout << refPt[l] << " ";
		}
		//std::cout << "\n";

		// call Cplex to compute the hyperplane
		S->lpSolved++;
		dummy = feasibilityCheck->solve(refPt, iteration);
		normalVector = feasibilityCheck->extractNormalVector();
		rhs = feasibilityCheck->extractConstant(*lp, branchDec);
		y[k] = rhs;
		H[k] = new Hyperplane(normalVector, rhs);
		facets.push_back(H[k]);
	}
	//if (iteration == DEBUG_IT)
		//std::cout << "stop\n";

	// compute the p + 1 first extreme points

	Point* ePt = new Point(y);
	std::vector<Point*> eRay(lp->get_p());
	std::vector<double> ray(lp->get_p());

	for (int k = 0; k < lp->get_p(); k++) {

		// compute the coordinates of the ray
		for (int l = 0; l < lp->get_p(); l++) {
			if (k == l) {
				ray[l] = antiIdealPoint[l]; // y[l] + 10 // 10 is an arbitrary value
			}
			else {
				ray[l] = y[l];
			}
		}
		eRay[k] = new Point(ray, k);

		// compute the active hyperplanes of ray k
		for (int l = 0; l < lp->get_p(); l++) {
			if (k != l) {
				eRay[k]->addActiveHyperplane(H[l]);
				H[l]->addVertex(eRay[k]);
			}
		}

		// update the active hyperplanes of ePt
		ePt->addActiveHyperplane(H[k]);
		H[k]->addVertex(ePt);
	}

	// compute adjacency lists of the extreme point and rays

	for (int k = 0; k < lp->get_p(); k++) {

		// adjacency with the unique extreme point
		ePt->connect(eRay[k]);
		ePt->receive_ray(k);
	}

	// add the extreme points to the list extrPoints

	ePt->becomesNonDegenerate();
	extrPoints.push_back(ePt);
	for (int k = 0; k < lp->get_p(); k++) {
		eRay[k]->becomesNonDegenerate();
		extrPoints.push_back(eRay[k]);
	}

	S->nbNewFacetsCurrent += lp->get_p();
}

/*! \brief This function computes the linear relaxation
 */
void LinearRelaxation::compute(UpperBoundSet* U, std::list<LocalUpperBound*> lubDomi) {

	call++;
	try {
		if (!isFeasible()) {
			status = INFEASIBLE;
		}
		else {

			bool triggerRecomputation = attemptedDeleteNewHpp;

			if (!warmstarted || attemptedDeleteNewHpp) { // || branchDec->depth % CORRECTION_WARMSTART == 0
				if (warmstarted) {
					//std::cout << "LP relax recomputed at iteration " << iteration << "\n";
					clear();
				}
				S->timeInitialization.StartTimer();
				initialize();
				S->timeInitialization.StopTimer();
			}

			bool feasible(true);
			std::vector<double> y;
			Hyperplane* H;
			Timer cpuD = Timer();

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
			bool hyperplaneDominated = false;

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
						feasible = feasibilityCheck->solve(*(*currentPoint)->get_objVector(), iteration);
						S->timeFeasibilityCheck.StopTimer();
						/*std::cout << " Cplex called on : ";
						(*currentPoint)->print();
						std::cout << "\n";*/
						S->lpSolved++;
					}
					else {
						if ((*currentPoint)->satisfyBranchingDecisions(branchDec)) {
							feasible = true;
							S->nbFeasVtxCurrent++;
						}
						else { // check for alternative pre-images
							S->timeFeasibilityCheck.StartTimer();
							feasible = feasibilityCheck->solve(*(*currentPoint)->get_objVector(), iteration);
							S->timeFeasibilityCheck.StopTimer();
							/*std::cout << " Cplex called on : ";
							(*currentPoint)->print();
							std::cout << "\n";*/
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
						S->nbNewFacetsCurrent++;
						bool isNewWnd = false;
						for (int k = 0; k != lp->get_p(); k++) {
							if (H->get_normalVector(k) == 0)
								isNewWnd = true;
						}
						if (!isNewWnd)
							S->newFacetsNoRayCurrent++;

						cpuD.StartTimer();
						S->timeUpdatePolyhedron.StartTimer();
						updatePolyhedron3(H,*currentPoint);
						S->timeUpdatePolyhedron.StopTimer();
						cpuD.StopTimer();

						if (ENABLE_EARLY_LB_DOMI) {
							S->timeDominanceTest.StartTimer();
							hyperplaneDominated = U->testWeightedSumValue(branchDec, *H->getNormalVector(), H->get_rhs());
							S->timeDominanceTest.StopTimer();

							if (hyperplaneDominated) {
								std::cout << " we stop early LB computation.\n";
							}
						}

						// if the previous checkpoint has been destroyed, search for a new one

						if (checkPointDestroyed) { 
							checkpoint = extrPoints.begin();
							while (checkpoint != extrPoints.end() && !shouldBeChecked(*checkpoint)) { //(*checkpoint)->isCheckedByCplex
								checkpoint++;
							}
							checkPointDestroyed = false;
							if (!(checkpoint == extrPoints.end())) {
								(*checkpoint)->becomesCheckpoint();
							}
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

			} while (currentPoint != extrPoints.end() && tps.CumulativeTime("sec") <= TIME_OUT_LB && !hyperplaneDominated && (!attemptedDeleteNewHpp || triggerRecomputation)); // && tps.CumulativeTime("sec") <= 1000


			if (attemptedDeleteNewHpp && !triggerRecomputation) {
				compute(U,lubDomi);
			}
			else {
				filterExtremePoints();
				status = SOLVED;
				//std::cout << " \n LP solved : " << S->lpSolved << "\n";
				/*if (iteration == 1) {
					print();
					std::cout << "oof\n";
				}*/
				//exportIteration();
				if (iteration == DEBUG_IT) { // DEBUG_IT
					//print();
					exportProblem();
					std::cout << "stop";
				}
				firstGeneration = false;
				/*if (iteration == 93) {
					throw std::string("Debug\n");
				}*/

				S->cpuCurrent = cpuD.CumulativeTime("sec");
				for (std::list<Point*>::iterator vtx = extrPoints.begin(); vtx != extrPoints.end(); vtx++) {
					if (!(*vtx)->is_ray())
						S->vtxOnlyCurrent++;
				}
				bool isWnd = false;
				for (std::list<Hyperplane*>::iterator f = facets.begin(); f != facets.end(); f++) {
					isWnd = false;
					for (int i = 0; i != lp->get_p(); i++) {
						if ((*f)->get_normalVector(i) == 0)
							isWnd = true;
					}
					if (!isWnd)
						S->facetsNoRayCurrent++;
				}
			}
			
		}
		//std::cout << "extr pts: " << extrPoints.size() << "\n";

		/*for (std::list<Hyperplane*>::iterator h = facets.begin(); h != facets.end(); h++) {
			std::cout << " facet rhs: " << (*h)->get_rhs() << "\n";
		}*/
		//showPercentageIntegrality();
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
void LinearRelaxation::updatePolyhedron3(Hyperplane* H, Point* ptsDiscardedByH) {

	int laput;
	if (iteration == DEBUG_IT) { // DEBUG_IT
		//print();
		exportIteration();
		//print();
		std::cout << "\n\n------------------- \n\n";
		H->print();
		//std::cout << " enter a nb to continue... ";
		//std::cin >> laput;
		std::cout << "   -> cuts out : ";
		ptsDiscardedByH->print();
		std::cout << " at " << ptsDiscardedByH << "\n\n";

		// print detail
		/*for (std::list<Point*>::iterator vtx = extrPoints.begin(); vtx != extrPoints.end(); vtx++) {
			std::cout << "\n new vtx : ";
			(*vtx)->print();
			std::cout << "\n Adj list :\n";
			std::list<Point*>* adj = (*vtx)->get_adjList();
			for (std::list<Point*>::iterator vtx2 = adj->begin(); vtx2 != adj->end(); vtx2++) {
				std::cout << "    -> ";
				(*vtx2)->print();
				std::cout << "\n";
			}
			std::cout << " Active hpp :\n";
			std::list<Hyperplane*>* actHpp = (*vtx)->get_activeHyperplanes();
			for (std::list<Hyperplane*>::iterator f = actHpp->begin(); f != actHpp->end(); f++) {
				std::cout << "    -> ";
				(*f)->print();
			}
			std::cout << "\n";
		}

		std::cout << "\n";
		for (std::list<Hyperplane*>::iterator f = facets.begin(); f != facets.end(); f++) {
			std::cout << " new facet : ";
			(*f)->print();
			std::list<Point*>* defpts = (*f)->get_defPts();
			std::cout << " Adj list :\n";
			for (std::list<Point*>::iterator y = defpts->begin(); y != defpts->end(); y++) {
				std::cout << "   -> ";
				(*y)->print();
				std::cout << "\n";
			}
			std::cout << "\n";
		}*/
	}

	// compute the new points

	double lambda;
	std::vector<double> sol;
	std::list<Point*> N(0), D(0), E(0), M(0), R(0); // new pts, degenerate pts, points to explore, modified feasible points, extreme rays
	std::list<Point*>::iterator v, w, v1, v2, v3;
	std::list<Point*>* adjacentVertex, * adjvtx1, * adjvtx2;
	Point* u, * newPts;
	E.push_back(ptsDiscardedByH);
	ptsDiscardedByH->becomesDiscarded();
	ptsDiscardedByH->becomesVisited();
	int itv = 0;

	while (E.size() != 0) {

		u = E.front();
		E.pop_front();

		if (iteration == DEBUG_IT) {
			std::cout << " \n\n -------- \n";
			u->print();
			std::cout << " is selected\n";
		}

		adjacentVertex = u->get_adjList();
		w = adjacentVertex->begin();

		while (w != adjacentVertex->end()) {

			v = w;
			w++;

			if (iteration == DEBUG_IT) { // call == PRINT_DEBUG
				std::cout << "   -> ";
				(*v)->print();
			}


			if ((*v)->isDegenerate(H)) { //rechg
				if (!(*v)->isVisited()) {
					D.push_back(*v);
					M.push_back(*v);
					(*v)->becomesModified();
					(*v)->notifyRays(H, M);
					(*v)->becomesVisited();
				}
				if (iteration == DEBUG_IT) // call == PRINT_DEBUG
					std::cout << " is degenerate";
			}

			else if ((*v)->isDiscarded(H)) {
				if (iteration == DEBUG_IT)
					std::cout << " is discarded ";
				if (!(*v)->isVisited() && !(*v)->is_ray()) {
					(*v)->becomesVisited();
					E.push_back(*v);
					if (iteration == DEBUG_IT) { // call == PRINT_DEBUG
						std::cout << " and added to E";
					}
				}
			}

			else if (!(*v)->isDegenerate()) {

				if (!(*v)->isModified())
					M.push_back(*v);
				(*v)->becomesModified();

				sol = H->edgeIntersection2(u, *v);
				newPts = new Point(u, *v, sol, H);
				N.push_back(newPts);
				if (iteration == DEBUG_IT) { // call == PRINT_DEBUG
					std::cout << "  => ";
					newPts->print();
					std::cout << " is created\n           -> from infeasible : ";
					u->print();
					std::cout << " " << u << "\n           -> to feasible : ";
					(*v)->print();
					std::cout << " " << *v << "\n";
				}
			}
			if (iteration == DEBUG_IT) { // call == PRINT_DEBUG
				std::cout << "\n";
			}
		}
	}

	// generate extreme rays

	for (int k = 0; k < lp->get_p(); k++) {
		if (H->get_normalVector(k) == 0) { // unbounded dimension    // H->get_normalVector(k) == 0

			// generate rays for the new points
			for (v = N.begin(); v != N.end(); v++) {
				if (!(*v)->has_ray(k)) {
					if (iteration == DEBUG_IT) {
						std::cout << " ... generation on ";
						(*v)->print();
						std::cout << " in direction " << k << "\n";
					}
					generateRay(*v, k, R);
				}
			}

			// generate rays for the degenerate vertices
			for (v = D.begin(); v != D.end(); v++) {
				if (!(*v)->has_ray(k)) {
					if (iteration == DEBUG_IT) {
						std::cout << " ... generation on ";
						(*v)->print();
						std::cout << " in direction " << k << "\n";
					}
					generateRay(*v, k, R);
				}
			}
		}
	}

	// clear the discarded points

	clearDiscardedPoints();

	// clear the non-facet-defining hyperplane
	
	clearRedundantFacets();

	// connect the new vertices and the degenerate vertices (update adjacency lists

	updateAjacencyLists(N, D); // extrPoints D

	// clear status of the concerned points before the next iteration

	clearStatus(N, R, M);

	// add the new facet and the new points & rays to the polyhedron description

	facets.push_back(H);
	for (v = N.begin(); v != N.end(); v++)
		extrPoints.push_back(*v);
	for (v = R.begin(); v != R.end(); v++)
		extrPoints.push_back(*v);

	if (iteration == DEBUG_IT && call == 1) {
		//print();
		exportIteration();
		std::cout << "out";
	}
	call++; // go to next iteration
}


/*! \brief Clear the status of all the concerned points before the next iteration.
 *
 * Note: D (degenerate points) is included in M (modified points).
 * \param N, list of Point*. The list of the polyhedron's new vertices.
 * \param M, list of Point*. The list of the polyhedron's degenerate vertices.
 * \param R, list of Point*. The list of the polyhedron's new extreme rays
 */
void LinearRelaxation::correctRedundancy(std::list<Point*>& N, std::list<Point*>& D) {

	std::list<Point*>::iterator v, v1, v2, v3;

	// merge with new points

	v = N.begin();
	while (v != N.end()) {
		v2 = v;
		v2++;
		while (v2 != N.end()) {
			v3 = v2;
			v2++;
			if ((*v)->isVeryCloseTo(*v3)) {
				(*v)->merge2(*v3);
				(*v3)->notifyDeletion();
				delete* v3;
				N.erase(v3);
			}
		}
		v++;
	}

	// merge with degenerate points

	for (v = D.begin(); v != D.end(); v++) {
		v1 = N.begin();
		while (v1 != N.end()) {
			v2 = v1;
			v1++;
			if ((*v)->isVeryCloseTo(*v2)) {
				std::cout << " bruh \n";
				(*v)->merge2(*v2);
				(*v2)->notifyDeletion();
				delete* v2;
				N.erase(v2);
			}
		}
	}
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

	if (y->isCheckedByCplex() || y->is_ray()) // isWeeklyDominated(y)
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

			if (iteration == DEBUG_IT) { //
				std::cout << " Hpp deleted : ";
				(*fCurrent)->print();
			}

			if ((*fCurrent)->is_new()) {
				attemptedDeleteNewHpp = true;
			}
			else {
				(*fCurrent)->notifyDeletion();
				delete* fCurrent;
				facets.erase(fCurrent);
			}
		}
	}
}

/*! \brief Clear the dsicarded points.
 */
void LinearRelaxation::clearDiscardedPoints() {

	std::list<Point*>::iterator v1 = extrPoints.begin(), v2;
	std::list<Point*>* adj;

	while (v1 != extrPoints.end()) {

		v2 = v1;
		v1++;

		if ((*v2)->isDiscarded()) {
			if (iteration == DEBUG_IT) {
				(*v2)->print();
				std::cout << " is discarded\n";
			}
			(*v2)->notifyDeletion();
			if ((*v2)->isCheckpoint())
				checkPointDestroyed = true;
			delete* v2;
			extrPoints.erase(v2);
		}
		/*else if ((*v2)->is_ray() && (*(*v2)->get_adjList()->begin())->isDegenerate()) {
			(*v2)->notifyDeletion();
			if ((*v2)->isCheckpoint())
				checkPointDestroyed = true;
			delete* v2;
			extrPoints.erase(v2);
		}*/
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
	bool it1, debug = false;

	// testing new vertices with ...
	for (v1 = N.begin(); v1 != N.end(); v1++) {

		/*if (iteration == 4143 && abs((*v1)->get_objVector(0) + 203.7) <= 0.1) {
			std::cout << "\n";
			(*v1)->print();
			std::cout << " is tested...\n";
		}*/
		// ... other new vertices
		it1 = true;
		for (v2 = v1; v2 != N.end(); v2++) {
			if (it1)
				it1 = false;
			else {
				/*if (iteration == 4143 && abs((*v1)->get_objVector(0) + 203.7) <= 0.1) {
					std::cout << "   -> with ";
					(*v2)->print();
					std::cout << "\n";
				}*/
				(*v1)->checkAdjacency(*v2, N, D, iteration); // D
			}
		}

		// ... degenerate vertices
		for (v2 = D.begin(); v2 != D.end(); v2++) {
			/*if (iteration == 4143 && abs((*v1)->get_objVector(0) + 203.7) <= 0.1) {
				std::cout << "   -> with ";
				(*v2)->print();
				std::cout << "\n";
			}*/
			(*v1)->checkAdjacency(*v2, N, D, iteration);
		} // D
	}

	// testing degenerate vertices with other degenerate vertices
	for (v1 = D.begin(); v1 != D.end(); v1++) {
		it1 = true;
		/*if (iteration == 4143 && abs((*v1)->get_objVector(0) + 203.7) <= 0.1) {
			std::cout << "\n";
			(*v1)->print();
			std::cout << " is tested...\n";
		}*/
		for (v2 = v1; v2 != D.end(); v2++) {
			if (it1)
				it1 = false;
			else {
				/*if (iteration == 4143 && abs((*v1)->get_objVector(0) + 203.7) <= 0.1) {
					std::cout << "   -> with ";
					(*v2)->print();
					std::cout << "\n";
				}*/
				(*v1)->checkAdjacency(*v2, N, D, iteration); // D
			}
		}
	}
}

/*! \brief Clear the status of all the concerned points before the next iteration.
 *
 * Note: D (degenerate points) is included in M (modified points).
 * \param N, list of Point*. The list of the polyhedron's new vertices.
 * \param M, list of Point*. The list of the polyhedron's degenerate vertices.
 */
void LinearRelaxation::clearStatus(std::list<Point*>& N, std::list<Point*>& M, std::list<Point*>& R) {

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
	
	// clear modified points (R)

	for (u = R.begin(); u != R.end(); u++) {
		(*u)->becomesNonDegenerate();
		(*u)->becomesNonVisited();
		(*u)->clearModifications();
	}
}

/*! \brief Generate an extreme ray in direction k at point w and filter redundant rays in R.
 *
 * \param w Point*. A pointer to the point attached to the ray computed.
 * \param k int. The unbounded direction of the ray.
 * \param R, list of Point*. The list of the polyhedron's new rays.
 */
void LinearRelaxation::generateRay(Point* w, int k, std::list<Point*>& R) {

	std::list<Hyperplane*>::iterator f, g;
	std::list<Hyperplane*>* hppW, * hppR;
	std::list<Point*>::iterator v, vSave;
	bool rayDeleted = false;

	std::vector<double> ref = { -24.5, -29.25, -27.25, -29.5, -34.25 };

	// create the ray
	
	Point* r = new Point(*w->get_objVector(), k);
	r->setObjVector(k, antiIdealPoint[k]); // r->get_objVector(k) + 10// 10 is an arbitrary value
	//if (call <= 33 && w->isVeryCloseTo(ref)) {
		//std::cout << "\n Ray generated at ";
		//w->print();
		//std::cout << " in direction " << k << " : \n";
	//}

	// compute the active hyperplanes of the ray

	hppW = w->get_activeHyperplanes();
	for (f = hppW->begin(); f != hppW->end(); f++) {
		if ((*f)->get_normalVector(k) == 0)
			r->addActiveHyperplane(*f);
	}

	// filter the set of extreme rays R for redundancy

	hppR = r->get_activeHyperplanes();
	
	if (hppR->size() <= lp->get_p() - 2) { // weaker but faster condition : a ray should define at least an edge of the polyhedron
		//if (call <= 1)
			//std::cout << " we entered size condition...\n";
		//r->print();
		//if (call <= 33 && w->isVeryCloseTo(ref)) {
			//std::cout << "     -> deleted by nb minimal hpp condition\n";
		//}
		delete r;
		rayDeleted = true;
	}
	else { // check for facet inclusion

		// with the new rays
		vSave = R.begin();
		while (!rayDeleted && vSave != R.end()) {
			v = vSave;
			vSave++;
			if (r->hasSameFacets(*v)) {
				rayDeleted = true;
				//if (call <= 1) { //  && w->isVeryCloseTo(ref)
					//std::cout << "     -> deleted by other point included condition\n";
				//}
				delete r;
			}
			else if ((*v)->hasSameFacets(r)) {
				(*v)->notifyDeletion();
				//if (call == 1) { //  && (*v)->isVeryCloseTo(ref)
					//(*v)->print();
					//std::cout << " is deleted by nb minimal hpp condition !!! BRUUUUUH\n";
				//}
				delete* v;
				R.erase(v);
			}
		}

		// and with the other rays
		vSave = extrPoints.begin();
		while (!rayDeleted && vSave != extrPoints.end()) {
			v = vSave;
			vSave++;
			if ((*v)->is_ray()) {
				if (r->hasSameFacets(*v)) {
					rayDeleted = true;
					//(*v)->print();
					//std::cout << " zbleh\n";
					delete r;
				}
				// causes segfault error, since it may be included in M
				/*else if ((*v)->hasSameFacets(r)) {
					if (iteration == DEBUG_IT) {
						std::cout << "    -> ";
						(*v)->print();
						std::cout << " is deleted !!\n";
					}
					(*v)->notifyDeletion();
					std::cout << " bruh \n";
					delete* v;
					R.erase(v);
				}*/
			}
		}
	}

	// connect the ray to its adjacent vertex and defining facets

	if (!rayDeleted) {
		r->connect(w);
		w->receive_ray(k);
		hppR = r->get_activeHyperplanes();
		for (f = hppR->begin(); f != hppR->end(); f++)
			(*f)->addVertex(r);
		R.push_back(r);
	}
	else {
		if (iteration == DEBUG_IT)
			std::cout << "   => r is deleted !\n";
	}
	//std::cout << "\n";
}




void LinearRelaxation::debug__resetAdjLists() {

	std::list<Point*> N(0);
	std::list<Point*>::iterator v, w;
	bool it1 = false;
	for (v = extrPoints.begin(); v != extrPoints.end(); v++) {
		(*v)->debug__clearAdjList();
	}

	for (v = extrPoints.begin(); v != extrPoints.end(); v++) {
		it1 = true;
		for (w = v; w != extrPoints.end(); w++) {
			if (it1)
				it1 = false;
			else
				(*v)->checkAdjacency(*w, N, extrPoints, iteration);
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

	bool iswnd = false;
	std::list<Point*>* adjList = y->get_adjList();
	std::list<Point*>::iterator ePts = adjList->begin();

	while (!y->isOnBoundingBox() && ePts != adjList->end()) { // check in adjacency list whether a point weakly dominates y
		if ((*ePts) != y && (*ePts)->dominates(y)) {
			//y->becomes_ray(); // field ray is true if the point is weakly dominated
			y->becomesOnBoundingBox();
			iswnd = true;
		}
		ePts++;
	}

	return iswnd;//y->is_ray();
}

/*! \brief Check whether the problem solved is feasible.
 *
 * \return true is it is feasible, false otherwise.
 */
bool LinearRelaxation::isFeasible() {

	S->lpSolved++;
	S->timeFeasibilityCheck.StartTimer();
	feasibilityCheck->solve(antiIdealPoint, iteration); // b
	S->timeFeasibilityCheck.StopTimer();

	return feasibilityCheck->getStatus();
}




void LinearRelaxation::clear() {

	for (std::list<Point*>::iterator v = extrPoints.begin(); v != extrPoints.end(); v++) {
		delete* v;
	}
	extrPoints.clear();

	for (std::list<Hyperplane*>::iterator f = facets.begin(); f != facets.end(); f++) {
		delete* f;
	}
	facets.clear();

	attemptedDeleteNewHpp = false;
}

/*! \brief Prints the lower bound set.
 */
void LinearRelaxation::print() {

	//exportProblem();

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
		if (!(*it2)->is_ray()) { //   (*it2)->is_ray()   (*it2)->isOnBoundingBox()
			++s;
		}
	}

	std::list<Point*>* adjList;
	std::list<Point*>::iterator v;
	std::cout << "\n Extreme points representation (" << s << " extreme points + " << extrPoints.size() - s << " artificial points):\n";
	for (it2 = extrPoints.begin(); it2 != extrPoints.end(); ++it2) {
		if (!(*it2)->is_ray()) { // is_ray     || true   isOnBoundingBox
			if ((*it2)->is_ray())
				std::cout << " RAY ";
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
			if ((*it2)->get_nbVar() != 0) {
				std::cout << " ( ";
				for (int k = 0; k < lp->get_n() - 1; k++) {
					std::cout << (*it2)->get_preImage(k) << " , ";
				}
				std::cout << (*it2)->get_preImage(lp->get_n() - 1) << " )\n";
			}
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
				S->avgIntegerVtxCurrent++;
				//std::cout << S->avgIntegerVtxCurrent << "\n";
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
void LinearRelaxation::applyDominanceTest(UpperBoundSet& U, Parameters* param, std::list<int>& ndLub, std::list<LocalUpperBound*>* domiLub) {

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
				S->nbDominanceTestLubHpp += 1.0;
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
	else if (param->objectiveBranching == FULL_OBJECTIVE_BRANCHING || param->objectiveBranching == CONE_OBJECTIVE_BRANCHING || param->objectiveBranching == LIMITED_OBJECTIVE_BRANCHING) {

		if (ENABLE_EARLY_LB_DOMI) domiLub->clear();

		std::list<int>::iterator nextId = ndLub.begin();
		std::list<int>::iterator jump;
		bool lubDominated = true;
		bool weCheckedAllKnownNdLubs = false;
		for (u = NU->begin(); u != NU->end(); u++) {

			//lubDominated = true; // by default, the lub is dominated
			//// detect if next lub is known as dominated
			//if (!weCheckedAllKnownNdLubs) {
			//	// we go through the list until an existing lub is found.
			//	while (nextId != ndLub.end() && u->get_id() > * nextId) { // while we hit a lub that does not exists, we delete it
			//		jump = nextId;
			//		jump++;
			//		ndLub.erase(nextId);
			//		nextId = jump;
			//	}
			//	// the actual test
			//	if (nextId == ndLub.end()) { // we went trough all the lus with a known status (new ones excluded)
			//		weCheckedAllKnownNdLubs = true; // we want to skip this test now because we checked all the lubs
			//	}
			//	else if (u->get_id() == *nextId) { // lub is non-dominated
			//		nextId++;
			//		lubDominated = false; // the next test is skipped
			//		//std::cout << " jump " << u->get_id() << "\n";
			//	}
			//}
			//
			//// we proceed to the dominance test if we don't know the status of u
			//if (lubDominated) {
			//	H = facets.begin();
			//	while (lubDominated && H != facets.end()) {
			//		if (!(u->above(*H, param))) { // if lub below an hyperplane H, not dominated
			//			lubDominated = false;
			//		}
			//		H++;
			//	}
			//	if (lubDominated) {
			//		lbDominated = false;
			//	}
			//	else { // insert sorted
			//		//nextId--;
			//		locInsert = nextId;//ndLub.begin();
			//		while (locInsert != ndLub.end() && *locInsert <= u->get_id()) {
			//			locInsert++;
			//			//std::cout << " buluuh\n";
			//		}
			//		//std::cout << u->get_id() << " placed before " << *nextId;
			//		if (locInsert == ndLub.end()) {
			//			ndLub.push_back(u->get_id());
			//		}
			//		else {
			//			ndLub.insert(locInsert,u->get_id());
			//		}
			//		/*std::cout << " lub nd : " << u->get_id() << std::endl;
			//		for (locInsert = ndLub.begin(); locInsert != ndLub.end(); locInsert++) {
			//			std::cout << " " << *locInsert;
			//		}
			//		std::cout << "\n";*/
			//	}
			//}

			H = facets.begin();
			lubDominated = true;
			while (lubDominated && H != facets.end()) {
				if (!(u->above(*H, param))) { // if lub below an hyperplane H, not dominated
					lubDominated = false;
				}
				S->nbDominanceTestLubHpp += 1.0;
				H++;
			}

			if (lubDominated) {
				lbDominated = false;
				domiLub->push_back(&(*u));
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
					vtx->receive_ray(k);
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
		(*pt)->becomesUnknownDominater();
		// OnBoundingBox
		if (!(*pt)->is_ray() && (*pt)->get_nbVar() != 0) { // last condition happen if timer threshold is reached during LB computation
			//(*pt)->printPreImage();
			if (slub.dominated(*pt)){
				//std::cout << " is considered";
				(*pt)->becomesDominater();
				nbIncludedPoints++;
				empty = false;
				for (int i = 0; i < lp->get_n(); i++) {
					val = (*pt)->get_preImage(i);
					if (val - trunc(val + 0.0000001) >= 0.0000001) { // for numerical instabilities
						counter[i]++;
						allInteger = false;
					}
					avgValue[i] += val;
				}
			}
			else {
				(*pt)->becomesNonDominater();
			}
			//std::cout << "\n";
		}
	}

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
		double minDiff = 10000;
		double refVal = 0;
		for (int i = 0; i < lp->get_n(); i++) {
			if (branchDec->ub[i] != branchDec->lb[i]) {
				refVal = double(branchDec->ub[i] - branchDec->lb[i]) / 2;
				//refVal = (branchDec->ub[i] - branchDec->lb[i]) / 2;
				avgValue[i] = avgValue[i] / nbIncludedPoints;
				avgValue[i] = avgValue[i] - trunc(avgValue[i]);
				//std::cout << i << " -> " << avgValue[i] << " vs " << refVal << " ";
				if (abs(avgValue[i] - refVal) < minDiff) {
					index = i;
					minDiff = abs(avgValue[i] - 0.5);
					//std::cout << " => minDiff = " << minDiff;
				}
				//std::cout << "\n";
			}
		}
	}
	else { // if there is no extreme point at all, or a unique integer extreme point
		int i = 0;
		while (i < lp->get_n()) { // index == -1 && 
			if (branchDec->ub[i] != branchDec->lb[i]) {
				index = i;
			}
			i++;
		}
	}

	if (index == -1) {
		std::cout << "\n -> at iteration " << iteration << "\n";
		std::cout << "stopere plz";
		throw std::string("Error: invalid splitting index in ComputeMostOftenFractionalIndex");
	}

	//std::cout << "index is: " << index << "\n\n";
	
	return index;
}

/*! \brief Computes the index of the most often fractional variable among the extreme points.
 *
 * \param slub SLUB. The slub that defines the part of the objective space to search in.
 * \return the index of the most often fractional variable, as an int.
 */
int LinearRelaxation::computeMostOftenFractionalIndex(SLUB& slub, BranchingDecisions* newbd) {

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
		(*pt)->becomesUnknownDominater();
		// OnBoundingBox
		if (!(*pt)->is_ray() && (*pt)->get_nbVar() != 0) { // last condition happen if timer threshold is reached during LB computation
			//(*pt)->printPreImage();
			if (slub.dominated(*pt)) {
				//std::cout << " is considered";
				(*pt)->becomesDominater();
				nbIncludedPoints++;
				empty = false;
				for (int i = 0; i < lp->get_n(); i++) {
					if (newbd->ub[i] != newbd->lb[i]) {
						val = (*pt)->get_preImage(i);
						if (val - trunc(val + 0.0000001) >= 0.0000001) { // for numerical instabilities
							counter[i]++;
							allInteger = false;
						}
						avgValue[i] += val;
					}
				}
			}
			else {
				(*pt)->becomesNonDominater();
			}
			//std::cout << "\n";
		}
	}

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
		for (int i = 0; i < lp->get_n(); i++) {
			//if (counter[i] == max) std::cout << "check\n";
		}
	}
	else if (nbIncludedPoints >= 2) { // if there exists at least two extreme point (integer only)  // !empty
		double minDiff = 10000;
		double refVal = 0;
		for (int i = 0; i < lp->get_n(); i++) {
			if (newbd->ub[i] != newbd->lb[i]) {
				refVal = double(branchDec->ub[i] - branchDec->lb[i]) / 2;
				//refVal = (branchDec->ub[i] - branchDec->lb[i]) / 2;
				avgValue[i] = avgValue[i] / nbIncludedPoints;
				avgValue[i] = avgValue[i] - trunc(avgValue[i]);
				//std::cout << i << " -> " << avgValue[i] << " vs " << refVal << " ";
				if (abs(avgValue[i] - refVal) < minDiff) {
					index = i;
					minDiff = abs(avgValue[i] - 0.5);
					//std::cout << " => minDiff = " << minDiff;
				}
				//std::cout << "\n";
			}
		}
	}
	else { // if there is no extreme point at all, or a unique integer extreme point
		int i = 0;
		while (i < lp->get_n()) { // index == -1 && 
			if (newbd->ub[i] != newbd->lb[i]) {
				index = i;
			}
			i++;
		}
	}

	if (index == -1) {
		std::cout << "\n -> at iteration " << iteration << "\n";
		std::cout << "stopere plz";
		throw std::string("Error: invalid splitting index in ComputeMostOftenFractionalIndex");
	}

	//std::cout << "index is: " << index << "\n\n";

	return index;
}

/*! \brief Computes the median value of variable $x_i$ among the extreme points of the linear relaxation.
 *
 * \param slub SLUB. The slub that defines the part of the objective space to search in.
 * \param i int. Index of the splitting variable.
 * \return the index of the most often fractional variable, as an int.
 */
int LinearRelaxation::computeMedianSplittingValue2(SLUB& slub, int i) {

	std::vector<double> valDec(0); // decimal values
	std::vector<double> valInt(0); // integer values
	std::vector<double> valAll(0);

	std::list<Point*>::iterator pt;
	double val;
	for (pt = extrPoints.begin(); pt != extrPoints.end(); pt++) {
		if (!(*pt)->is_ray() && slub.dominated(*pt) && (*pt)->get_nbVar() != 0) {
			val = (*pt)->get_preImage(i);
			if (val - trunc(val + 0.0000001) >= 0.0000001) { // for numerical instabilities -> if integer
				valDec.push_back(val);
			}
			else {
				valInt.push_back(val);
			}
			valAll.push_back(val);
		}
	}

	int splittingValue = 0;
	//if (valInt.size() == 0 && valDec.size() == 0) { // we have no point in the cone
	//	splittingValue = branchDec->lb[i];
	//}
	//else if (valDec.size() == 0) { // no decimal values
	//	std::sort(valInt.begin(),valInt.end());
	//	splittingValue = valInt[static_cast<int>(floor(valInt.size()/2))];
	//}
	//else {
	//	std::sort(valDec.begin(), valDec.end());
	//	splittingValue = valDec[static_cast<int>(floor(valDec.size() / 2))];
	//}

	//if (splittingValue == branchDec->ub[i]) { // see if that does not creates an indesired behaviour, e.g. infinite loop in branching // lp->getUb(i)
	//	splittingValue -= 1;
	//}

	if (valAll.size() == 0) {
		splittingValue = branchDec->lb[i];
	}
	else {
		std::sort(valAll.begin(), valAll.end());
		int middle = static_cast<int>(floor(valAll.size() / 2));
		if (valAll.size() % 2 == 0) {
			splittingValue = floor((valAll[middle - 1] + valAll[middle]) / 2);
		}
		else {
			splittingValue = floor(valAll[middle]);
		}
	}

	return splittingValue;
}

/*! \brief Computes the median value of variable $x_i$ among the extreme points of the linear relaxation.
 *
 * \param slub SLUB. The slub that defines the part of the objective space to search in.
 * \param i int. Index of the splitting variable.
 * \return the index of the most often fractional variable, as an int.
 */
int LinearRelaxation::computeMedianSplittingValue(SLUB& slub, int i, int maxUb) {

	std::vector<double> valDec(0); // decimal values
	//std::vector<double> valInt(0); // integer values

	std::list<Point*>::iterator pt;
	double val;
	for (pt = extrPoints.begin(); pt != extrPoints.end(); pt++) {
		if (!(*pt)->is_ray() && slub.dominated(*pt) && (*pt)->get_nbVar() != 0) {
			val = (*pt)->get_preImage(i);
			valDec.push_back(val);
		}
	}

	int splittingValue = 0;
	if (valDec.size() == 0) { // we have no point in the cone
		splittingValue = branchDec->lb[i];
	}
	else {
		std::sort(valDec.begin(), valDec.end());
		//splittingValue = valDec[static_cast<int>(floor(valDec.size() / 2))];
		int middle = static_cast<int>(floor(valDec.size() / 2));
		if (valDec.size() % 2 == 0) {
			splittingValue = floor((valDec[middle - 1] + valDec[middle]) / 2);
		}
		else {
			splittingValue = floor(valDec[middle]);
		}
	}

	if (splittingValue == branchDec->ub[i]) { // see if that does not creates an indesired behaviour, e.g. infinite loop in branching // lp->getUb(i)
		splittingValue -= 1;
	}

	//std::cout << "split val: " << splittingValue << "\n";

	return min(splittingValue, maxUb - 1);
}

/*! \brief Computes most often fractional value of variable $x_i$ among the extreme points of the linear relaxation.
 *
 * \param slub SLUB. The slub that defines the part of the objective space to search in.
 * \param i int. Index of the splitting variable.
 * \return the index of the most often fractional variable, as an int.
 */
int LinearRelaxation::computeMostOftenFractionalSplittingValue(SLUB& slub, int i, int maxUb) {

	int maxVal = 0;
	int avgVal = 0;
	int nbPts = 0;

	// get max val

	std::list<Point*>::iterator pt;
	double val;
	for (pt = extrPoints.begin(); pt != extrPoints.end(); pt++) {
		//(*pt)->becomesUnknownDominater();
		if (!(*pt)->is_ray() && (*pt)->get_nbVar() != 0) {
			if (slub.dominated(*pt)) {
				//(*pt)->becomesDominater();
				val = (*pt)->get_preImage(i);
				if (val >= maxVal)
					maxVal = val;
				//std::cout << " " << val;
				avgVal += val;
				nbPts++;
			}
			else {
				//(*pt)->becomesNonDominater();
			}
		}
	}
	avgVal /= nbPts;

	// get frational val

	std::vector<int> nbFracVal(floor(maxVal) + 1, 0);
	for (pt = extrPoints.begin(); pt != extrPoints.end(); pt++) {
		if (!(*pt)->isOnBoundingBox() && slub.dominated(*pt) && (*pt)->get_nbVar() != 0) {
			val = (*pt)->get_preImage(i);
			if (val - trunc(val + 0.0000001) >= 0.0000001) // if not integer
				nbFracVal[floor(val)]++; // each cell correspond to one integer value, and we look at the floor of values of extr pts to get what is the int value
		}
	}

	// search for the actual value

	int vl = 0;
	if (nbFracVal.size() == 0) {
		vl = branchDec->lb[i];
	}
	else {
		vl = std::distance(nbFracVal.begin(), std::max_element(nbFracVal.begin(), nbFracVal.end())); // we get the most occuring value
		if (nbFracVal[vl] == 0) { // no decimal value => call median rule
			//vl = computeMedianSplittingValue(slub, i);
			//vl = floor(avgVal);
			vl = computeRandomSplittingValue(slub, i, branchDec->ub[i]);
		}
	}

	return min(vl, maxUb - 1);
}

/*! \brief Computes most often fractional value of variable $x_i$ among the extreme points of the linear relaxation.
 *
 * \param slub SLUB. The slub that defines the part of the objective space to search in.
 * \param i int. Index of the splitting variable.
 * \return the index of the most often fractional variable, as an int.
 */
int LinearRelaxation::computeRandomSplittingValue(SLUB& slub, int i, int maxUb) {

	int valSplit = -1;
	int nbPts = 0;
	int maxVal = branchDec->lb[i];
	int minVal = branchDec->ub[i] - 1;

	// get min & max val

	std::list<Point*>::iterator pt;
	double val;
	for (pt = extrPoints.begin(); pt != extrPoints.end(); pt++) {
		if (!(*pt)->is_ray() && (*pt)->get_nbVar() != 0) {
			if (slub.dominated(*pt)) {
				nbPts++;
				val = (*pt)->get_preImage(i);
				if (floor(val) > maxVal)
					maxVal = floor(val);
				if (floor(val) < minVal)
					minVal = floor(val);
				//std::cout << " " << val;
			}
		}
	}

	if (nbPts == 0) {
		maxVal = branchDec->ub[i] - 1;
		minVal = branchDec->lb[i];
	}
	
	//srand(time(NULL)); // time(NULL)
	valSplit = rand() % (maxVal - minVal + 1) + minVal;

	//std::cout << " range : " << minVal << " , " << maxVal << "\n";

	if (valSplit == branchDec->ub[i])
		valSplit--;

	return min(valSplit, maxUb - 1);
}

void LinearRelaxation::exportIteration() {

	// computing values of parameters

	int nz = 0;
	int i = 0;
	std::list<Point*>::iterator v;
	std::list<Hyperplane*>::iterator f;

	for (f = facets.begin(); f != facets.end(); f++) {
		for (int k = 0; k < lp->get_p(); k++) {
			if ((*f)->get_normalVector(k) != 0)
				nz++;
		}
	}

	// writing to the file

	std::ofstream file;
	file.open("C:/Users/au643334/Documents/PhD/Code/Bensolve/home/au643334/bensolve-2.1.0/ex/iteration.vlp"); // append instead of overwrite

	// vlp carac
	file << "p vlp min " << facets.size() << " " << lp->get_p() << " " << nz << " " << lp->get_p() << " " << lp->get_p() << "\n";

	// constraints
	i = 0;
	for (f = facets.begin(); f != facets.end(); f++) {
		++i;
		for (int k = 0; k < lp->get_p(); k++) {
			if ((*f)->get_normalVector(k) != 0)
				file << "a " << i << " " << k + 1 << " " << (*f)->get_normalVector(k) << "\n";
		}
	}

	// objectives
	for (int k = 0; k < lp->get_p(); k++) {
		file << "o " << k + 1 << " " << k + 1 << " 1\n";
	}

	// rhs constraints
	i = 0;
	for (f = facets.begin(); f != facets.end(); f++) {
		++i;
		file << "i " << i << " l " << (*f)->get_rhs() << "\n";
	}

	// variables bounds
	for (int k = 0; k < lp->get_p(); k++) {
		file << "j " << k + 1 << " f\n";
	}

	file << "e";
}

void LinearRelaxation::exportProblem() {

	// computing values of parameters

	int nz = 0;
	int i = 0;
	std::list<Point*>::iterator v;
	std::list<Hyperplane*>::iterator f;

	for (int i = 0; i < lp->get_m(); i++) {
		for (int j = 0; j < lp->get_n(); j++) {
			if (lp->get_constraint(i, j) != 0) {
				nz++;
			}
		}
	}

	int ctr = 0;
	for (int i = 0; i != lp->get_n(); i++) {
		if (branchDec->lb[i] != 0) {
			nz++;
			ctr++;
		}
		if (branchDec->ub[i] != 100) {
			nz++;
			ctr++;
		}
	}

	// writing to the file

	std::ofstream file;
	file.open("C:/Users/au643334/Documents/PhD/Code/Bensolve/home/au643334/bensolve-2.1.0/ex/iteration.vlp"); // append instead of overwrite

	// vlp carac
	file << "p vlp min " << lp->get_m() + ctr << " " << lp->get_n() << " " << nz << " " << lp->get_p() << " " << lp->get_p() * lp->get_n() << "\n";

	// constraints
	for (int i = 0; i < lp->get_m(); i++) {
		for (int j = 0; j < lp->get_n(); j++) {
			if (lp->get_constraint(i, j) != 0) {
				file << "a " << i + 1 << " " << j + 1 << " " << lp->get_constraint(i, j) << "\n";
			}
		}
	}
	
	int ctesBranching = 0;
	for (int i = 0; i != lp->get_n(); i++) {
		if (branchDec->lb[i] != 0) {
			++ctesBranching;
			file << "a " << lp->get_m() + ctesBranching << " " << i + 1 << " 1\n";
		}
		if (branchDec->ub[i] != 100) {
			++ctesBranching;
			file << "a " << lp->get_m() + ctesBranching << " " << i + 1 << " 1\n";
		}
	}

	// objectives
	for (int k = 0; k < lp->get_p(); k++) {
		for (int i = 0; i < lp->get_n(); i++) {
			file << "o " << k + 1 << " " << i + 1 << " " << lp->get_objective(k, i) << "\n";
		}
	}

	// rhs constraints
	std::string signCte;
	for (int j = 0; j < lp->get_m(); j++) {
		if (lp->get_signCte(j) == 0) signCte = " l ";
		else if (lp->get_signCte(j) == 1) signCte = " u ";
		else signCte = " s ";

		file << "i " << j + 1 << signCte << lp->get_rhs(j) << "\n";
	}

	ctesBranching = 0;
	for (int i = 0; i != lp->get_n(); i++) {
		if (branchDec->lb[i] != 0) {
			++ctesBranching;
			file << "i " << lp->get_m() + ctesBranching << " l " << branchDec->lb[i] << "\n";
		}
		if (branchDec->ub[i] != 100) {
			++ctesBranching;
			file << "i " << lp->get_m() + ctesBranching << " u " << branchDec->ub[i] << "\n";
		}
	}
	

	// variables bounds
	for (int i = 0; i < lp->get_n(); i++) {
		file << "j " << i + 1 << " d 0 100\n";
	}

	file << "e";
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

/* \brief Clear all the cuts generated in the model used to comptue the linear relaxation.
 */
void LinearRelaxation::clearCuts() {
	feasibilityCheck->clearCuts();
}

/*! \brief Add the cover cut described in cc to the model.
 *
 * \param cc vector of vector of int. Represent the cover cuts. Each row is a new cut, column 0 is the rhs, and the other columns are the indices of the variables to add to the cut.
 */
void LinearRelaxation::applyCoverCuts(std::vector<std::vector<int>>& cc) {
	feasibilityCheck->applyCoverCuts(cc);
}

/*! \brief Class the indices of the variables from the most to the least fractional among the extreme points of the lower bound set.
 *
 * \param cc vector int. Indices are stored in this vector.
 */
void LinearRelaxation::computeFractionalProportion(std::vector<int>* I, SLUB& s) {

	std::vector nbFrac = std::vector<int>(I->size(), 0); // count how many times each variable is fractional
	std::list<Point*>::iterator y;
	double val;

	// compute the proportion of fractional extreme points for each variable

	for (y = extrPoints.begin(); y != extrPoints.end(); y++) {
		if (!(*y)->is_ray() && (*y)->get_nbVar() != 0 && s.dominated(*y)) { // if not a ray, has a pre-image, and is in the dominance cone of s
			for (int i = 0; i < lp->get_n(); i++) {
				val = (*y)->get_preImage(i);
				if (val - trunc(val + 0.0000001) >= 0.0000001) { // variable is fractional for y
					nbFrac[i]++;
				}
			}
		}
	}

	// sort variables from most to least often fractional

	int max, idx = -1;
	for (int i = 0; i < lp->get_n(); i++) {
		max = -1;
		for (int i2 = 0; i2 < lp->get_n(); i2++) {
			if (nbFrac[i2] > max) {
				idx = i2;
				max = nbFrac[i2];
			}
		}
		(*I)[i] = idx;
		nbFrac[idx] = -1;
	}

	// print

	/*std::cout << " I : ";
	for (int i = 0; i < lp->get_n(); i++) {
		std::cout << (*I)[i] << " ";
	}
	std::cout << "\n";*/
}

std::vector<double> LinearRelaxation::getNadirPoint() {

	std::list<Point*>::iterator y;
	std::vector<double> yN(lp->get_p(), -1000000);

	// calculation

	for (y = extrPoints.begin(); y != extrPoints.end(); y++) {
		if (!(*y)->is_ray()) {
			for (int k = 0; k < lp->get_p(); k++) {
				if (yN[k] < (*y)->get_objVector(k)) {
					yN[k] = (*y)->get_objVector(k);
				}
			}
		}
	}

	// print
	//std::cout << " LP nadir point : ";
	//for (int k = 0; k < lp->get_p(); k++) {
		//std::cout << yN[k] << " ";
	//}
	//std::cout << "\n";

	return yN;
}

/*! \brief Check for each free variable wether each extreme value occurs among extreme points of the lower bound set in the region described in bd. Values occuring are notified in h;
 *
 * \param bd BranchingDecisions*. A pointer to the branching decisions used for reference.
 * \param h vector of vector of bool. If a variable take an extreme value, its corresponding cell is set to false.
 */
void LinearRelaxation::checkFeasibleValues(BranchingDecisions* bd, std::vector<std::vector<bool>>* h) {

	std::list<Point*>::iterator y;
	SLUB s = SLUB(bd->slub);

	for (y = extrPoints.begin(); y != extrPoints.end(); y++) {
		(*y)->becomesNonRefVarFix();
		if (!(*y)->is_ray() && s.dominated(*y)) { // if y is in the region of the obj space described by s
			//std::cout << " new point tested : ";
			for (int i = 0; i < lp->get_n(); i++) {
				//std::cout << (*y)->get_preImage(i) << " ";
				if (bd->lb[i] != bd->ub[i]) { // if variable i is free
					if ((*y)->get_preImage(i) == 0) {
						(*h)[i][0] = false;
						(*y)->becomesRefVarFix();
					}
					else if ((*y)->get_preImage(i) == 1) {
						(*h)[i][1] = false;
						(*y)->becomesRefVarFix();
					}
				}
			}
			//std::cout << "\n";
		}
	}

}

/* \brief Check whether new values for fixing variables should be tested afted variable i was fixed to v.
 *
 * \param toTest, vector of vector of bool. Where the values to check for each variables are recorded.
 * \param int i. The index of the variable that was fixed.
 * \param int v. The value to which x_i was fixed.
 */
void LinearRelaxation::updateFeasibleValues(std::vector<std::vector<bool>>* toTest, BranchingDecisions* bd, int i, int v) {

	std::list<Point*>::iterator y;
	int cptFeas = 0;

	// reset ref pts
	for (y = extrPoints.begin(); y != extrPoints.end(); y++) {
		if ((*y)->isRefVarFix()) {
			if ((*y)->get_preImage(i) == v) {
				cptFeas++;
			}
			else {
				(*y)->becomesNonRefVarFix();
			}
		}
	}

	// check values of variables
	std::vector<int> cptr0(lp->get_n(), 0);
	std::vector<int> cptr1(lp->get_n(), 0);
	for (y = extrPoints.begin(); y != extrPoints.end(); y++) {
		if ((*y)->isRefVarFix()) {
			//(*y)->print();
			for (int l = 0; l < lp->get_n(); l++) {
				if ((*y)->get_preImage(l) == 0) {
					cptr0[l]++;
				}
				else if ((*y)->get_preImage(l) == 1) {
					cptr1[l]++;
				}
			}
		}
	}

	// update tests to be performed
	for (int l = 0; l < lp->get_n(); l++) {
		if (l == 9) {
			//std::cout << "after Lb set check: x9 in [" << bd->lb[l] << " , " << bd->ub[l] << "]\n";
		}
		if (bd->lb[l] != bd->ub[l]) {
			if (!(*toTest)[l][0] && cptr0[l] == 0 && bd->lb[l] != 1) { // no extr pts has value 0 for x_l and x_l is not fixed to 1
				(*toTest)[l][0] = true;
				std::cout << " /!! we now need to test x[" << l << "] = 0\n";
			}
			if (!(*toTest)[l][1] && cptr1[l] == 0 && bd->ub[l] != 0) { // no extr pts has value 1 for x_l and x_l is not fixed to 0
				(*toTest)[l][1] = true;
				std::cout << " /!! we now need to test x[" << l << "] = 1\n";
			}
		}
	}
}

/*! \brief Compute the index that violates the most often the branching decisions among the extreme points of the lower bound set.
 *
 * \param BranchingDecisions* newbd. A pointer to the branching decisions used for violation test.
 * \return the index of the variable chosen, as an int.
 */
int LinearRelaxation::computeMostViolatingIndex(BranchingDecisions* newbd) {

	int index = -1;
	std::list<Point*>::iterator y;
	std::vector<int> nbViolations(lp->get_n(), 0);
	std::vector<bool> isViolated;

	// compute the number of violations before branching
	//for (y = extrPoints.begin(); y != extrPoints.end(); y++) {
	//	if (!(*y)->is_ray()) { // we look at extreme points only, not rays
	//		for (int i = 0; i < lp->get_n(); i++) {
	//			if (newbd->lb[i] != newbd->ub[i]) { // we look at free variables only
	//				if ((*y)->get_preImage(i) != newbd->lb[i]) nbViolations[i]++;
	//				if ((*y)->get_preImage(i) != newbd->ub[i]) nbViolations[i]++;
	//			}
	//		}
	//	}
	//}
	for (y = extrPoints.begin(); y != extrPoints.end(); y++) {
		if (!(*y)->is_ray()) { // we look at extreme points only, not rays
			for (int i = 0; i < lp->get_n(); i++) {
				if (newbd->lb[i] != newbd->ub[i]) { // we look at free variables only
					if (!(*y)->satisfyBranchingDecisions(newbd) || (*y)->get_preImage(i) != newbd->lb[i]) nbViolations[i]++;
					if (!(*y)->satisfyBranchingDecisions(newbd) || (*y)->get_preImage(i) != newbd->ub[i]) nbViolations[i]++;
				}
			}
		}
	}

	// get the most violating one
	int maxVal = 0;
	//std::cout << "\n new node : ";
	for (int i = 0; i < lp->get_n(); i++) {
		//std::cout << nbViolations[i] << " ";
		if (nbViolations[i] > maxVal) {
			index = i;
			maxVal = nbViolations[i];
		}
	}

	return index;
}

/* Get the minimal difference between the rhs of an hyperplane and the w.s. of local upper bound u using the normal vector of the hyperplane
 * as the weigth vector.
 *
 * \param LocalUpperBound u. The local upper bound used for reference.
 */
double LinearRelaxation::getMinFacetGap(LocalUpperBound& u) {

	std::list<Hyperplane*>::iterator f;
	double minGap = 10000000;
	double currentGap = 10000000;
	std::vector<double> l(lp->get_p());

	for (f = facets.begin(); f != facets.end(); f++) {
		for (int k = 0; k < lp->get_p(); k++) {
			l[k] = (*f)->get_normalVector(k);
		}
		currentGap = u.getWeightedSum(l) - (*f)->get_rhs();
		if (currentGap < minGap) {
			minGap = currentGap;
		}
	}

	return minGap;
}

/* Retrieve the w.s. value of the extreme point with the smallest w.s. value, using weights l.
 *
 * \param vector of double, l. The weight vector.
 */
double LinearRelaxation::getSmallestWeightedSumValue(std::vector<double>& l) {

	double wsMin = 1000000, wsVal;
	std::list<Point*>::iterator y;

	for (y = extrPoints.begin(); y != extrPoints.end(); y++) {
		if (!(*y)->is_ray()) {
			wsVal = (*y)->getWeightedSumValue(l);
			if (wsVal < wsMin) {
				wsMin = wsVal;
			}
		}
	}

	return wsMin;
}

/* Retrieve the variables that are not fixed and that takes at least once a fractional value among the extreme points of the LB set.
 * If none takes a fractional value, returns the variables that take at least once different values.
 * If no extreme point exists in the area defined by the subproblem, the set of all free variables is returned.
 *
 * \param isCandidate, vector of bool. True if variable i is candidate, false otherwise.
 * \param BranchingDecision* bd. The branching decisions considered.
 */
void LinearRelaxation::getCandidateSet(std::vector<bool>& isCandidate, BranchingDecisions* bd) {

	int n = isCandidate.size();
	std::list<Point*>::iterator y;
	SLUB s(bd->slub);

	// search for free & decimal variables

	//std::cout << "Stage 1 is trigered. ";
	for (int i = 0; i < n; i++) {
		isCandidate[i] = bd->ub[i] != bd->lb[i]; // check if it is a free variable
		if (isCandidate[i]) { // if it is, search for fractional extreme points
			isCandidate[i] = false;
			y = extrPoints.begin();
			while (!isCandidate[i] && y != extrPoints.end()) {
				//(*y)->print(); std::cout << " vs "; s.print();
				//std::cout << " domi : " << s.dominatedFix(*y) << "\n";
				if (!(*y)->is_ray() && s.dominatedFix(*y)) {
					if ((*y)->get_preImage(i) - trunc((*y)->get_preImage(i)) > EPS_INT) {
						isCandidate[i] = true;
					}
				}
				y++;
			}

		}
	}

	// if there is no candidate, search for the integer ones taking different values by looking at their average value
	if (std::count(isCandidate.begin(), isCandidate.end(), true) == 0) {

		//std::cout << "Stage 2 is trigered... ";
		std::vector<double> avg(n, 0);
		int nbPts = 0;

		for (y = extrPoints.begin(); y != extrPoints.end(); y++) {
			if (!(*y)->is_ray() && s.dominatedFix(*y)) {
				for (int i = 0; i < n; i++) {
					avg[i] += (*y)->get_preImage(i);
				}
				nbPts++;
			}
		}
		for (int i = 0; i < n; i++) {
			avg[i] /= double(nbPts);
			isCandidate[i] = (bd->ub[i] != bd->lb[i]) && (avg[i] > EPS_INT && avg[i] - 1 < EPS_INT);
		}

		// if there is still no candidate, put all free variables as candidates
		if (std::count(isCandidate.begin(), isCandidate.end(), true) == 0) {
			
			for (int i = 0; i < n; i++) {
				isCandidate[i] = bd->ub[i] != bd->lb[i];
			}
			//std::cout << "Stage 3 is trigered !";

		}
	}
	//std::cout << "\n";
}

void LinearRelaxation::getFractionalCandidateSet(std::vector<bool>& isCandidate, BranchingDecisions* bd) {

	std::list<Point*>::iterator pt;
	std::vector<int> counter(lp->get_n(), 0);
	std::vector<double> avgValue(lp->get_n(), 0);
	SLUB s(bd->slub);
	double val;
	int nbIncludedPoints = 0;
	bool allInteger = true;
	bool empty = true;

	for (int i = 0; i < lp->get_n(); i++) isCandidate[i] = false;

	//std::cout << "\n------- splitting -------\n";
	for (pt = extrPoints.begin(); pt != extrPoints.end(); pt++) {
		//(*pt)->print();
		//std::cout << "\n";
		(*pt)->becomesUnknownDominater();
		// OnBoundingBox
		if (!(*pt)->is_ray() && (*pt)->get_nbVar() != 0) { // last condition happen if timer threshold is reached during LB computation
			//(*pt)->printPreImage();
			if (s.dominated(*pt)) {
				//std::cout << " is considered";
				(*pt)->becomesDominater();
				nbIncludedPoints++;
				empty = false;
				for (int i = 0; i < lp->get_n(); i++) {
					if (bd->ub[i] != bd->lb[i]) {
						val = (*pt)->get_preImage(i);
						if (val - trunc(val + 0.0000001) >= 0.0000001) { // for numerical instabilities
							counter[i]++;
							allInteger = false;
						}
						avgValue[i] += val;
					}
				}
			}
			else {
				(*pt)->becomesNonDominater();
			}
			//std::cout << "\n";
		}
	}

	/*std::cout << " Fractionality of variables: \n";
	for (int i = 0; i < lp->get_n(); i++) {
		std::cout << i << " -> " << counter[i] << "\n";
	}*/

	int index = -1;
	if (!allInteger) { // if there exists at least one fractional variable in one of the extreme points
		//std::cout << "here\n";
		int max = 0;
		for (int i = 0; i < lp->get_n(); i++) {
			if (counter[i] > max) {
				max = counter[i];
			}
		}
		for (int i = 0; i < lp->get_n(); i++) {
			//if (counter[i] == max) std::cout << "check\n";
			if (counter[i] == max) isCandidate[i] = true;
		}
	}
	else if (nbIncludedPoints >= 2) { // if there exists at least two extreme point (integer only)  // !empty
		double minDiff = 10000;
		double refVal = 0;
		for (int i = 0; i < lp->get_n(); i++) {
			if (bd->ub[i] != bd->lb[i]) {
				//refVal = double(bd->ub[i] - bd->lb[i]) / 2;
				////refVal = (branchDec->ub[i] - branchDec->lb[i]) / 2;
				//avgValue[i] = avgValue[i] / nbIncludedPoints;
				//avgValue[i] = avgValue[i] - trunc(avgValue[i]);
				////std::cout << i << " -> " << avgValue[i] << " vs " << refVal << " ";
				//if (abs(avgValue[i] - refVal) < minDiff) {
				//	index = i;
				//	minDiff = abs(avgValue[i] - 0.5);
				//	//std::cout << " => minDiff = " << minDiff;
				//}
				////std::cout << "\n";
				avgValue[i] = avgValue[i] / nbIncludedPoints;
				if (avgValue[i] > EPS_INT && avgValue[i] < 1 - EPS_INT) {
					isCandidate[i] = true;
				}
			}
		}
	}
	else { // if there is no extreme point at all, or a unique integer extreme point
		for (int i = 0; i < lp->get_n(); i++) {
			if (bd->ub[i] != bd->lb[i]) {
				isCandidate[i] = true;
			}
		}
	}

	/*if (index == -1) {
		std::cout << "\n -> at iteration " << iteration << "\n";
		std::cout << "stopere plz";
		throw std::string("Error: invalid splitting index in ComputeMostOftenFractionalIndex");
	}*/

	std::cout << "Candidate list : ";
	for (int i = 0; i < lp->get_n(); i++) {
		if (isCandidate[i]) std::cout << i << " ";
	}
	std::cout << "\n";

	//std::cout << "index is: " << index << "\n\n";

	//return index;
}

bool LinearRelaxation::filterNdLubs(std::list<LocalUpperBound*>& lubDomi, Hyperplane* h) {

	std::list<LocalUpperBound*>::iterator u, uNext;
	
	u = lubDomi.begin();
	uNext = u;
	while (u != lubDomi.end()) {
		u = uNext;
		uNext++;
		//if ()
	}

	return lubDomi.size() == 0;
}

bool LinearRelaxation::takeValue(BranchingDecisions* bd, int i, int v) {
	
	bool valueFound = false;
	bool feasible;
	std::list<Point*>::iterator y = extrPoints.begin();
	SLUB s = SLUB(bd->slub);
	int i2;

	while (!valueFound && y != extrPoints.end()) {

		if (!(*y)->is_ray() && (*y)->get_preImage(i) == v && s.dominated(*y)) {

			feasible = true;
			i2 = 0;
			while (feasible && i2 < lp->get_n()) {
				if (!((*y)->get_preImage(i2) >= bd->lb[i2] && (*y)->get_preImage(i2) <= bd->ub[i])) {
					feasible = false;
				}
				i2++;
			}

			if (feasible) {
				valueFound = true;
				//std::cout << "  x[" << i << "] take value " << v << " in ";
				//(*y)->printPreImage();
			}
		}
		y++;
	}

	return valueFound;
}

std::vector<double> LinearRelaxation::getAverageNormalVector(BranchingDecisions* bd) {

	std::vector<double> avgNV(lp->get_p(), 0);
	SLUB s = SLUB(bd->slub);
	std::list<Hyperplane*>::iterator h;
	std::list<Point*>::iterator y;
	std::list<Point*>* defPts;
	bool facetInSubpb;
	int nbFacetsInSubpb = 0;

	for (h = facets.begin(); h != facets.end(); h++) {
		facetInSubpb = false;
		defPts = (*h)->get_defPts();
		for (y = defPts->begin(); y != defPts->end(); y++) {
			if (s.strictlyDominated(*y)) {
				facetInSubpb = true;
			}
		}

		if (facetInSubpb) {
			for (int k = 0; k < lp->get_p(); k++) {
				avgNV[k] += (*h)->get_normalVector(k);
			}
			nbFacetsInSubpb++;
		}
	}

	if (true || nbFacetsInSubpb == 0) {
		for (int k = 0; k < lp->get_p(); k++) {
			avgNV[k] = 1;
		}
	}
	else {
		for (int k = 0; k < lp->get_p(); k++) {
			avgNV[k] /= nbFacetsInSubpb;
		}
	}

	return avgNV;
}

double LinearRelaxation::getPercentageIntegrality() {

	std::list<Point*>::iterator y;
	double pct = 0;
	int nbExtrPts = 0, nbInt = 0;

	if (extrPoints.size() != 0) {

		for (y = extrPoints.begin(); y != extrPoints.end(); y++) {
			if (!(*y)->is_ray()) {
				if ((*y)->isInteger()) { //  && (*y)->satisfyBranchingDecisions(branchDec)
					nbInt++;
				}
				nbExtrPts++;
			}
		}

		pct = double(nbInt) / double(nbExtrPts);
	}
	else {
		pct = -1; // pct to -1 if we can't decide => don't change branching rule
	}
	//std::cout << pct << "% of the LB set is integer\n";

	return pct;
}

double LinearRelaxation::getPercentageFeasibility() {

	std::list<Point*>::iterator y;
	double pct = 0;
	int nbExtrPts = 0, nbFeas = 0;

	if (extrPoints.size() != 0) {

		for (y = extrPoints.begin(); y != extrPoints.end(); y++) {
			if (!(*y)->is_ray()) {
				if ((*y)->satisfyBranchingDecisions(branchDec)) {
					nbFeas++;
				}
				nbExtrPts++;
			}
		}

		pct = double(nbFeas) / double(nbExtrPts);
	}
	else {
		pct = -1; // pct to -1 if we can't decide => don't change branching rule
	}
	//std::cout << 100 * pct << "% of the LB set is still feasible ( " << nbFeas << " / " << nbExtrPts << " )\n";

	return pct;
}

void LinearRelaxation::checkViolatingCut() {

	int s;

	for (std::list<Point*>::iterator y = extrPoints.begin(); y != extrPoints.end(); y++) {
		if (!(*y)->is_ray()) {
			for (int k = 0; k < lp->get_p(); k++) {

				s = branchDec->slub[k];
				for (int i = 0; i < lp->get_n(); i++) {
					if ((*y)->get_preImage(i) >= 1 - EPS_INT) {
						s -= int((*y)->get_preImage(i));
					}
				}
				for (int i = 0; i < lp->get_n(); i++) {
					if ((*y)->get_preImage(i) < 1 - EPS_INT && (*y)->get_preImage(i) > EPS_INT) {
						if (lp->get_objective(k, i) > s) {
							std::cout << " cover cut possible on objective " << k << "\n";
						}
					}
				}
			}
		}
	}
}

void LinearRelaxation::checkViolatingCut(BranchingDecisions* bd) {

	SLUB sl(bd->slub);
	int s;
	std::vector<int> cutIdx(0), newCut(0);
	std::vector<bool> covered(lp->get_n(), false);
	std::vector<std::vector<int>> cc(0);

	sl.print();
	for (std::list<Point*>::iterator y = extrPoints.begin(); y != extrPoints.end(); y++) {
		if (!(*y)->is_ray() && sl.dominated(*y)) {

			covered = std::vector(lp->get_n(), false);
			for (int k = 0; k < lp->get_p(); k++) {

				if (lp->get_objDir(k) == -1) { // minimization problems only for now

					cutIdx.clear();
					// set initial rhs & add variables that take value 1 to the cut
					s = min(bd->slub[k], 1000000) - 1;
					for (int i = 0; i < lp->get_n(); i++) {
						if ((*y)->get_preImage(i) >= 1 - EPS_INT && lp->get_objective(k, i) > 0) {
							s -= lp->get_objective(k, i); // int((*y)->get_preImage(i));
							cutIdx.push_back(i);
							covered[i] = true;
						}
					}

					// search for fractional variables which result in solution exceeding an OB constraint if set to 1
					for (int i = 0; i < lp->get_n(); i++) {
						if (!covered[i] && (*y)->get_preImage(i) < 1 - EPS_INT && (*y)->get_preImage(i) > EPS_INT) { // if var fractional
							if (lp->get_objective(k, i) > s) {
								(*y)->print();
								std::cout << " => cover cut possible on objective " << k << " (cost of " << lp->get_objective(k, i) << " )" << "\n";
								
								// generate a new cover cut
								newCut = std::vector<int>(cutIdx);
								newCut.push_back(i);
								newCut.push_back(cutIdx.size()); // rhs of cover cut at the end of the vector
								cc.push_back(std::vector<int>(newCut));
								covered[i] = true;
							}
						}
					}
				}	
			}
		}
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




// =======================================
// VERSIONS FULL POLYHEDON

/*! \brief This function update the representations of the linear relaxation by adding a new hyperplane.
 *
 * USING THE MERGING METHOD
 * \param H Hyperplane*. A pointer to the new hyperplane added to the representation.
 * \param ptsDiscardedByH Point*. A pointer to the point that is cut out of the lower bound set by the new hyperplane H.
 */
void LinearRelaxation::updatePolyhedronFull(Hyperplane* H, Point* ptsDiscardedByH) {

	/*if (iteration == 0) {
		std::cout << "\n\n------------------- \n\n";
		H->print();
		std::cout << "   -> cuts out : ";
		ptsDiscardedByH->print();
		std::cout << " at " << ptsDiscardedByH << "\n\n";
	}*/

	// compute the new points

	double lambda;
	std::list<Point*> N(0), D(0), E(0), M(0), R(0); // new pts, degenerate pts, points to explore, modified feasible points, extreme rays
	std::list<Point*>::iterator v, w, v1, v2;
	std::list<Point*>* adjacentVertex;
	Point* u, * newPts;
	E.push_back(ptsDiscardedByH);
	ptsDiscardedByH->becomesDiscarded();
	ptsDiscardedByH->becomesVisited();

	while (E.size() != 0) {

		u = E.front();
		E.pop_front();

		//if (call == 55) {
			//std::cout << "\n ----- ";
			//u->print();
			//std::cout << " is selected\n";
		//}

		adjacentVertex = u->get_adjList();
		w = adjacentVertex->begin();

		while (w != adjacentVertex->end()) {

			//if (call == 55) {
				//std::cout << " adj list is : \n";
				//for (v1 = adjacentVertex->begin(); v1 != adjacentVertex->end(); v1++) {
					//std::cout << "   -> " << *v1 << " : ";
					//(*v1)->print();
					//std::cout << "\n";
				//}
				//std::cout << " selected is : ";
				//(*w)->print();
				//std::cout << "\n\n";
			//}

			v = w;
			w++;

			//(*v)->print();
			if ((*v)->isDiscarded(H)) {
				if (!(*v)->isVisited()) {
					(*v)->becomesVisited();
					E.push_back(*v);
					//if (call == 55) {
						//(*v)->print();
						//std::cout << " is added to E\n";
					//}
				}
				//else if ((*v)->is_ray()) {
					//if (call == 55)
						//std::cout << " infeasible ray lol\n";
				//}
				//else {
					//if (call == 55)
						//std::cout << " already visited infeasibe vertex\n";
				//}
			}

			else if (!(*v)->isDegenerate()) {

				if (!(*v)->isModified())
					M.push_back(*v);
				(*v)->becomesModified();

				lambda = H->edgeIntersection(u, *v);
				//std::cout << "    -> l = " << lambda << "\n";
				if (abs(lambda - 1) <= 0 + EPS_LAMBDA) { // 0.00001
					//if (call == 55) {
						//u->print();
						//std::cout << " IS DEGENERATE !!!\n";
					//}
					u->becomesDegenerate(&N, H);
					D.push_back(u);
					M.push_back(*v);
					(*v)->becomesModified();
					//u->notifyRays(H, M);
					break; // the unfeasible point is actually degenerate, we can leave the adjacent-vertex loop since we will only generate duplicates
				}
				else if ((abs(lambda) <= 0 + EPS_LAMBDA)) { // 0.00001
					(*v)->becomesDegenerate(&N, H);
					D.push_back(*v);
					//(*v)->notifyRays(H, M);
					//if (call == 55) {
						//(*v)->print();
						//std::cout << " IS DEGENERATE !\n";
					//}
				}
				else {
					newPts = new Point(u, *v, lambda, H);
					//if (call == 55) {
						//newPts->print();
						//std::cout << " is created";
					//}
					N.push_back(newPts);
				}
			}
			//if (call == 55) {
			//	//(*v)->print();
			//	//std::cout << " ????????? \n\n";
			//	std::cout << "\n";
			//}
		}
	}

	// clear the discarded points

	clearDiscardedPoints();

	// clear the non-facet-defining hyperplane

	clearRedundantFacets();

	// connect the new vertices and the degenerate vertices (update adjacency lists

	updateAjacencyLists(N, D); // extrPoints D

	// clear status of the concerned points before the next iteration

	clearStatus(N, R, M);

	// add the new facet and the new points & rays to the polyhedron description

	facets.push_back(H);
	for (v = N.begin(); v != N.end(); v++)
		extrPoints.push_back(*v);
	for (v = R.begin(); v != R.end(); v++)
		extrPoints.push_back(*v);

	//debug__resetAdjLists();


	//print();
	//if (call == 5)
		//std::cout << "stop\n";
	call++; // go to next iteration
}



/*! \brief Initialise the linear relaxation by building a simplex that contains it.
  *
  * \param linprog MathematicalModel*. A pointer to the problem this lower bound set refers to.
  */
void LinearRelaxation::initializeFull() { //MathematicalModel& lp

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

	pt = new Point(antiIdealPoint);
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


/*! \brief Return true if the point should be checked.
 *
 * This function checks whether (1) the point has already been checked by Cplex and (2) if not, if it is weakly dominated by
 * another existing point. If none of these conditions are met, the point should be checked for feasibility and the function
 * returns true.
 * \param y Point*. A pointer to the point being tested.
 * \return true if the point should be checked.
 */
bool LinearRelaxation::shouldBeCheckedFull(Point* y) {

	// if the point y has already been checked by cplex or if it is weakly dominated by another point, there is no need to check it

	if (y->isCheckedByCplex() || isWeeklyDominated(y)) // 
		return false;
	else
		return true;
}



/*! \brief This function computes the linear relaxation
 */
void LinearRelaxation::computeFull() {

	call++;
	try {
		if (!isFeasible()) {
			status = INFEASIBLE;
		}
		else {

			if (!warmstarted) {
				S->timeInitialization.StartTimer();
				initializeFull();
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
				else if (shouldBeCheckedFull(*currentPoint)) {

					//if (iteration == 1) { // 17239
					//	std::cout << "\n\n ררררררררררררררררררררררר\n point checked : ";
					//	(*currentPoint)->print();
					//	std::cout << " (at " << *currentPoint << ") ";
					//}

					if ((*currentPoint)->isNew() || (*currentPoint)->get_nbVar() == 0) {
						S->timeFeasibilityCheck.StartTimer();
						feasible = feasibilityCheck->solve(*(*currentPoint)->get_objVector(), iteration);
						S->timeFeasibilityCheck.StopTimer();
						std::cout << " Cplex called on : ";
						(*currentPoint)->print();
						std::cout << "\n";
						S->lpSolved++;
					}
					else {
						if ((*currentPoint)->satisfyBranchingDecisions(branchDec)) {
							feasible = true;
						}
						else { // check for alternative pre-images
							S->timeFeasibilityCheck.StartTimer();
							feasible = feasibilityCheck->solve(*(*currentPoint)->get_objVector(), iteration);
							S->timeFeasibilityCheck.StopTimer();
							std::cout << " Cplex called on : ";
							(*currentPoint)->print();
							std::cout << "\n";
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
						updatePolyhedronFull(H, *currentPoint);
						S->timeUpdatePolyhedron.StopTimer();

						// if the previous checkpoint has been destroyed, search for a new one

						if (checkPointDestroyed) {
							checkpoint = extrPoints.begin();
							while (checkpoint != extrPoints.end() && !shouldBeCheckedFull(*checkpoint)) { //(*checkpoint)->isCheckedByCplex
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
			std::cout << " \n LP solved : " << S->lpSolved << "\n";
			/*if (iteration == 1) {
				print();
				std::cout << "oof\n";
			}*/
			//print();
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







// ======================================================
// OLDER VERSIONS

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
	std::vector<double> sol;
	std::list<Point*> N(0), D(0), E(0), M(0), R(0); // new pts, degenerate pts, points to explore, modified feasible points, extreme rays
	std::list<Point*>::iterator v, w, v1, v2;
	std::list<Point*>* adjacentVertex;
	Point* u, * newPts;
	E.push_back(ptsDiscardedByH);
	ptsDiscardedByH->becomesDiscarded();
	ptsDiscardedByH->becomesVisited();
	int itv = 0;

	while (E.size() != 0) {

		u = E.front();
		E.pop_front();

		if (call <= 1) {
			std::cout << " \n\n -------- \n";
			u->print();
			std::cout << " is selected\n";
		}

		adjacentVertex = u->get_adjList();
		w = adjacentVertex->begin();

		while (w != adjacentVertex->end()) {

			v = w;
			w++;

			if (call <= 0) {
				std::cout << "   -> ";
				(*v)->print();
			}

			if ((*v)->isDiscarded(H)) {
				if (!(*v)->isVisited() && !(*v)->is_ray()) {
					(*v)->becomesVisited();
					E.push_back(*v);
					if (call <= 0) {
						std::cout << " is added to E";
					}
				}
			}

			else if (!(*v)->isDegenerate()) {

				if (!(*v)->isModified())
					M.push_back(*v);
				(*v)->becomesModified();

				//lambda = H->edgeIntersection(u, *v);
				sol = H->edgeIntersection2(u, *v);

				if (u->isVeryCloseTo(sol)) { // 0.00001 // abs(lambda - 1) <= 0 + EPS_LAMBDA
					u->becomesDegenerate(&N, H);
					D.push_back(u);
					M.push_back(*v);
					(*v)->becomesModified();
					u->notifyRays(H, M);

					if (call <= 0) {
						std::cout << " has made u degenerate";
					}
					break; // the unfeasible point is actually degenerate, we can leave the adjacent-vertex loop since we will only generate duplicates
				}
				else if ((*v)->isVeryCloseTo(sol) && !(*v)->is_ray()) { // 0.00001 // (abs(lambda) <= 0 + EPS_LAMBDA)
					(*v)->becomesDegenerate(&N, H);
					D.push_back(*v);
					(*v)->notifyRays(H, M);
					if (call <= 0) {
						std::cout << " is degenerate";
					}
				}
				else {
					newPts = new Point(u, *v, sol, H);
					N.push_back(newPts);
					if (call <= 0) {
						std::cout << "  => ";
						newPts->print();
						std::cout << " is created";
					}
				}
			}
			if (call <= 0) {
				std::cout << "\n";
			}
		}
	}


	// generate extreme rays

	for (int k = 0; k < lp->get_p(); k++) {
		if (H->get_normalVector(k) == 0) { // unbounded dimension    // H->get_normalVector(k) == 0

			// generate rays for the new points
			for (v = N.begin(); v != N.end(); v++) {
				if (!(*v)->has_ray(k)) {
					//std::cout << " ... generation on ";
					//(*v)->print();
					//std::cout << " in direction " << k << "\n";
					generateRay(*v, k, R);
				}
			}

			// generate rays for the degenerate vertices
			for (v = D.begin(); v != D.end(); v++) {
				if (!(*v)->has_ray(k)) {
					//std::cout << " ... generation on ";
					//(*v)->print();
					//std::cout << " in direction " << k << "\n";
					generateRay(*v, k, R);
				}
			}
		}
	}

	// clear the discarded points

	clearDiscardedPoints();

	// clear the non-facet-defining hyperplane

	clearRedundantFacets();

	// connect the new vertices and the degenerate vertices (update adjacency lists

	updateAjacencyLists(N, D); // extrPoints D

	// clear status of the concerned points before the next iteration

	clearStatus(N, R, M);

	// add the new facet and the new points & rays to the polyhedron description

	facets.push_back(H);
	for (v = N.begin(); v != N.end(); v++)
		extrPoints.push_back(*v);
	for (v = R.begin(); v != R.end(); v++)
		extrPoints.push_back(*v);

	//debug__resetAdjLists();

	// debug

	std::list<Hyperplane*>* hpp1, * hpp2;
	std::list<Hyperplane*>::iterator f1, f2;
	bool inIntersection, allIn = true;
	for (v1 = extrPoints.begin(); v1 != extrPoints.end(); v1++) {

		// debug 1
		for (v2 = extrPoints.begin(); v2 != extrPoints.end(); v2++) {
			if (*v1 != *v2 && (*v1)->isVeryCloseTo(*v2, call)) {
				(*v1)->print();
				std::cout << " is duplicate to ";
				(*v2)->print();
				std::cout << "\n";
				//print();
				std::cout << " duplicates ??\n";
			}
		}

		// debug 2
		hpp1 = (*v1)->get_activeHyperplanes();
		for (v2 = extrPoints.begin(); v2 != extrPoints.end(); v2++) {
			if (*v1 != *v2 && !(*v1)->is_ray() && !(*v2)->is_ray()) {
				hpp2 = (*v2)->get_activeHyperplanes();
				allIn = true;
				for (f1 = hpp1->begin(); f1 != hpp1->end(); f1++) {
					inIntersection = false;
					for (f2 = hpp2->begin(); f2 != hpp2->end(); f2++) {
						if (*f1 == *f2) {
							inIntersection = true;
							break;
						}
					}
					if (!inIntersection) {
						allIn = false;
						break;
					}
				}
				if (allIn) {
					(*v1)->print();
					std::cout << " is in the face of ";
					(*v2)->print();
					std::cout << "\n inclusion of facet-definition !!\n";
				}
			}
		}

		// debug 3
		if ((*v1)->get_adjList()->size() > lp->get_p())
			std::cout << " greater adj list !\n";
		else if (!(*v1)->is_ray() && (*v1)->get_adjList()->size() < lp->get_p())
			std::cout << " lower adj list !!\n";
	}


	//print();
	//if (call == 10)
		//std::cout << "stop\n";
	if (call == 19) {
		print();
		exportIteration();
		std::cout << "out";
	}
	call++; // go to next iteration
}




int LinearRelaxation::get_numberFacets() {
	return facets.size();
}

int LinearRelaxation::get_numberVertices() {
	return extrPoints.size();
}