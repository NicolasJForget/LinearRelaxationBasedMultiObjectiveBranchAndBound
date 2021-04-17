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
		pts[k] = new Point(antiIdealPoint, true);
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
	for (int k = 0; k < lp->get_p(); k++) {
		extrPoints.push_back(pts[k]);
		pts[k]->print();
		std::cout << "\n";
	}
}



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

		// adjust the point used for computation
		for (int l = 0; l < lp->get_p(); l++) {
			if (k == l) {
				refPt[l] = interiorPoint[l];
			}
			else {
				refPt[l] = antiIdealPoint[l];
			}
		}

		// call Cplex to compute the hyperplane
		S->lpSolved++;
		dummy = feasibilityCheck->solve(refPt);
		normalVector = feasibilityCheck->extractNormalVector();
		rhs = feasibilityCheck->extractConstant(*lp, branchDec);
		y[k] = rhs;
		H[k] = new Hyperplane(normalVector, rhs);
		facets.push_back(H[k]);
		//H[k]->print();
	}

	// compute the 2^p + 1 first extreme points

	Point* ePt = new Point(y, false);
	std::vector<Point*> eRay(lp->get_p());
	std::vector<double> ray(lp->get_p());

	for (int k = 0; k < lp->get_p(); k++) {

		// compute the coordinates of the ray
		for (int l = 0; l < lp->get_p(); l++) {
			if (k == l) {
				ray[l] = antiIdealPoint[l];
			}
			else {
				ray[l] = y[l];
			}
		}
		eRay[k] = new Point(ray, true);

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

	// compute the bounding hyperplanes

	for (int k = 0; k < lp->get_p(); k++) {
		for (int l = 0; l < lp->get_p(); l++) {
			if (k == l) normalVector[l] = 1;
			else normalVector[l] = 0;
		}
		rhs = antiIdealPoint[k];
		hpp = new Hyperplane(normalVector, rhs);
		boundingBox.push_back(hpp);
	}

	// compute adjacency lists of the extreme point and rays

	ePt->receive_ray();
	for (int k = 0; k < lp->get_p(); k++) {
		// adjacency with the unique extreme point
		ePt->addAdjacentPoint(eRay[k]);
		eRay[k]->addAdjacentPoint(ePt);

		// adjecency between extreme rays & adjust value of unbounded directions to meet the bounding hyperplane
		for (int l = 0; l < lp->get_p(); l++) {
			if (k != l) eRay[k]->addAdjacentPoint(eRay[l]);
			//if (k != l) eRay[k]->addAdjacentPoint(eRay[l]);
			//else eRay[k]->setObjVector(k, eRay[k]->findMissingCoordinate(*hpp, l));
		}

		// add the bounding box the the active hyperplanes of the rays
		eRay[k]->addActiveHyperplane(boundingBox[k]);
		boundingBox[k]->addVertex(eRay[k]);
	}

	// add the extreme points to the list extrPoints

	ePt->becomesNonDegenerate();
	extrPoints.push_back(ePt);
	for (int k = 0; k < lp->get_p(); k++) {
		eRay[k]->becomesNonDegenerate();
		extrPoints.push_back(eRay[k]);
	}

	//print();
	//std::cout << "lol";
}

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

		// adjust the point used for computation
		for (int l = 0; l < lp->get_p(); l++) {
			if (k == l) {
				refPt[l] = interiorPoint[l];
			}
			else {
				refPt[l] = antiIdealPoint[l];
			}
		}

		// call Cplex to compute the hyperplane
		S->lpSolved++;
		dummy = feasibilityCheck->solve(refPt);
		normalVector = feasibilityCheck->extractNormalVector();
		rhs = feasibilityCheck->extractConstant(*lp, branchDec);
		y[k] = rhs;
		H[k] = new Hyperplane(normalVector, rhs);
		facets.push_back(H[k]);
		//H[k]->print();
	}

	// compute the 2^p + 1 first extreme points

	Point* ePt = new Point(y, false);
	std::vector<Point*> eRay(lp->get_p());
	std::vector<double> ray(lp->get_p());

	for (int k = 0; k < lp->get_p(); k++) {

		// compute the coordinates of the ray
		for (int l = 0; l < lp->get_p(); l++) {
			if (k == l) {
				ray[l] = antiIdealPoint[l];
			}
			else {
				ray[l] = y[l];
			}
		}
		eRay[k] = new Point(ray, true);

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

	// compute the bounding hyperplanes

	for (int k = 0; k < lp->get_p(); k++) {
		for (int l = 0; l < lp->get_p(); l++) {
			if (k == l) normalVector[l] = 1;
			else normalVector[l] = 0;
		}
		rhs = antiIdealPoint[k];
		hpp = new Hyperplane(normalVector, rhs);
		boundingBox.push_back(hpp);
	}

	// compute adjacency lists of the extreme point and rays

	ePt->receive_ray();
	for (int k = 0; k < lp->get_p(); k++) {
		// adjacency with the unique extreme point
		ePt->addAdjacentPoint(eRay[k]);
		eRay[k]->addAdjacentPoint(ePt);

		// adjecency between extreme rays & adjust value of unbounded directions to meet the bounding hyperplane
		for (int l = 0; l < lp->get_p(); l++) {
			if (k != l) eRay[k]->addAdjacentPoint(eRay[l]);
			//if (k != l) eRay[k]->addAdjacentPoint(eRay[l]);
			//else eRay[k]->setObjVector(k, eRay[k]->findMissingCoordinate(*hpp, l));
		}

		// add the bounding box the the active hyperplanes of the rays
		eRay[k]->addActiveHyperplane(boundingBox[k]);
		boundingBox[k]->addVertex(eRay[k]);
	}

	// add the extreme points to the list extrPoints

	ePt->becomesNonDegenerate();
	extrPoints.push_back(ePt);
	for (int k = 0; k < lp->get_p(); k++) {
		eRay[k]->becomesNonDegenerate();
		extrPoints.push_back(eRay[k]);
	}

	//print();
	//std::cout << "lol";
}




// compute the bounding hyperplanes

double sum;
rhs = -100000000;
for (int k = 0; k < lp->get_p(); k++) {
	sum = 0;
	normalVector[k] = 1;
	for (int l = 0; l < lp->get_p(); l++) {
		if (k == l) sum += ePt->get_objVector(l);
		else sum += antiIdealPoint[l];
	}
	if (sum > rhs) {
		rhs = sum;
	}
}
hpp = new Hyperplane(normalVector, rhs);
boundingBox.push_back(hpp);

void LinearRelaxation::updatePolyhedron2(Hyperplane* H, Point* ptsDiscardedByH) {

	std::list<Point*>::iterator vertex;
	std::list<Point*>::iterator vertex2;
	std::list<Point*>::iterator vertex3;
	std::list <Hyperplane*>::iterator f;
	std::list <Hyperplane*>::iterator fCurrent;
	std::list<Hyperplane*>* hpp;
	std::list<Point*>* adjacentVertex;
	Point* newPts;
	Point* newRay;
	std::vector<Point*> allNewVertices(0);
	std::vector<Point*> allNewRays(0);
	std::vector<Point*> allDegenerateVertices(0);
	std::vector<int> unboundedDim(0);
	std::vector<int> boundedDim(0);
	std::list<Point*>* av;
	std::vector<double> y(lp->get_p());
	bool hppIsBoxDefining = false;

	/*if (iteration == 0) {
		std::cout << "\n\n------------------- \n\n";
		H->print();
	}*/

	// check the unbounded dimensions of this hyperplane for the extreme rays

	for (int k = 0; k < lp->get_p(); k++) {
		//std::cout << " dimension " << k << " is ";
		if (H->get_normalVector(k) == 0) {
			unboundedDim.push_back(k);
			//std::cout << "unbounded\n";
		}
		else {
			boundedDim.push_back(k);
			//std::cout << "bounded\n";
		}
	}
	if (boundedDim.size() == 1) {
		hppIsBoxDefining = true;
	}

	// determine the status of each extreme point

	vertex = extrPoints.begin();
	while (vertex != extrPoints.end()) {
		(*vertex)->becomesNonDegenerate();
		if (*vertex == ptsDiscardedByH || (*vertex)->below2(*H)) { // we are sure to at least discard the point cut out by H
			(*vertex)->becomesDiscarded();
			/*if (iteration == 0) {
				(*vertex)->print();
				std::cout << " is discarded\n";
			}*/
		}
		else {
			/*if (iteration == 0) {
				(*vertex)->print();
				std::cout << " is not discarded\n";
			}*/
		}
		++vertex;
	}
	//if (iteration == 1)	std::cout << "\n";

	// compute the new extreme points

	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); vertex++) {
		if ((*vertex)->isDiscarded()) { // we have a discarded vertex

			adjacentVertex = (*vertex)->get_adjList();
			for (vertex2 = adjacentVertex->begin(); vertex2 != adjacentVertex->end(); vertex2++) {
				if (!(*vertex2)->isDiscarded() && !(*vertex2)->isDegenerate() && !(*vertex)->isDegenerate()) { // we have a non-discarded adjacent vertex

					newPts = new Point((*vertex2)->edgeIntersection(**vertex, *H));

					// is the new point degenerate ?
					if (newPts->isVeryCloseTo(*vertex)) { // possible opti : go out of vertex2 loop here
						(*vertex)->becomesNonDiscarded(); // we don't want to delete it anymore
						(*vertex)->becomesDegenerate();
						(*vertex)->addActiveHyperplane(H);
						H->addVertex(*vertex);
						if (!(*vertex)->is_ray()) allDegenerateVertices.push_back(*vertex);
						delete newPts;
					}
					else if (newPts->isVeryCloseTo(*vertex2)) {
						(*vertex2)->becomesDegenerate();
						(*vertex2)->addActiveHyperplane(H);
						H->addVertex(*vertex2);
						if (!(*vertex2)->is_ray()) allDegenerateVertices.push_back(*vertex2);
						delete newPts;
					}
					else {

						// create the new vertex 
						(*vertex2)->replaceAdjVertex(*vertex, newPts);
						newPts->addAdjacentPoint(*vertex2); // we know that at least
						newPts->updateActiveHyperplanes(**vertex, **vertex2, H); // update the active hyperplanes of the new point
						allNewVertices.push_back(newPts);

						/*newPts->print();
						std::cout << " is created (point)\n";*/

						// check the extreme rays
						// computed outside of the loop
					}
				}
			}

		}
	}

	/*if (iteration == 0) {
		for (vertex = extrPoints.begin(); vertex != extrPoints.end(); vertex++) {
			if ((*vertex)->isDegenerate()) {
				(*vertex)->print();
				std::cout << " is degenerate\n";
			}
		}
	}*/

	// compute extreme rays

	double maxVal;
	double maxVal2;
	Point* maxPt = NULL;
	bool linkedToDegenerate;
	std::vector<bool> existingDirections;
	int ctr;
	std::vector<Point*> attachedPts;

	if (hppIsBoxDefining) {
		for (auto l : boundedDim) {
			//maxVal
		}
	}

	//std::cout << "\n\n ----------- \n";
	for (auto k : unboundedDim) {
		//std::cout << " -> unbounded dimension : " << k << "\n";
		attachedPts = std::vector<Point*>(0);
		for (auto l : boundedDim) {
			//std::cout << "     . search for max value in direction " << l << "\n";

			// search for max
			maxVal = interiorPoint[l] - 1;
			maxVal2 = interiorPoint[k] - 1;
			linkedToDegenerate = false;
			for (auto vtx : allNewVertices) { // among the new vertices...
				//vtx->print();
				//std::cout << " is tested ...";
				if (vtx->get_objVector(l) > maxVal) {
					maxPt = vtx;
					maxVal = maxPt->get_objVector(l);
					maxVal2 = maxPt->get_objVector(k);
					//std::cout << " and accepted !";
				}
				else if (vtx->get_objVector(l) == maxVal) {
					if (vtx->get_objVector(k) >= maxVal2) {
						maxPt = vtx;
						maxVal2 = maxPt->get_objVector(k);
						//std::cout << " and accepted !";
					}
				}
				//std::cout << "\n";
			}
			for (auto vtx : allDegenerateVertices) { // ...but also among the degenerate vertices (and not rays) !
				//vtx->print();
				//std::cout << " is tested ...";
				if (vtx->get_objVector(l) > maxVal) {
					maxPt = vtx;
					maxVal = maxPt->get_objVector(l);
					maxVal2 = maxPt->get_objVector(k);
					linkedToDegenerate = true;
					//std::cout << " and accepted !";
				}
				else if (vtx->get_objVector(l) == maxVal) {
					if (vtx->get_objVector(k) >= maxVal2) {
						maxPt = vtx;
						maxVal2 = maxPt->get_objVector(k);
						linkedToDegenerate = true;
						//std::cout << " and accepted !";
					}
				}
				//std::cout << "\n";
			}

			// it may happen that a ray in one of the unbounded dimensions already exist at a degenerate point
			if (linkedToDegenerate) {
				existingDirections = maxPt->getRayDirections(antiIdealPoint);
			}

			if (!linkedToDegenerate || existingDirections[l] == false) {
				// compute coordinates of the ray
				for (int t = 0; t < lp->get_p(); t++) {
					if (t == k) y[t] = antiIdealPoint[t]; // infinite in the unbounded direction we are currently looking at
					else y[t] = maxPt->get_objVector(t);
				}
				newRay = new Point(y, true);

				//newRay->print();
				//std::cout << " is generated\n\n";

				// connect the new ray
				maxPt->receive_ray();
				maxPt->addAdjacentPoint(newRay);
				newRay->addAdjacentPoint(maxPt);
				newRay->updateActiveHyperplanes(*maxPt, H, l);

				allNewRays.push_back(newRay);
			}
		}
	}


	// clear the discarded points

	vertex2 = extrPoints.begin();
	while (vertex2 != extrPoints.end()) {
		vertex = vertex2;
		vertex2++;

		// if an extreme ray is adjacent to a degenerate vertex, it is also degenerate
		if ((*vertex)->is_ray() && (*vertex)->isDiscarded()) {
			adjacentVertex = (*vertex)->get_adjList();
			if (adjacentVertex->size() != 0 && (*adjacentVertex->begin())->isDegenerate()) {
				(*vertex)->becomesNonDiscarded();
				(*vertex)->becomesDegenerate();
				//std::cout << " one saved\n";
			}
		}

		// discard the vertex
		if ((*vertex)->isDiscarded()) {
			(*vertex)->notifyDeletion();
			if ((*vertex)->isCheckpoint()) checkPointDestroyed = true;
			delete* vertex;
			extrPoints.erase(vertex);
		}
	}

	// clear the non-facet-defining hyperplane

	f = facets.begin();
	while (f != facets.end()) {
		fCurrent = f;
		f++;
		(*fCurrent)->checkRedundancy();
		if ((*fCurrent)->isRedundant()) {
			(*fCurrent)->notifyDeletion();
			delete* fCurrent;
			facets.erase(fCurrent);
		}
	}

	// connect the new vertices and the degenerate vertices (update adjacency lists)

	if (allNewVertices.size() >= 2) {
		for (int i = 0; i < allNewVertices.size() - 1; i++) { // new vertices together
			for (int j = i + 1; j < allNewVertices.size(); j++) {
				if (allNewVertices[i]->sizeIntersectionActiveHyperplanes(*allNewVertices[j]) >= lp->get_p() - 1) {
					allNewVertices[i]->addAdjacentPoint(allNewVertices[j]);
					allNewVertices[j]->addAdjacentPoint(allNewVertices[i]);
				}
			}
		}
	}

	if (allDegenerateVertices.size() >= 1 && allNewVertices.size() >= 1) {
		for (int j = 0; j < allDegenerateVertices.size(); j++) { // degenerate vertices with new vertices
			for (int i = 0; i < allNewVertices.size(); i++) {
				if (allDegenerateVertices[j]->sizeIntersectionActiveHyperplanes(*allNewVertices[i]) >= lp->get_p() - 1) {
					allDegenerateVertices[j]->addAdjacentPoint(allNewVertices[i]);
					allNewVertices[i]->addAdjacentPoint(allDegenerateVertices[j]);
				}
			}
		}
	}

	if (allDegenerateVertices.size() >= 2) {
		for (int i = 0; i < allDegenerateVertices.size() - 1; i++) { // degenerate vertices together
			for (int j = i + 1; j < allDegenerateVertices.size(); j++) {
				if (allDegenerateVertices[i]->sizeIntersectionActiveHyperplanes(*allDegenerateVertices[j]) >= lp->get_p() - 1) {
					allDegenerateVertices[i]->addAdjacentPoint(allDegenerateVertices[j]);
					allDegenerateVertices[j]->addAdjacentPoint(allDegenerateVertices[i]);
				}
			}
		}
	}

	// add the new facet and the new points & rays to the polyhedron description

	facets.push_back(H);
	for (int i = 0; i < allNewVertices.size(); i++) {
		allNewVertices[i]->becomesNonDegenerate();
		extrPoints.push_back(allNewVertices[i]);
	}
	for (int i = 0; i < allNewRays.size(); i++) {
		allNewRays[i]->becomesNonDegenerate();
		extrPoints.push_back(allNewRays[i]);
	}

	//print();

	/*if (tamerlaput) {
		std::cout << " non ? \n";
	}*/

	//std::cout << "ok lol\n";
	//std::cout << "tmr";
}



for (auto k : boundedDim) {

	//std::cout << "\n -------------\n search for new point:\n";
	// for each bounded dimension, we search for the most extreme values
	maxVal = interiorPoint[k] - 1;
	linkedToDegenerate = false;
	for (auto vtx : allNewVertices) { // among the new vertices...
		if (vtx->get_objVector(k) >= maxVal) {
			maxPt = vtx;
			maxVal = maxPt->get_objVector(k);
		}
	}
	for (auto vtx : allDegenerateVertices) { // ...but also among the degenerate vertices (and not rays) !
		if (vtx->get_objVector(k) >= maxVal) { // !vtx->is_ray() && 
			maxPt = vtx;
			maxVal = maxPt->get_objVector(k);
			linkedToDegenerate = true;
		}
	}

	// filtering out unecessary

	if (linkedToDegenerate) { // it may happen that a ray in one of the unbounded dimensions already exist at a degenerate point
		existingDirections = maxPt->getRayDirections(antiIdealPoint);
	}

	// we now attach one extreme ray for each unbounded dimension at the point found
	for (auto l : unboundedDim) {
		if (!linkedToDegenerate || existingDirections[l] == false) {

			// compute the coordinates of the new ray
			for (int t = 0; t < lp->get_p(); t++) {
				if (t == l) y[t] = antiIdealPoint[t];
				else y[t] = maxPt->get_objVector(t);
			}
			newRay = new Point(y, true);

			// connect the new ray
			maxPt->receive_ray();
			maxPt->addAdjacentPoint(newRay);
			newRay->addAdjacentPoint(maxPt);
			newRay->updateActiveHyperplanes(*maxPt, H, l);

			allNewRays.push_back(newRay);
			newRay->print();
			std::cout << " is created (ray) at ";
			maxPt->print();
			std::cout << "\n";
		}
	}
}






if ((*vertex)->has_ray()) {
	for (vertex3 = adjacentVertex->begin(); vertex3 != adjacentVertex->end(); vertex3++) {
		if ((*vertex3)->isDiscarded() && (*vertex3)->is_ray()) { // we have a discarded extreme ray

			// compute the coordinates of the new ray
			for (int k = 0; k < lp->get_p(); k++) {
				if ((*vertex3)->get_objVector(k) == antiIdealPoint[k]) { // the unbounded dimension
					y[k] = antiIdealPoint[k];
				}
				else {
					y[k] = newPts->get_objVector(k);
				}
			}
			newRay = new Point(y, true);

			// connect the new ray
			newPts->receive_ray();
			newPts->addAdjacentPoint(newRay);
			newRay->addAdjacentPoint(newPts);
			newRay->updateActiveHyperplanes(**vertex3, *newPts, H);

			// an extreme ray should only be adjacent to one point, and we already know which one,
			// so there is no need to add the new ray to the list of allNewVerticies.
			allNewRays.push_back(newRay);
		}
	}
}




/*! \brief This function update the representations of the linear relaxation by adding a new hyperplane.
 *
 * USING THE MERGING METHOD
 * \param H Hyperplane*. A pointer to the new hyperplane added to the representation.
 * \param ptsDiscardedByH Point*. A pointer to the point that is cut out of the lower bound set by the new hyperplane H.
 */
void LinearRelaxation::updatePolyhedron2(Hyperplane* H, Point* ptsDiscardedByH) {

	std::list<Point*>::iterator vertex;
	std::list<Point*>::iterator vertex2;
	std::list <Hyperplane*>::iterator f;
	std::list <Hyperplane*>::iterator fCurrent;
	std::list<Hyperplane*>* fln;
	std::list<Point*>* adjacentVertex;

	int nbNewPts = 0;

	if (iteration == 0) {
		std::cout << "\n\n------------------- \n\n";
		H->print();
	}

	/*if (iteration == 1) {
		std::cout << "STAUHP LAUWL";
	}*/

	//facets.push_back(H);
	bool stop = false; // for debugging

	// search for infeasible and degenerated points
	vertex = extrPoints.begin();
	while (vertex != extrPoints.end()) {
		(*vertex)->becomesNonDegenerate();
		if (*vertex == ptsDiscardedByH) { // we are sure to discard the infeasible point
			(*vertex)->becomesDiscarded();
			(*vertex)->print();
			std::cout << " is discarded\n";
		}
		else if ((*vertex)->below2(*H)) { // !(*vertex)->isCheckpoint() && 
			(*vertex)->becomesDiscarded();
			(*vertex)->print();
			std::cout << " is discarded\n";
		}
		++vertex;
	}
	std::cout << "\n";

	// search for new extreme points (comparison with edges, by adjacency lists)
	Point* newPts;
	std::vector<Point*> allNewVertices(0);
	std::list<Point*>* av;
	std::list<Point*>::iterator vertex3;

	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); ++vertex) { // for each extreme point
		if ((*vertex)->isDiscarded()) { // if it is discarded
			adjacentVertex = (*vertex)->get_adjList(); // we look at its adjacency list
			for (vertex2 = adjacentVertex->begin(); vertex2 != adjacentVertex->end(); ++vertex2) {
				if (!(*vertex2)->isDiscarded()) { // for each non-discarded adjacent vertex, we compute the new point

					newPts = new Point((*vertex2)->edgeIntersection(**vertex, *H)); // new pts as intersection of H and edge
					if (iteration == 0) {
						newPts->print();
						std::cout << " inherited from ";
						(*vertex)->print();
						std::cout << " and ";
						(*vertex2)->print();
						std::cout << " at it " << nbNewPts << "\n";

						av = newPts->get_adjList();
						for (vertex3 = av->begin(); vertex3 != av->end(); vertex3++) {
							std::cout << "    -> ";
							(*vertex3)->print();
							std::cout << "\n";
						}
					}

					(*vertex2)->replaceAdjVertex((*vertex)->get_adress(), newPts); // update adj vertex for feasible one
					newPts->addAdjacentPoint((*vertex2)->get_adress()); // init adj vertex for new one
					newPts->updateActiveHyperplanes(**vertex, **vertex2, H); // update the active hyperplanes for the new point // , H
					//extrPoints.push_back(newPts); //new pts added to list of extreme points => DONE LATER
					nbNewPts++;
					allNewVertices.push_back(newPts); // remember the new vertex

				}
			}
		}
	}

	int s = allNewVertices.size();
	std::list<Hyperplane*>* hpp;
	bool alreadyAdj = false;
	std::list <Hyperplane*>* listHpp;
	std::list<Point*>* defPtsHpp;

	// filter the NEW points very close to each other and destroy the duplicates
	int i = 0;
	int j;
	while (i < s - 1) {
		j = i + 1;
		while (j < s) {
			if (allNewVertices[i]->isVeryCloseTo(allNewVertices[j])) {

				if (iteration == 0) {
					std::cout << "\n";
					allNewVertices[i]->print();
					std::cout << " is used a replacement for ";
					allNewVertices[j]->print();
				}

				// merge the two points: i is kept while j will be deleted
				allNewVertices[i]->merge(allNewVertices[j]);

				// notify active hyperplane that j will be deleted
				listHpp = allNewVertices[j]->get_activeHyperplanes();
				for (f = listHpp->begin(); f != listHpp->end(); f++) {
					(*f)->removeVertex(allNewVertices[j]);
				}

				// notify the adjacent vertex of j that j is replaced by i
				//std::cout << "\n   -> with adj list:\n";
				adjacentVertex = allNewVertices[j]->get_adjList();
				for (vertex = adjacentVertex->begin(); vertex != adjacentVertex->end(); vertex++) {
					/*(*vertex)->print();
					std::cout << "\n";*/
					(*vertex)->replaceAdjVertex(allNewVertices[j], allNewVertices[i]);
				}

				// delete j
				delete allNewVertices[j];
				allNewVertices.erase(allNewVertices.begin() + j);
				j -= 1;
				s -= 1;
			}
			j++;
		}
		i++;
	}

	// Filter the new points with the old points marked as degenerate, meaning that there is a duplicate of the degenerate point
	s = allNewVertices.size();
	std::list<Point*>* nik;
	std::list<Point*>::iterator tamer;
	bool mergedAtLeastOnce;
	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); vertex++) {
		if (true) { //(*vertex)->isDegenerate()
			j = 0;
			mergedAtLeastOnce = false;
			while (j < s) {
				if ((*vertex)->isVeryCloseTo(allNewVertices[j])) {

					if (iteration == 0) {
						std::cout << "\n";
						(*vertex)->print();
						std::cout << " (fixed pts) is used a replacement for ";
						allNewVertices[j]->print();
					}

					// merge the two points: vertex is kept while j will be deleted
					(*vertex)->merge(allNewVertices[j]);

					// notify active hyperplane that j will be deleted
					listHpp = allNewVertices[j]->get_activeHyperplanes();
					for (f = listHpp->begin(); f != listHpp->end(); f++) {
						(*f)->removeVertex(allNewVertices[j]);
					}

					// notify the adjacent vertex of j that j is replaced by i
					//std::cout << "\n   -> with adj list:\n";
					adjacentVertex = allNewVertices[j]->get_adjList();
					for (vertex2 = adjacentVertex->begin(); vertex2 != adjacentVertex->end(); vertex2++) {
						/*(*vertex2)->print();
						std::cout << "\n";*/
						(*vertex2)->replaceAdjVertex(allNewVertices[j], *vertex);
					}

					/*adjacentVertex = (*vertex)->get_adjList();
					for (vertex2 = adjacentVertex->begin(); vertex2 != adjacentVertex->end(); vertex2++) {
						std::cout << "\n     -> ";
						(*vertex2)->print();
					}
					std::cout << "\n";*/

					// delete j
					delete allNewVertices[j];
					allNewVertices.erase(allNewVertices.begin() + j);
					j -= 1;
					s -= 1;

					// here we keep an old vertex, that did not generated a new point, meaning that its discarded
					// neighboor has not been replaced with a new point. Thus, we disconnect it to the other discarded vertices,
					// which will be deleted later
					//adjacentVertex = (*vertex)->get_adjList();
					//tamer = adjacentVertex->begin();
					//while (tamer != adjacentVertex->end()) {
					//	vertex2 = tamer;
					//	tamer++;
					//	if ((*vertex2)->isDiscarded()) { // && !(*vertex2)->isDegenerate()
					//		(*vertex)->removeAdjacentPoint(*vertex2);
					//	}
					//}
					mergedAtLeastOnce = true; // remember that we merged this old vertex at least once

					// debug
					/*std::cout << "\n ררררר \n";
					adjacentVertex = (*vertex)->get_adjList();
					for (vertex2 = adjacentVertex->begin(); vertex2 != adjacentVertex->end(); vertex2++) {
						std::cout << "\n New adj vertex : ";
						(*vertex2)->print();
						std::cout << "\n";
						nik = (*vertex2)->get_adjList();
						for (tamer = nik->begin(); tamer != nik->end(); tamer++) {
							std::cout << "      -> at " << *tamer << ", we have ";
							(*tamer)->print();
							std::cout << "\n";
						}
					}*/
				}
				j++;
			}

			// if the degenerate vertex is merged at least once, add it to the list of new vertices to make the connetions
			// for the adjacency lists later.
			if (mergedAtLeastOnce) {
				//allNewVertices.push_back(*vertex);
				(*vertex)->becomesDegenerate();
				(*vertex)->becomesNonDiscarded();
				/*(*vertex)->print();
				std::cout << " becomes degenerate\n";*/
			}
			//(*vertex)->becomesNonDegenerate();		
		}
		(*vertex)->purge(); // CORRECTION LAUNCHED HERE !!!!!
	}

	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); vertex++) {
		if ((*vertex)->isDegenerate()) {
			/*if (iteration == 1) {
				std::cout << "\n ---- degenerate vtx : ";
				(*vertex)->print();
				std::cout << "\n";
			}*/
			// here we keep an old vertex, that did not generated a new point, meaning that its discarded
			// neighboor has not been replaced with a new point. Thus, we disconnect it to the other discarded vertices,
			// which will be deleted later
			adjacentVertex = (*vertex)->get_adjList();
			tamer = adjacentVertex->begin();
			while (tamer != adjacentVertex->end()) {
				vertex2 = tamer;
				tamer++;
				/*if (iteration == 1) {
					std::cout << "   -> ";
					(*vertex2)->print();
				}*/
				if ((*vertex2)->isDiscarded() && !(*vertex2)->isDegenerate()) {
					(*vertex)->removeAdjacentPoint(vertex2);
					/*if (iteration == 1) {
						std::cout << " is discarded !";
					}*/
				}
				/*if (iteration == 1) {
					std::cout << "\n";
					if (tamer != adjacentVertex->end()) {
						std::cout << " (tamer is : ";
						(*tamer)->print();
						std::cout << " )\n";
					}
				}*/
			}
		}
	}

	// destroy the discarded points
	vertex = extrPoints.begin();
	while (vertex != extrPoints.end()) { // explore the extreme points
		//(*vertex)->becomesNonDegenerate(); // reset degenerancy status for each extreme point
		vertex2 = vertex;
		++vertex;
		if ((*vertex2)->isDiscarded()) { // if a vertex is discarded

			// if the discarded vertex was the checkpoint, notify that a new one should be found
			if ((*vertex2)->isCheckpoint()) {
				checkPointDestroyed = true;
			}

			// notify its active hyperplanes
			listHpp = (*vertex2)->get_activeHyperplanes();
			for (f = listHpp->begin(); f != listHpp->end(); f++) {
				(*f)->removeVertex(*vertex2);
			}

			// each discarded point should have been replaced by a feasible (for this iteration) point in adj lists
			// when building new points.
			/*if (iteration == 36) {
				std::cout << " Destroyed vertex at " << *vertex2 << " : ";
				(*vertex2)->print();
				std::cout << "\n";
			}*/
			//(*vertex2)->print();
			//std::cout << "\n";
			//nik = (*vertex2)->get_adjList();
			//for (tamer = nik->begin(); tamer != nik->end(); tamer++) {
			//	std::cout << "   -> ";
			//	//(*tamer)->print();
			//	std::cout << "\n";
			//}

			// destroy the point
			if (iteration == 0) {
				(*vertex2)->print();
				std::cout << " is deleted\n";
			}
			delete* vertex2;
			extrPoints.erase(vertex2);
		}
	}

	// updates the defining points of the active hyperplanes of the new points & add new points to the list of extreme points
	s = allNewVertices.size();
	for (i = 0; i < s; i++) {

		// update the hyperplanes
		listHpp = allNewVertices[i]->get_activeHyperplanes();
		for (f = listHpp->begin(); f != listHpp->end(); f++) {
			(*f)->addVertex(allNewVertices[i]);
		}

		// add the point to the list of extreme points

		/*if (!allNewVertices[i]->isDegenerate()) {
			extrPoints.push_back(allNewVertices[i]);
		}*/
		//extrPoints.push_back(allNewVertices[i]);
	}

	// pre-filter redundants hyperplanes
	for (f = facets.begin(); f != facets.end(); f++) {
		(*f)->checkRedundancy();
	}

	f = facets.begin();
	while (f != facets.end()) {
		fCurrent = f;
		++f;

		// quickly check redundancy (based on number of defining points only => necessary but not sufficient condition)
		if ((*fCurrent)->isRedundant()) { // (*fCurrent)->quickCheckRedundancy()

			// notify its defining points that this hyperplane is redundant
			defPtsHpp = (*fCurrent)->get_defPts();
			for (vertex = defPtsHpp->begin(); vertex != defPtsHpp->end(); vertex++) {
				(*vertex)->discardHyperplane(*fCurrent);
			}

			// destroy the hyperplane
			delete* fCurrent;
			facets.erase(fCurrent);
		}
	}

	/*for (vertex = extrPoints.begin(); vertex != extrPoints.end(); vertex++) {
		//(*vertex)->purge();
		adjacentVertex = (*vertex)->get_adjList();
		std::cout << "\n\n Point ";
		(*vertex)->print();
		std::cout << " (" << *vertex << ") with adj list:";
		for (vertex2 = adjacentVertex->begin(); vertex2 != adjacentVertex->end(); vertex2++) {
			std::cout << "\n      -> ";
			(*vertex2)->print();
			std::cout << " from " << *vertex2;
		}
	}*/

	// compute adjacency lists
	// May generate more adjacency than expected due to a partial but faster cleanup in the hyperplanes.
	// This is corrected when the lower bound set is fully computed.
	int cas;
	//std::cout << " \n New vertices : \n";
	i = 0;
	j = 0;
	s = allNewVertices.size();
	//for (int i = 0; i < s; i++) {
	while (i < s - 1) {
		j = i + 1;
		//allNewVertices[i]->print();
		//if (allNewVertices[i]->isDegenerate()) std::cout << " is degenerate";
		//std::cout << "\n";
		//for (int j = i + 1; j < s; j++) {
		while (j < s) {
			//if (abs(allNewVertices[i]->degreeOfBounding(antiIdealPoint) - allNewVertices[j]->degreeOfBounding(antiIdealPoint)) <= 1) {
			if (allNewVertices[i]->sizeIntersectionActiveHyperplanes(*allNewVertices[j]) >= H->get_dim()) {
				//if (cas == 0 || cas == 1) { // everything is ok, they are adjacent
				allNewVertices[i]->addAdjacentPoint(allNewVertices[j]);
				allNewVertices[j]->addAdjacentPoint(allNewVertices[i]);
				/*allNewVertices[i]->print();
				std::cout << " is adj to ";
				allNewVertices[j]->print();
				std::cout << "\n";
				hpp = allNewVertices[i]->get_activeHyperplanes();
				std::cout << "\n hpp list 1: \n";
				for (f = hpp->begin(); f != hpp->end(); f++) {
					std::cout << "    -> " << *f << " : ";
					(*f)->print();
				}
				hpp = allNewVertices[j]->get_activeHyperplanes();
				std::cout << "\n hpp list 2: \n";
				for (f = hpp->begin(); f != hpp->end(); f++) {
					std::cout << "    -> " << *f << " : ";
					(*f)->print();
				}*/
				cas = allNewVertices[i]->isOnEdge(allNewVertices[j]);
				//}
				//else if (cas == 2) { // allNewVertices[i] is located on an existing edge and should be corrected (i.e. discarded)
				//	// we notify its hyperplanes
				//	listHpp = allNewVertices[i]->get_activeHyperplanes();
				//	for (f = listHpp->begin(); f != listHpp->end(); f++) {
				//		(*f)->removeVertex(allNewVertices[i]);
				//	}
				//	// we notify its adjacent points
				//	adjacentVertex = allNewVertices[i]->get_adjList();
				//	for (vertex3 = adjacentVertex->begin(); vertex3 != adjacentVertex->end(); vertex3++) {
				//		(*vertex3)->replaceAdjVertex(allNewVertices[i], allNewVertices[j]);
				//	}
				//	delete allNewVertices[i];
				//	allNewVertices.erase(allNewVertices.begin() + i);
				//	i--;
				//	j--;
				//	s--;
				//}
				//else if (cas == 1) { // The adjacent point is in the middle of the edge defined by these two new vertices. Treated internally.
				//	allNewVertices[i]->addAdjacentPoint(allNewVertices[j]);
				//	allNewVertices[j]->addAdjacentPoint(allNewVertices[i]);
				//}
				//else if (cas == 3) { // allNewVertices[j] is in the interior of the edge defined by aNV[i] and its adjacent point.
				//	// we notify its hyperplanes
				//	listHpp = allNewVertices[j]->get_activeHyperplanes();
				//	for (f = listHpp->begin(); f != listHpp->end(); f++) {
				//		(*f)->removeVertex(allNewVertices[j]);
				//	}
				//	// we notify its adjacent points
				//	adjacentVertex = allNewVertices[j]->get_adjList();
				//	for (vertex3 = adjacentVertex->begin(); vertex3 != adjacentVertex->end(); vertex3++) {
				//		(*vertex3)->removeAdjacentPoint(allNewVertices[j]);
				//	}
				//	delete allNewVertices[j];
				//	allNewVertices.erase(allNewVertices.begin() + j);
				//	j--;
				//	s--;
				//}
			}
			j++;
			//}
		}
		i++;
	}

	// check with degenerates vertices, as they may be adjacent to new vertices.
	std::list<Point*>::iterator intermediate;
	bool firstAdjSecond;
	bool secondAdjFirst;
	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); vertex++) {
		if ((*vertex)->isDegenerate()) {
			// check with new points
			for (int i = 0; i < s; i++) {
				//if (abs((*vertex)->degreeOfBounding(antiIdealPoint) - allNewVertices[i]->degreeOfBounding(antiIdealPoint)) <= 1) {
				if ((*vertex)->sizeIntersectionActiveHyperplanes(*allNewVertices[i]) >= H->get_dim()) {
					//if ((*vertex)->isOnEdge(allNewVertices[i])) std::cout << "";
					//if (cas == 0 || cas == 1) {
					allNewVertices[i]->addAdjacentPoint(*vertex);
					(*vertex)->addAdjacentPoint(allNewVertices[i]);
					/*allNewVertices[i]->print();
					std::cout << " is adj to ";
					(*vertex)->print();
					std::cout << "\n";
					hpp = allNewVertices[i]->get_activeHyperplanes();
					std::cout << "\n hpp list 1: \n";
					for (f = hpp->begin(); f != hpp->end(); f++) {
						std::cout << "    -> " << *f << " : ";
						(*f)->print();
					}
					hpp = (*vertex)->get_activeHyperplanes();
					std::cout << "\n hpp list 2: \n";
					for (f = hpp->begin(); f != hpp->end(); f++) {
						std::cout << "    -> " << *f << " : ";
						(*f)->print();
					}*/
					cas = (*vertex)->isOnEdge(allNewVertices[i]);

					//}
				}
				//}
			}

			// check with degenerate vertices also, we may have two degenerate vertices that becomes adjacent, when an hyerplane
			// crosses a facet f through two non-adjacent extreme points of f, making them adjacent.
			intermediate = vertex;
			intermediate++;
			for (vertex2 = intermediate; vertex2 != extrPoints.end(); vertex2++) {
				if ((*vertex2)->isDegenerate()) {
					/*firstAdjSecond = (*vertex)->isAdjacent(*vertex2);
					secondAdjFirst = (*vertex2)->isAdjacent(*vertex);
					if (firstAdjSecond && !secondAdjFirst) {
						(*vertex2)->addAdjacentPoint(*vertex);
					}
					else if (!firstAdjSecond && secondAdjFirst) {
						(*vertex)->addAdjacentPoint(*vertex2);
					}
					else if (!firstAdjSecond && !secondAdjFirst) {
						if ((*vertex)->sizeIntersectionActiveHyperplanes(**vertex2) >= H->get_dim()) {
							(*vertex2)->addAdjacentPoint(*vertex);
							(*vertex)->addAdjacentPoint(*vertex2);
						}
					}*/
					//if (abs((*vertex)->degreeOfBounding(antiIdealPoint) - (*vertex2)->degreeOfBounding(antiIdealPoint)) <= 1) {
					if ((*vertex)->sizeIntersectionActiveHyperplanes(**vertex2) >= H->get_dim()) {
						//if ((*vertex)->isOnEdge(*vertex2)) std::cout << ""; // oof!\n
						//if (cas == 0 || cas == 1) {
						(*vertex2)->addAdjacentPoint(*vertex);
						(*vertex)->addAdjacentPoint(*vertex2);
						/*(*vertex)->print();
						std::cout << " is adj to ";
						(*vertex2)->print();
						std::cout << "\n";
						hpp = (*vertex)->get_activeHyperplanes();
						std::cout << "\n hpp list 1: \n";
						for (f = hpp->begin(); f != hpp->end(); f++) {
							std::cout << "    -> " << *f << " : ";
							(*f)->print();
						}
						hpp = (*vertex2)->get_activeHyperplanes();
						std::cout << " hpp list 2: \n";
						for (f = hpp->begin(); f != hpp->end(); f++) {
							std::cout << "    -> " << *f << " : ";
							(*f)->print();
						}*/
						cas = (*vertex)->isOnEdge(*vertex2);
						//}
					}
					//}
				}
			}
			(*vertex)->becomesNonDegenerate();
		}
	}

	// discard new pts from the correction
	//i = 0;
	//s = allNewVertices.size();
	//while (i < s) {
	//	if (allNewVertices[i]->isDiscarded()) {
	//		// we notify its hyperplanes
	//		listHpp = allNewVertices[i]->get_activeHyperplanes();
	//		for (f = listHpp->begin(); f != listHpp->end(); f++) {
	//			(*f)->removeVertex(allNewVertices[i]);
	//		}
	//		// we notify its adjacent points
	//		adjacentVertex = allNewVertices[i]->get_adjList();
	//		for (vertex3 = adjacentVertex->begin(); vertex3 != adjacentVertex->end(); vertex3++) {
	//			(*vertex3)->removeAdjacentPoint(allNewVertices[i]);
	//		}
	//		delete allNewVertices[i];
	//		allNewVertices.erase(allNewVertices.begin() + i);
	//		i--;
	//		s--;
	//	}
	//	i++;
	//}
	//vertex2 = extrPoints.begin();
	//while (vertex2 != extrPoints.end()) {
	//	vertex = vertex2;
	//	vertex2++;
	//	if ((*vertex)->isDiscarded()) {
	//		// we notify its hyperplanes
	//		listHpp = (*vertex)->get_activeHyperplanes();
	//		for (f = listHpp->begin(); f != listHpp->end(); f++) {
	//			(*f)->removeVertex(*vertex);
	//		}
	//		// we notify its adjacent points
	//		adjacentVertex = (*vertex)->get_adjList();
	//		for (vertex3 = adjacentVertex->begin(); vertex3 != adjacentVertex->end(); vertex3++) {
	//			(*vertex3)->removeAdjacentPoint(*vertex);
	//		}
	//		if ((*vertex)->isCheckpoint()) {
	//			checkPointDestroyed = true;
	//			//std::cout << "that shit just happened putin...\n";
	//		}
	//		delete* vertex;
	//		extrPoints.erase(vertex);
	//	}
	//}

	// add new hpp to list of facets and new pts to list of extreme points
	facets.push_back(H);
	for (int i = 0; i < allNewVertices.size(); i++) {
		allNewVertices[i]->becomesNonDegenerate();
		extrPoints.push_back(allNewVertices[i]);
		/*std::cout << "New point: ";
		allNewVertices[i]->print();
		std::cout << "\n";*/
	}

	/*bool beech = false;
	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); vertex++) {
		if ((*vertex)->isDegenerate()) {
			std::cout << " DEGENERATE VERTEX REMAINING : ";
			(*vertex)->print();
			std::cout << "\n";
			beech = true;
		}
	}

	if (beech) {
		throw std::string("Debug: degenerate vertex status reset");
	}*/

	int pause;
	//if (iteration == 93) {
	//	/*for (f = facets.begin(); f != facets.end(); f++) {
	//		adjacentVertex = (*f)->get_defPts();
	//		std::cout << " HPP is : ";
	//		(*f)->print();
	//		for (vertex = adjacentVertex->begin(); vertex != adjacentVertex->end(); vertex++) {
	//			std::cout << "    -> ";
	//			(*vertex)->print();
	//			std::cout << "\n";
	//		}
	//	}*/
	////	std::list<Hyperplane*>* ahl;
	//	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); vertex++) {

	////		/*ahl = (*vertex)->get_activeHyperplanes();
	////		if (ahl->size() <= 2) {
	////			std::cout << "AH!\n";
	////		}*/

	////		//(*vertex)->purge();
	//		adjacentVertex = (*vertex)->get_adjList();
	//		//if (adjacentVertex->size() >= 6) {
	//			std::cout << "\n\n Point ";
	//			(*vertex)->print();
	//			std::cout << " (" << *vertex << ") with adj list:";
	//			for (vertex2 = adjacentVertex->begin(); vertex2 != adjacentVertex->end(); vertex2++) {
	//				std::cout << "\n      -> ";
	//				(*vertex2)->print();
	//				std::cout << " from " << *vertex2;
	//				/*intermediate = vertex2;
	//				intermediate++;
	//				for (vertex3 = intermediate; vertex3 != adjacentVertex->end(); vertex3++) {
	//					if (*vertex2 == *intermediate) {
	//						std::cout << "tamerlaput\n";
	//					}
	//				}*/
	//			}
	//			//std::cin >> pause;
	//		//}
	////		/*if (adjacentVertex->size() <= 2) {
	////			std::cout << "Here we have ze perber.\n";
	////		}*/

	//	}
	//	//std::cin >> pause;
	//}
}

/*! \brief This function update the representations of the linear relaxation by adding a new hyperplane.
 *
 * OLD VERSION
 * \param H Hyperplane*. A pointer to the new hyperplane added to the representation.
 * \param ptsDiscardedByH Point*. A pointer to the point that is cut out of the lower bound set by the new hyperplane H.
 */
void LinearRelaxation::updatePolyhedron(Hyperplane* H, Point* ptsDiscardedByH) {

	std::list<Point*>::iterator vertex;
	std::list<Point*>::iterator vertex2;
	std::list <Hyperplane*>::iterator f;
	std::list <Hyperplane*>::iterator fCurrent;
	std::list<Hyperplane*>* fln;
	std::list<Point*>* adjacentVertex;

	int nbNewPts = 0;

	facets.push_back(H);
	bool stop = false; // for debugging

	//print();
	//std::cout << " \n\n ---- new it ----\n";
	//if (iteration == 7) { // 17239
	//	H->print();
	////////	std::cout << "\n";
	////////	//print();
	//}
	// search for infeasible and degenerated points
	vertex = extrPoints.begin();
	//std::cout << "\n dominance test : \n";
	while (vertex != extrPoints.end()) {

		//(*vertex)->print();
		//adjacentVertex = (*vertex)->get_adjList();
		//for (vertex2 = adjacentVertex->begin(); vertex2 != adjacentVertex->end(); ++vertex2) {
			//std::cout << "\n     adj vertex : ";
			//(*vertex2)->print();
		//}
		//std::cout << "\n";
		(*vertex)->becomesNonDegenerate();
		if (*vertex == ptsDiscardedByH) {
			(*vertex)->becomesDiscarded();
		}
		else if ((*vertex)->locatedOn(*H)) {
			(*vertex)->becomesDegenerate();
			/*if (iteration == 2817) {
				(*vertex)->print();
				std::cout << " is degenerated\n";
			}*/
		}
		else if (!(*vertex)->isCheckpoint() && (*vertex)->below(*H)) { // !(*vertex)->isCheckpoint() && 
			(*vertex)->becomesDiscarded();
			/*if (iteration == 2817) {
				(*vertex)->print();
				std::cout << " is discarded\n";
			}*/
		}
		else {
			/*if (iteration == 2817) {
				(*vertex)->print();
				std::cout << " is feasible\n";
			}*/
		}
		++vertex;
	}

	// search for new extreme points (comparison with edges, by adjacency lists)
	Point* newPts;
	std::vector<Point*> allNewVertices(0);

	for (vertex = extrPoints.begin(); vertex != extrPoints.end(); ++vertex) { // for each extreme point

		if ((*vertex)->isDiscarded()) { // if it is discarded

			adjacentVertex = (*vertex)->get_adjList(); // we look at its adjacency list

			/*if (iteration == 2817) {
				std::cout << " \n --> Discarded point: ";
				(*vertex)->print();
			}*/

			for (vertex2 = adjacentVertex->begin(); vertex2 != adjacentVertex->end(); ++vertex2) {

				/*if (iteration == 2817) {
					std::cout << "\n     check with ";
					(*vertex2)->print();
				}*/
				if (!(*vertex2)->isDiscarded() && !(*vertex2)->isDegenerate()) { // for each non-discarded adjacent vertex, we compute the new point

					//(*vertex)->print();
					//(*vertex2)->print();
					newPts = new Point((*vertex2)->edgeIntersection(**vertex, *H)); // new pts as intersection of H and edge
					(*vertex2)->replaceAdjVertex((*vertex)->get_adress(), newPts); // update adj vertex for feasible one
					newPts->addAdjacentPoint((*vertex2)->get_adress()); // init adj vertex for new one
					newPts->updateActiveHyperplanes(**vertex, **vertex2, H); // add a new active hyperplanes [ICI] !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					extrPoints.push_back(newPts); //new pts added to list of extreme points
					nbNewPts++;
					if (iteration == 7) {
						std::cout << "\n\n --> New pts " << nbNewPts << ": ";
						newPts->print();
						/*std::cout << " inherited from ";
						(*vertex)->print();
						std::cout << " and ";
						(*vertex2)->print();*/
						std::cout << "\n";
						/*std::cout << "\n Defining hpp for V1:\n";
						fln = (*vertex)->get_activeHyperplanes();
						for (f = fln->begin(); f != fln->end(); f++) {
							(*f)->print();
						}*/
						/*std::cout << "\n Defining hpp for V2:\n";
						fln = (*vertex2)->get_activeHyperplanes();
						for (f = fln->begin(); f != fln->end(); f++) {
							(*f)->print();
						}*/

						//std::cout << *vertex << " and " << *vertex2 << "\n"; // " inherited from " <<
						//if (nbNewPts == 7) throw std::string("stop debug");
					}
					//std::cout << "\nat intersection of " << *vertex << " and " << *vertex2;
					//std::cout << "\n --------------------\n";
					//std::cout << "   ->    New pts: ";
					//newPts->print();
					//std::cout << " inherited from ";
					//(*vertex)->print();
					//std::cout << " and ";
					//(*vertex2)->print();
					//std::cout << "\n";
					//for (int i = 0; i < allNewVertices.size(); i++) {
					//	if (*allNewVertices[i] == *newPts) {
					//		allNewVertices[i]->print();
					//		std::cout << " -> i  = " << i;
					//		newPts->print();
					//		stop = true;
					//		std::cout << "\n";
					//	}
					//}

					allNewVertices.push_back(newPts); // remember the new vertex	
				}
				else if ((*vertex2)->isDegenerate()) {
					(*vertex2)->removeAdjacentPoint(*vertex);
					//std::cout << "  -> Degenerate case";
				}
				//std::cout << "\n";
			}
		}
		else if ((*vertex)->isDegenerate()) {
			(*vertex)->addActiveHyperplane(H);
			H->addVertex(*vertex); //new point to hyperplane H
			allNewVertices.push_back(*vertex); // remember the degenerate vertex
			//(*vertex)->print();
			//std::cout << "\n";
		}
	}

	// Update the adjacency lists for the new vertices + counters hyperplanes
	int s = allNewVertices.size();
	std::list<Hyperplane*>* hpp;
	bool alreadyAdj = false;

	for (int i = 0; i < s; i++) {
		//for (int j = i + 1; j < s; j++) {
		//	alreadyAdj = false;
		//	if (allNewVertices[i]->isDegenerate() && allNewVertices[j]->isDegenerate()) {
		//		adjacentVertex = allNewVertices[i]->get_adjList();
		//		vertex = adjacentVertex->begin();
		//		while (!alreadyAdj && vertex != adjacentVertex->end()) {
		//			if (allNewVertices[j] == *vertex) {
		//				alreadyAdj = true;
		//			}
		//			vertex++;
		//		}
		//	}
		//	if (!alreadyAdj) {
		//		if (allNewVertices[i]->sizeIntersectionActiveHyperplanes(*allNewVertices[j]) == H->get_dim()) { // is this true for p >= 3 ? yes.
		//			allNewVertices[i]->addAdjacentPoint(allNewVertices[j]);
		//			allNewVertices[j]->addAdjacentPoint(allNewVertices[i]);
		//		}
		//	}
		//}
		if (!allNewVertices[i]->isDegenerate()) {
			hpp = allNewVertices[i]->get_activeHyperplanes();
			for (f = hpp->begin(); f != hpp->end(); f++) {
				(*f)->addVertex(allNewVertices[i]);
			}
		}
	}

	//std::cout << " step 1 ... ";

	// delete non-feasible vertices & update hyperplanes
	std::list <Hyperplane*>* listHpp;
	vertex = extrPoints.begin();

	while (vertex != extrPoints.end()) { // [actually not anymore] relies on the fact that the first extreme point is the antiIdeal and is thus always feasible. This simplifies the list management.
		vertex2 = vertex;
		++vertex;
		if ((*vertex2)->isDiscarded()) {
			listHpp = (*vertex2)->get_activeHyperplanes();
			for (f = listHpp->begin(); f != listHpp->end(); f++) {
				(*f)->removeVertex(*vertex2);
			}
			//(*vertex2)->~Point(); // destroy this point
			delete* vertex2;
			extrPoints.erase(vertex2);
		}
	}

	//std::cout << " step 2 ... ";

	// delete redundants hyperplanes
	f = facets.begin();
	while (f != facets.end()) {
		(*f)->checkRedundancy();
		f++;
	}

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

	//std::cout << "\n And a new purge comes in ...\n";
	f = facets.begin();
	while (f != facets.end()) {
		fCurrent = f;
		f++;
		if ((*fCurrent)->isRedundant()) {
			delete* fCurrent;
			facets.erase(fCurrent);
		}
	}

	// ADDED TO SEE IF THAT DEBUG 4OBJ
	int ittt = 0;
	for (int i = 0; i < s; i++) {
		for (int j = i + 1; j < s; j++) {
			alreadyAdj = false;
			if (allNewVertices[i]->isDegenerate() && allNewVertices[j]->isDegenerate()) {
				adjacentVertex = allNewVertices[i]->get_adjList();
				vertex = adjacentVertex->begin();
				while (!alreadyAdj && vertex != adjacentVertex->end()) {
					if (allNewVertices[j] == *vertex) {
						alreadyAdj = true;
					}
					vertex++;
				}
			}
			if (!alreadyAdj) {
				if (allNewVertices[i]->sizeIntersectionActiveHyperplanes(*allNewVertices[j]) == H->get_dim()) { // is this true for p >= 3 ? yes.
					/*ittt++;
					if (iteration == 7) {
						std::cout << "inter " << ittt << ": ";
						allNewVertices[i]->print();
						std::cout << " and ";
						allNewVertices[j]->print();
						std::cout << " are adjecent.\n";
					}*/
					allNewVertices[i]->addAdjacentPoint(allNewVertices[j]);
					allNewVertices[j]->addAdjacentPoint(allNewVertices[i]);
				}
			}
		}
		/*if (!allNewVertices[i]->isDegenerate()) {
			hpp = allNewVertices[i]->get_activeHyperplanes();
			for (f = hpp->begin(); f != hpp->end(); f++) {
				(*f)->addVertex(allNewVertices[i]);
			}
		}*/
	}

}






// add the new hyperplan to the LB set

Hyperplane* h = new Hyperplane(normalVector, ws);
facets.push_back(h);

// update the set of points
Point* newPts = new Point(antiIdealPoint);
newPts->becomesNonDegenerate();
extrPoints.push_back(newPts);
for (int k = 0; k < lp->get_p(); k++) {
	extrPoints.push_back(new Point(antiIdealPoint));
	extrPoints.back()->setObjVector(k, extrPoints.back()->findMissingCoordinate(*h, k));
	extrPoints.back()->addActiveHyperplane(h);
	extrPoints.back()->becomesNonDegenerate();
	h->addVertex(extrPoints.back());
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

// add active artificial hyperplanes to the initial points
int obj = 0;
for (it1 = extrPoints.begin(); it1 != extrPoints.end(); ++it1) {
	for (int k = 1; k <= lp->get_p(); k++) {
		if (obj != k) {
			(*it1)->addActiveHyperplane(boundingBox[k - 1]);
		}
	}
	obj++;
}

// add defining vertices to the artificial hyperplanes
std::list<Point*>::iterator vertex;
for (int k2 = 0; k2 < lp->get_p(); k2++) {
	vertex = extrPoints.begin();
	boundingBox[k2]->addVertex(*vertex); // add the anti-ideal point
	vertex++;
	while (vertex != extrPoints.end()) { // add the other points
		if ((*vertex)->get_objVector(k2) == antiIdealPoint[k2]) {
			boundingBox[k2]->addVertex(*vertex);
		}
		vertex++;
	}
}

//stat.timeInitialization.StopTimer();
