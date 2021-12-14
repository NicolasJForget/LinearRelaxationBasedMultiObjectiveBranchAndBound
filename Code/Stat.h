#pragma once
#include "timer.hpp"

struct Statistics {
	Statistics() : solved(0), totalTime(), timeLBComputation(), timeDominanceTest(), timeUpdateUB(), timeComputeOB(), timeVarSel(), timeNodeSel(), nbNodes(0), nbFathomedInfeasibility(0), nbFathomedDominance(0), nbFathomedOptimality(0), 
		minDepth(10000), maxDepth(0), avgDepth(0), lpSolved(0), timeUpdatePolyhedron(), timeFeasibilityCheck(), timeInitialization(), timeCopyLB(), avgFacets(0), avgVertices(0), maxFacets(0), maxVertices(0), minFacets(0), minVertices(0), nbNodesAtStage(0), 
		nbInfeasibleNodes(0), timeUpdatePoly(0), nbFeasVtx(0), nbNewFacets(0), vtxOnly(0), facetsNoRay(0), newFacetsNoRay(0), avgIntegerVtx(0), cpuCurrent(0.0), nbFeasVtxCurrent(0), nbNewFacetsCurrent(0), vtxOnlyCurrent(0), facetsNoRayCurrent(0), 
		newFacetsNoRayCurrent(0), avgIntegerVtxCurrent(0), cpuDepth(0), lpProbing(0), varFixedProbing(0), nbNodesClosedProbing(0),
		nbProbingInfeas(0), nbProbingOpti(0), nbFixingLp(0), cpuProbing(), cpuProbingCplex(), cpuProbingManual(), avgDepthOB(0), minDepthOB(100000), maxDepthOB(0), 
		avgSubPbOB(0), minSubPbOB(100000), maxSubPbOB(0), nbOccurencesOB(0), avgNbLubs(0), maxNbLubs(0), nbDominanceTestLubHpp(0.0), nodesOBPerDepth(0), avgSubPbOBPerDepth(0), minSubPbOBPerDepth(0), maxSubPbOBPerDepth(0) , profiler(0), debug() {} // , cpuSlubsOB(), cpuAdditionalDomiTestOB()

	// performance
	int solved; //!< 1 if the instance is solved to optimality, 0 otherwise
	Timer totalTime; //!< total cpu time
	Timer timeLBComputation; //!< cpu time spent in computing LB sets
	Timer timeDominanceTest; //!< cpu time spent in applying dominance tests
	Timer timeUpdateUB; //!< cpu time spent in updating the UB set
	Timer timeComputeOB; //!< cpu time spent in computing objective branching
	Timer timeVarSel; //!< cpu time spent in selecting a variable (and a value) to split a node
	Timer timeNodeSel; //!< cpu time spent in selecting a node to explore

	// BB tree stat
	int nbNodes; //!< number of nodes explored
	int nbFathomedInfeasibility; //!< number of node fathomed by infeasibility
	int nbFathomedDominance; //!< number of node fathomed by dominance
	int nbFathomedOptimality; //!< number of node fathomed by optimality
	int minDepth; //!< depth of the least deep leaf node in the tree
	int maxDepth; //!< average depth of the leaf nodes
	double avgDepth; //!< depth of the deepest lead node in the tree

	// LP part
	int lpSolved; //!< total number of LP solved
	Timer timeUpdatePolyhedron; //!< time spent in updating the polyhedron
	Timer timeFeasibilityCheck; //!< time spent in checking whether extreme points are feasible and computing pre-images
	Timer timeInitialization; //!< time spent in initializing the polyhedron
	Timer timeCopyLB; //!< time spent in copying lower bound sets (warmstart only)

	// Per depth stat
	std::vector<double> avgFacets; //!< avg number of facet at each depth
	std::vector<double> avgVertices; //!< avg number of vertices at each depth
	std::vector<double> maxFacets; //!< max number of facets at each depth
	std::vector<double> maxVertices; //!< max number of vertices at each depth
	std::vector<double> minFacets; //!< min number of facets at each depth
	std::vector<double> minVertices; //!< min number of vertices at each depth
	std::vector<int> nbNodesAtStage; //!< number of nodes existing at each depth
	std::vector<int> nbInfeasibleNodes; //!< number of infeasible nodes at each depth
	std::vector<double> timeUpdatePoly; //!< cpu time spent in updating polyhedra at each depth
	std::vector<int> nbFeasVtx; //!< number of feasible vertices from the father node
	std::vector<int> nbNewFacets; //!< number of facets generated in the current node
	std::vector<int> vtxOnly; //!< number of vertices, excluding extreme rays
	std::vector<int> facetsNoRay; //!< number of facets, excluding wnd facets
	std::vector<int> newFacetsNoRay; //!< number of new facets facets, excluding wnd facets
	std::vector<int> avgIntegerVtx; //!< avg number of integer vertices
	double cpuCurrent; //!< just an intermediate for timeUpdatePoly
	int nbFeasVtxCurrent; //!< jut an intermediate for nbFeasVtx
	int nbNewFacetsCurrent; //!< just an intermediate for nbNewFacetsCurrent
	int vtxOnlyCurrent; //!< just an intermediate for vtxOnly
	int facetsNoRayCurrent; //!< just an intermediate for vtxOnly
	int newFacetsNoRayCurrent; //!< just an intermediate for vtxOnly
	int avgIntegerVtxCurrent;
	std::vector<double> cpuDepth; //!< total cpu time spent at each depth

	// Probing
	int lpProbing;
	int varFixedProbing;
	int nbNodesClosedProbing;
	int nbProbingInfeas;
	int nbProbingOpti;
	int nbFixingLp;
	Timer cpuProbing;
	Timer cpuProbingCplex;
	Timer cpuProbingManual;

	// OB -> refer to when 2 or more subPb are created in objective space
	double avgDepthOB; //!< average depth where OB is applied
	int minDepthOB; //!< min depth where OB is applied
	int maxDepthOB; //!< max depth where OB is applied
	double avgSubPbOB; //!< average number of sub-problem created
	int minSubPbOB; //!< min number of sub-problem created (valid only if OB is applied at least once, i.e. nbOccurencesOB >= 1)
	int maxSubPbOB; //!< max number of sub-problems created
	int nbOccurencesOB; //!< number of nodes where OB is applied
	//Timer cpuSlubsOB; //!< total cpu time spent in computing SLUBs
	//Timer cpuAdditionalDomiTestOB; //!< total cpu time spent in additional dominance test (when non-dominance of the node is already known, but we need the status of each LUB)
	double avgNbLubs; //!< average number of local upper bounds per node
	int maxNbLubs; //!< maximum number of local upper bounds
	double nbDominanceTestLubHpp; //!< nb of test of dominance between hyperplanes and lubs
	std::vector<int> nodesOBPerDepth;
	std::vector<double> avgSubPbOBPerDepth;
	std::vector<int> minSubPbOBPerDepth;
	std::vector<int> maxSubPbOBPerDepth;

	double profiler;
	Timer debug;
};



struct StatProbing {
	StatProbing() : lpUsed(0), varFixed(0), nodesConcerned(0), ruleInfeas(0), ruleOpti(0), nbNodes(0), cpuProbing(), cpuTotal() {}

	int lpUsed;
	int varFixed;
	int nodesConcerned;
	int ruleInfeas;
	int ruleOpti;

	int nbNodes;
	Timer cpuProbing;
	Timer cpuTotal;
};