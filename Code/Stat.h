#pragma once
#include "timer.hpp"

struct Statistics {
	Statistics() : solved(0), totalTime(), timeLBComputation(), timeDominanceTest(), timeUpdateUB(), timeComputeOB(), timeVarSel(), timeNodeSel(), nbNodes(0), nbFathomedInfeasibility(0), nbFathomedDominance(0), nbFathomedOptimality(0), minDepth(10000), maxDepth(0), avgDepth(0), lpSolved(0), timeUpdatePolyhedron(), timeFeasibilityCheck(), timeInitialization(), profiler(0) {}

	// performance
	int solved; //!< X 1 if the instance is solved to optimality, 0 otherwise
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

	double profiler;
};