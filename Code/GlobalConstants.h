#pragma once
#include <stdlib.h>
//#include "timer.hpp"

struct Statistics {
	Statistics() : totalTime(0), timeLBComputation(0), timeDominanceTest(0), timeUpdateUB(0), nbNodes(0), nbFathomedInfeasibility(0), nbFathomedDominance(0), nbFathomedOptimality() {}

	double totalTime;
	double timeLBComputation;
	double timeDominanceTest;
	double timeUpdateUB;

	int nbNodes;
	int nbFathomedInfeasibility;
	int nbFathomedDominance;
	int nbFathomedOptimality;
};

struct Parameters {
	Parameters() : LBset(0), nodeSelection(0), variableSelection(0), objectiveBranching(0), GCD(0) {}

	int LBset;
	int nodeSelection;
	int variableSelection;
	int objectiveBranching;
	std::vector<int> GCD; //!< greatest common divisor for each objective
};

struct BranchingDecisions {
	BranchingDecisions() {
		lb = std::vector<int>(0);
		ub = std::vector<int>(0);
		slub = std::vector<int>(0);
	}
	BranchingDecisions(int n, int p) {
		lb = std::vector<int>(n);
		ub = std::vector<int>(n);
		slub = std::vector<int>(p);
	}

	std::vector<int> lb; //!< a vector that describes the lower bounds on the variables in this node. It thus contains some of the previous branching decisions.
	std::vector<int> ub; //!< a vector that describes the upper bounds on the variables in this node. It thus contains some of the previous branching decisions.
	std::vector<int> slub; //!< a vector that describes the search cone of this node in the objective space. It corresponds to the right-hand side of the objective branching constraints.
};

// status
const int UNSOLVED = 0;
const int SOLVED = 1;
const int INFEASIBLE = 2;
const int OPTIMAL = 3;
const int UNBOUNDED = 4;
const int DOMINATED = 5;

// others
const int IS_LB = 10;
const int IS_UB = 11;

// LB sets
const int LP_RELAX = 100; // implemented
const int WARMSTARTED_LP_RELAX = 101;
const int IMPLICIT_LP_RELAX = 102;

// node selection
const int DEPTH_FIRST = 200; // implemented
const int BREADTH_FIRST = 201; // implemented

// variable selection
const int FIRST_INDEX = 300;
const int MOST_OFTEN_FRACTIONAL = 301;

// objective branching
const int NO_OBJECTIVE_BRANCHING = 400;
const int CONE_OBJECTIVE_BRANCHING = 401;
const int FULL_OBJECTIVE_BRANCHING = 402;
