#pragma once
#include <stdlib.h>

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
		lastSplittedIndex = -1;
		depth = 0;
	}
	BranchingDecisions(int n, int p) {
		lb = std::vector<int>(n);
		ub = std::vector<int>(n);
		slub = std::vector<int>(p);
		lastSplittedIndex = -1;
		depth = -1;
	}

	std::vector<int> lb; //!< a vector that describes the lower bounds on the variables in this node. It thus contains some of the previous branching decisions.
	std::vector<int> ub; //!< a vector that describes the upper bounds on the variables in this node. It thus contains some of the previous branching decisions.
	std::vector<int> slub; //!< a vector that describes the search cone of this node in the objective space. It corresponds to the right-hand side of the objective branching constraints.
	int lastSplittedIndex; //!< index of the last splitted index. -1 if there is no such thing (e.g. in root node)
	int depth;
};

//struct LPRelaxStat {
//	LPRelaxStat() : lpSolved(0), feasibilityCheckSolved(0), furthestFeasbilePointSolved(0), dualBensonSolved(0), totalTime(), timeUpdatePolyhedron(), timeFeasibilityCheck(), timeFurthestFeasiblePoint(), timeInitialization(), profiler() {
//	}
//
//	int lpSolved; //!< total number of LP solved
//	int feasibilityCheckSolved; //!< number of calls to FeasibilityCheck model
//	int furthestFeasbilePointSolved; //!< number of calls to FurthestFeasiblePoint model
//	int dualBensonSolved; //!< number of calls to dualBenson model
//	Timer totalTime;
//	Timer timeUpdatePolyhedron; //!< time spent in updating the polyhedron
//	Timer timeFeasibilityCheck; //!< time spent in checking whether extreme points are feasible and computing pre-images
//	Timer timeFurthestFeasiblePoint; //!< time spent in computing the furthest feasible point given an infeasible point
//	Timer timeDualBenson; //!< time spent in solving the dual problem for finding the new hypeplanes
//	Timer timeInitialization; //!< time spent in initializing the polyhedron
//
//	Timer profiler;
//};

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
const int WARMSTARTED_LP_RELAX = 101; // implemented
const int IMPLICIT_LP_RELAX = 102;

// node selection
const int DEPTH_FIRST = 200; // implemented
const int BREADTH_FIRST = 201; // implemented

// variable selection
const int FIRST_INDEX = 300; // implemented
const int MOST_OFTEN_FRACTIONAL = 301; // implemented

// objective branching
const int NO_OBJECTIVE_BRANCHING = 400; // implemented
const int CONE_OBJECTIVE_BRANCHING = 401;
const int FULL_OBJECTIVE_BRANCHING = 402; // implemented but /!\ mistake -> this is actually cone OB !!!

// extreme rays
const double M = 999999989;
const double EPS_LAMBDA = 0.00001;
const double EPS_PROXIMITY = 0.000001;

// debug
const double PRINT_DEBUG = -1;
const double DEBUG_IT = -1; // 3198
const int CORRECTION_WARMSTART = 7;

// timers
const int TIME_OUT = 600;