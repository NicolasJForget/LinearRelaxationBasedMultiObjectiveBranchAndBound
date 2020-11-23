// project2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <list>
#include <vector>
#include "BranchAndBound.h"
//#include "LowerBoundSet.h"
//#include "LinearProgram.h"
//#include "LB.h"
//#include "Model.h"
//#include "UB.h"


int main()
{
    //std::cout << "Reading file...\n";
    //MathematicalModel lp = MathematicalModel();
    //lp.fill("C:/Users/au643334/source/repos/LinearRelaxationBasedMultiObjectiveBranchAndBound/Code/instances/KP10-5.txt");
    //lp.fill("C:/Users/au643334/Documents/MOrepo-Forget20/instances/raw/Forget20-UFLP_7_3_1-1000_1-100_spheredown_1_1.raw"); // Forget20-UFLP_7_3_1-1000_1-100_spheredown_1_1   Forget20-AP_13_3_1-1000_spheredown_1_1
    //lp.fill("C:/Users/au643334/source/repos/LinearRelaxationBasedMultiObjectiveBranchAndBound/Code/instances/small.txt");
    //lp.formateToMin();


    //std::cout << " Done !\nBuilding Cplex models...";

    //CplexModel db = CplexModel();
    //db.buildDualBenson(&lp);
    //CplexModel bvp = CplexModel();
    //bvp.buildBestValidPoint(&lp);
    //CplexModel feas = CplexModel();
    //feas.buildFeasibility(&lp);
    //CplexModel ws = CplexModel();
    //ws.buildWeightedSumScalarization(&lp);
    //WeightedSumModel ws = WeightedSumModel();
    //ws.build(lp);
    //FeasibilityCheckModel fc = FeasibilityCheckModel();
    //fc.build(lp);
    //DualBensonModel db = DualBensonModel();
    //db.build(lp);
    //FurthestFeasiblePointModel ffp = FurthestFeasiblePointModel();
    //ffp.build(lp);

    //std::cout << " Done !\nInitializing the lower bound set...";

    //LowerBoundSet* LB = new LinearRelaxation(&lp, &ws, &fc, &db, &ffp);
    //LinearRelaxation* LP = (LinearRelaxation*) LB; // need to convert LB to LinearRelaxation to call the copy constructor
    //LowerBoundSet* LB2 = new LinearRelaxation( *(LinearRelaxation*)LB );

    //LB->print();
    //LB2->print();

    //LowerBoundSet LB = LowerBoundSet(lp,&ws,&feas,&db,&bvp);

    //std::cout << " Done !\nSolving ...";

    //LB->compute();
    //LB2->compute();

    //std::cout << "Done !\n";

    //LB.print();
    //LB->print();
    //LB2->print();

    //UpperBoundSet U = UpperBoundSet(lp);
    //LB->gatherIntegerSolutions(U);
    //U.print();
    //U.printLub();

    std::string inst = "C:/Users/au643334/source/repos/LinearRelaxationBasedMultiObjectiveBranchAndBound/Code/instances/KP10-3.txt";
    //std::string inst = "C:/Users/au643334/Documents/MOrepo-Forget20/instances/raw/Forget20-KP_15_3_1-1000_spheredown_3_1.raw";
    BranchAndBound B = BranchAndBound(inst);

    // run 1

    //B.run(LP_RELAX, DEPTH_FIRST, FIRST_INDEX, NO_OBJECTIVE_BRANCHING);
    //B.printYN();
    //B.printStatistics();

    // run 2

    B.run(LP_RELAX, BREADTH_FIRST, FIRST_INDEX, NO_OBJECTIVE_BRANCHING);
    B.printStatistics();

    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
