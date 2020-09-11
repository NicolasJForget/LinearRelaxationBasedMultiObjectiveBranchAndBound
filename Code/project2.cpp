// project2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <list>
#include <vector>
//#include "Point.h"
//#include "Hyperplane.h"
#include "LowerBoundSet.h"
#include "LinearProgram.h"


int main()
{
    std::cout << "Reading file...\n";
    //LinearProgram lp = LinearProgram("C:/Users/au643334/Documents/MOrepo-Forget20/instances/raw/KP/Forget20-KP_10_3_1000-2000_random_3_3.raw");

    //std::cout << "ok on est good la\n";
    //LowerBoundSet LB = LowerBoundSet(lp);
    //LB.print();

    //std::vector<double> normalVec{ 1,0,0 };
    //std::vector<double> pts{ -3,-3,-3 };
    //Hyperplane H = Hyperplane(normalVec, pts, LB.newId());
    //LB.updatePolytope(H);
    //LB.print();

    MathematicalFormulation lp2 = MathematicalFormulation("C:/Users/au643334/Documents/MOrepo-Forget20/instances/raw/Forget20-KP_10_3_1-1000_random_1_1.raw");
    lp2.formate_to_min();

    //lp2.print_objective();
    //lp2.print_constraints();

    // building cplex models

    CplexModel M = CplexModel();
    M.build_dualBenson(&lp2);
    CplexModel M2 = CplexModel();
    M2.build_bestValidPoint(&lp2);

    // hyperplane

    std::vector<double> y = { 12,10,2 };
    //hpp = 
    M.solveDualBenson(y);

    // best valid pts

    
    std::vector<double> s = { 7,6,8 };
    std::vector<double> phat = { 10,10,10 };
    M2.solveBestValidPoint(s, phat);


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
