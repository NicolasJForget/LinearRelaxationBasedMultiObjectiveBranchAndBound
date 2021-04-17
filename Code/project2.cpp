// project2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <list>
#include <vector>
#include <filesystem>
#include "BranchAndBound.h"
//#include "LowerBoundSet.h"
//#include "LinearProgram.h"
//#include "LB.h"
//#include "Model.h"
//#include "UB.h"

void expe() {

    // getting the list of instance files
    namespace fs = std::filesystem;
    std::vector<std::string> filenames;
    std::string path = "C:/Users/au643334/source/repos/LinearRelaxationBasedMultiObjectiveBranchAndBound/Code/instances/expeData/newInstances/";
    const fs::directory_iterator end{};
    for (fs::directory_iterator iter{ path }; iter != end; ++iter) { 
        filenames.push_back(iter->path().string());
    }

    // parameters to test
    std::vector<int> LB = { LP_RELAX , WARMSTARTED_LP_RELAX };
    std::vector<int> nodeSel = { BREADTH_FIRST, DEPTH_FIRST }; // DEPTH_FIRST , 
    std::vector<int> objectiveBranching = { CONE_OBJECTIVE_BRANCHING, FULL_OBJECTIVE_BRANCHING, NO_OBJECTIVE_BRANCHING }; // NO_OBJECTIVE_BRANCHING , FULL_OBJECTIVE_BRANCHING , 

    // run expe on each file
    for (auto instance : filenames) {
        
        // building the BB for the instance
        std::cout << "\n\n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
        std::cout << " Currently running : " << instance;
        BranchAndBound B = BranchAndBound(instance);

        // run with the various test parameters
        for (auto configLB : LB) {
            for (auto configNodeSel : nodeSel) {
                for (auto configOB : objectiveBranching) {
                    B.run(configLB, configNodeSel, MOST_OFTEN_FRACTIONAL, configOB);
                    B.printStatistics();
                    B.writeStatistics();
                }
            }
        }

    }
}


int main()
{
    
    std::cout.precision(4);
    expe();

    //std::string inst = "C:/Users/au643334/source/repos/LinearRelaxationBasedMultiObjectiveBranchAndBound/Code/instances/Kirlik14-ILP_p-3_n-10_m-5_ins-3.dat";
    //std::string inst = "C:/Users/au643334/source/repos/LinearRelaxationBasedMultiObjectiveBranchAndBound/Code/instances/debug_3obj_int.dat";
    //std::string inst = "C:/Users/au643334/source/repos/LinearRelaxationBasedMultiObjectiveBranchAndBound/Code/instances/expeData/waiting/Forget21-UFLP_4_6_1-1000_1-100_spheredown_1_1.txt";
    //std::string inst = "C:/Users/au643334/source/repos/LinearRelaxationBasedMultiObjectiveBranchAndBound/Code/instances/UFLP6-3.txt";
    //std::string inst = "C:/Users/au643334/source/repos/LinearRelaxationBasedMultiObjectiveBranchAndBound/Code/instances/checkLub4obj.txt";
    //std::string inst = "C:/Users/au643334/source/repos/LinearRelaxationBasedMultiObjectiveBranchAndBound/Code/instances/KP10-6_int.txt"; // KP10-4Debug2
    //std::string inst = "C:/Users/au643334/source/repos/LinearRelaxationBasedMultiObjectiveBranchAndBound/Code/instances/counterProofHyperplanes.txt";
    //std::string inst = "C:/Users/au643334/Documents/MOrepo-Forget20/instances/raw/Forget20-UFLP_6_3_1-1000_1-100_spheredown_1_1.raw"; //Forget20-KP_20_3_1-1000_spheredown_3_1.raw
    //std::string inst = "C:/Users/au643334/source/repos/LinearRelaxationBasedMultiObjectiveBranchAndBound/Code/instances/KP20-debug.raw";
    //BranchAndBound B = BranchAndBound(inst);
    

    // run 1

    //B.run(WARMSTARTED_LP_RELAX, BREADTH_FIRST, MOST_OFTEN_FRACTIONAL, NO_OBJECTIVE_BRANCHING); //FULL_OBJECTIVE_BRANCHING
    //B.printYN();
    //B.printStatistics();

    // run 2

    //B.run(LP_RELAX, DEPTH_FIRST, MOST_OFTEN_FRACTIONAL, CONE_OBJECTIVE_BRANCHING);
    //B.printYN();
    //B.printStatistics();
    //B.writeStatistics();

    // run 3

    //B.run(WARMSTARTED_LP_RELAX, DEPTH_FIRST, MOST_OFTEN_FRACTIONAL, NO_OBJECTIVE_BRANCHING);
    //B.printYN();
    //B.printStatistics();
    //B.writeStatistics();

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
