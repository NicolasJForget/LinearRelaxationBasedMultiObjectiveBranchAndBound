# LinearRelaxationBasedMultiObjectiveBranchAndBound

The C++ code for the multi-objective branch-and-bound algorithm is contained in folder Code. The main function is located in project2.cpp. The code uses CPLEX 12.10 to solve linear programs. It can solve linear programs with two or more objective functions, and with integer variables only.

To solve an instance, one can run the main function with the following parameters:

* Parameter 1: path to the data file of the instance.
* Parameter 2: lower bound set computation
  - 100 to compute the linear relaxation from scratch at each node.
  - 101 to warmstart the computation of the linear relaxation using the lower bound set from the father node.
* Parameter 3: node selection rule
  - 200 to use a depth-first strategy
  - 201 to use a breadth-first strategy
* Parameter 4: splitting value rule v, when branching on a variable by creating sub-problems x_i <= v and x_i >= v + 1
  - 500 to use the median value as a splitting rule
  - 501 to use the most often fractional value as a splitting rule
  - 502 to select the value randomly
* Parameter 5: time limit, expressed in seconds.
  
Note that for Parameter 4, all rules will perform similarly in case binary variables only are present in the instance solved.
To run the best configuration in average from "Warm-starting lower bound set computations for branch-and-bound algorithmsfor multi objective integer linear programs" (N. Forget, S.L. Gadegaard, L.R. Nielsen) on instance Forget21-UFLP_6_3_1-1000_1-100_spheredown_1_3.txt (located in Code\instances), one should run the code with the following parameters:
* With one hour time limit: instance\Forget21-UFLP_6_3_1-1000_1-100_spheredown_1_3.txt 101 201 501 3600
* With four hours time limit: instance\Forget21-UFLP_6_3_1-1000_1-100_spheredown_1_3.txt 101 201 501 14400

All test instance can be found [here](https://github.com/MCDMSociety/MOrepo-Forget21) for PPP and UFLP, and [here](https://github.com/MCDMSociety/MOrepo-Kirlik14) for KP and ILP.
