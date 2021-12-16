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
* Parameter 4: objective branching
  - 400 to use no objective branching (noOB)
  - 401 to use cone objective branching (coneOB)
  - 402 to use full objective branching (fullOB)
  - 403 to use objective branching with an upper bound on the number of sub-problems that can be created (limitedOB)
* Parameter 5: splitting value rule v, when branching on a variable by creating sub-problems x_i <= v and x_i >= v + 1
  - 500 to use the median value as a splitting rule
  - 501 to use the most often fractional value as a splitting rule
  - 502 to select the value randomly
* Parameter 6: limit on the sub-problems created. If parameter 4 is 400, 401, or 402, the value of this parameter does not matter and will be ignored. If parameter 4 is 403 (limitedOB), parameter 6 is the maximum number of sub-problems that can be created at a given node.
* Parameter 7: time limit, expressed in seconds.
  
Note that for Parameter 4, all rules will perform similarly in case binary variables only are present in the instance solved.
To run the branch-and-bound from "Branch-and-bound and objective branching with three or more objective" (N. Forget, S.L. Gadegaard, K.Klamroth, L.R. Nielsen, A.Przybylski) on instance Forget21-UFLP_6_3_1-1000_1-100_spheredown_1_3.txt (located in Code\instances), one should run the code with the following parameters (1 hour time limit):
* With noOB: instance\Forget21-UFLP_6_3_1-1000_1-100_spheredown_1_3.txt 101 201 400 501 3600 9999
* With coneOB: instance\Forget21-UFLP_6_3_1-1000_1-100_spheredown_1_3.txt 101 201 401 501 3600 9999
* With fullOB: instance\Forget21-UFLP_6_3_1-1000_1-100_spheredown_1_3.txt 101 201 402 501 3600 9999
* With limitedOB (5 sub-problems max): instance\Forget21-UFLP_6_3_1-1000_1-100_spheredown_1_3.txt 101 201 403 501 3600 5

All test instance can be found [here](https://github.com/MCDMSociety/MOrepo-Forget21) for PPP and UFLP, and [here](https://github.com/MCDMSociety/MOrepo-Kirlik14) for KP and ILP.
