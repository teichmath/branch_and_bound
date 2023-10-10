# branch_and_bound
Solves the Travelling Salesman Problem for any graph.
Graph: A double, (V, E); V: a set of vertices; E: a set of weights assigned to pairs of vertices in V.
Much of the logic here could be used to solve any binary integer program, in particular to solve 
graph problems in addition to the TSP. One of these is the Max Cut problem, which will be included in an update.

Instructions:
The main function is tsp_driver in tsp_driver.py. The graph argument for this function should be a path to a .txt file with the following format:

n m                 // # of nodes and # of edges
end1 end2 weight    // for edge(0)
...
end1 end2 weight    // for edge(m-1)

Solving Overview:
No self-contained, graph-based combinatorial algorithm exists as a general solution to the TSP. Instead, we formulate the problem as a binary integer program, with the edge set indication vector as the variable of the program. We then explore a tree of subproblems, in which each tree node restricts the domain of the edge set vector by equality to 0 or 1 for some subset of edges. Each of these subproblems is then solved as an LP relaxation, for which the optimum could be non-integral. We also solve special sub-graph problems that provide us with cutting planes, which are inequalities that further restrict our feasible region, thus accelerating the search.
As we explore the tree, we keep track of the known lower and upper bounds on the optimal objective value. The lower bound is determined by taking the minimum value of the subproblem results at all live tree nodes; the upper bound is determined by taking the minimum value of all verified vertex tours, found either heuristically (by some greedy algorithm) or discovered at a node of the tree itself.
A branch of the tree can terminate at a node either because the subproblem at that node is infeasible, or its solution is a vertex tour, or its solution is greater than the known upper bound. The process ends when all tree branches have terminated, or when our upper bound - lower bound < 1.




