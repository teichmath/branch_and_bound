# branch_and_bound
Solvers for NP-Hard Graph Problems
2018-2023 A.Teich

*** Instructions ***
For these solvers, you need a graph argument, which should be a path to a .txt file with the following format:

n m                 // # of nodes and # of edges
end1 end2 weight    // for edge(0)
...
end1 end2 weight    // for edge(m-1)

For the Traveling Salesman Problem:
example:
    from solvers import tsp_solver
    graph_name = "TSP/gr21.txt"
    tsp_solver(graph_name)

For the Max Cut problem:
example:
    from solvers import max_cut_solver
    max_cut_solver(graph_name)

For a different branch-and-cut problem (see Discussion, below):
    from bnb_driver import bnb_driver
    from my_own_module import some_function_1, some_function_2, some_function_3, some_function_4, some_function_5, some_function_6 
    search_config = { 
        'problem name': 'Some Graph Problem',
        'graph name': graph_name,
        'optimization type': 'minimize',                #or 'maximize'
        'valid edge set type': 'buncha nice edges',     #whatever you'd call a valid edge set in your problem
        'edge set validator': some_function_1,              #function that checks edge set for validity 
        'greedy heuristic': some_function_2,                #function that gets a quick valid edge set
        'improvement heuristic': some_function_3,           #function that can improve on a given valid edge set
        'initial boundary finder': some_function_4,         #function that provides a cutting plane, even before the first subproblem
        'cutting plane finder': some_function_5,            #function that finds cutting planes during the rounds of branch and cut
        'vertex set description': "Optimal tour vertex sequence: ",     #whatever you'd call a vertex representation of your valid solution
        'vertex set getter': some_function_6                #function that finds a vertex representation from your edge representation of a valid edge set
    }
    solution = bnb_driver(graph_name, search_config, verbosity=1)


*** Discussion ***
Graph: A double, (V, E); V: a set of vertices, |V|=n; E: a set of weights assigned to pairs of vertices in V, |E|=m.

With our method bnb_driver, we have an implementation of the branch-and-cut strategy (see citation below) that can be adapted to create a solver for the Traveling Salesman Problem, the Max Cut problem, or any other graph problem satisying these criteria:
    1.)  The solution is in the form of an edge set. (It can be expressed as an indicator vector in R_m.)
    2.)  There exists a set A of edge sets such that the solution is in the set and is either the maximally or minimally weighted of the set.
    3.)  In polynomial time, we can determine whether a given edge set belongs to A.
    4.)  Our set A excludes many of our graph's possible edge sets, and we can express this in the form of separating hyperplanes.

In the case of the TSP, our set A is all vertex tours, our solution is the minimally weighted vertex tour, and we can find separating hyperplanes by noting that any edge set that is a vertex tour must have a minimum graph cut of 2 (thus, some combination of variables in our solution is >= 2).

To make a solver for a problem meeting the above criteria, you need to procure functions to perform the tasks listed above in the Instructions section. At the very least, you need an edge set validator and a cutting plane finder. A greedy heuristic is not strictly necessary, but not providing an initial bound on the solution could greatly extend the running time of the algorithm. The other components are nice-to-haves, for which you could provide pass functions. 

The idea of branch-and-bound, or branch-and-cut (largely interchangeable terms), is that we are dealing with problems for which no self-contained, graph-based combinatorial algorithm exists as a general solution. Instead, we formulate the problem as a binary integer program, with the edge set indicator vector as the variable of the program. We then explore a tree of subproblems, in which each tree node restricts the domain of the edge set vector by equality to 0 or 1 for some subset of edges. Each of these subproblems is then solved as an LP relaxation, for which the optimum could be non-integral. We also solve special sub-graph problems that provide us with cutting planes, which are inequalities that further restrict our feasible region, thus accelerating the search.

As we explore the tree, we keep track of the known lower and upper bounds on the optimal objective value. The lower bound is determined by taking the minimum value of the subproblem results at all live tree nodes; the upper bound is determined by taking the minimum value of all verified vertex tours, found either heuristically (by some greedy algorithm) or discovered at a node of the tree itself.

A branch of the tree can terminate at a node either because the subproblem at that node is infeasible, or its solution is a vertex tour, or its solution is greater than the known upper bound. The process ends when all tree branches have terminated, or when our upper bound - lower bound < 1.

Citation:
 M. Padberg and G. Rinaldi, "A branch-and-cut algorithm for the resolution of large-scale symmetric traveling salesman problems," SIAM Review 33 (1991) 60-100.
 




