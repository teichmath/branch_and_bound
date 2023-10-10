"""Branch and cut logic, plus helper functions."""
from collections import deque
import copy
import numpy as np
import random
import time
from heuristics import my_2_opt, makeSortedTourVerticesFromSolutionVector
from problem_def import get_single_data, Solution, TSPGraph, TSPLP
from subgraph_ops import check_if_subgraph_is_tour, find_random_subgraph_cut, global_min_cut
from lp_choice import LPChoiceManager

def tsp_driver(graph_name, verbose=False):
    """
    main function of the TSP solver. Initializes the problem and performs logical operations of the branch-and-cut strategy.
    Inputs:
        graph_name (string): path to text file containing graph description. See README.md.
        verbose (boolean): When True, output is provided to show logic. When False, output is limited to lower bound and upper bound updates.
    Outputs:
        tour (list):  list of vertex indices, sorted in order they are visited, for the optimal vertex tour.
        current_UB (double): total weight of the optimal vertex tour.    
    """
    start = time.time()
    vprint = vPrint(verbose)

    #get input
    input_data = get_single_data(graph_name)

    tsp_graph = TSPGraph(input_data)
    #g, g_inv, g_inv_sb, optimal, n, m = make_graph(input_data)

    #initialize the lp
    tsp_lp = TSPLP(tsp_graph)
    best_solution = Solution(tsp_graph, name="initial")
    _ = best_solution.try_improvement(limit=5e6)

    initial_upper_bound = best_solution.get_objective_value()
    
    print("got initial upper bound: "+str(initial_upper_bound))

    #this queue stores the equality constrained values of each variable in each subproblem. 
    subprob_queue = deque([[]])  #[0],[1],[0,0],[0,1], etc.
    cutting_plane_binder = []
    global_min_cuts = []
    current_UB = initial_upper_bound
    current_LB = -np.inf
    visited_tree_node_list = []
    live_tree_node_list = []
    pruned_tree_node_list = []
    solution_corral = [copy.deepcopy(best_solution)]
    visited_best_x = False
    choice_manager = LPChoiceManager(["favor LB", "random"], [0.5, 0.02], 100)
    tour = None

    while True:
        give_update_nonverbose = False
        #strategy management
        choice_manager.rotate_strategy()
        if choice_manager.new_cycle():
            best_sol = update_solution_corral_and_find_best(solution_corral)
            if best_sol.get_objective_value() < current_UB:
                best_solution = copy.deepcopy(best_sol)
                current_UB = best_sol.get_objective_value()
                prune_from_new_UB(current_UB, live_tree_node_list, pruned_tree_node_list, vprint)
                visited_best_x = False
                give_update_nonverbose = True
            solution_corral.append(Solution(tsp_graph, name=str(len(visited_tree_node_list))))
        
        #Choose an LP
        sub_lp = choice_manager.get_current_strategy()(subprob_queue, live_tree_node_list)
        visited_tree_node_list.append(sub_lp)

        vprint("starting subproblem:")
        vprint(sub_lp)

        #check if this tree node is pruned, then proceed
        is_pruned = False
        for pruned_node in pruned_tree_node_list:
            if len(pruned_node) <= len(sub_lp):
                if all(x1 == x2 for x1, x2 in zip(pruned_node, sub_lp[:len(pruned_node)])):
                    is_pruned = True
                    vprint("This tree node is pruned, continuing.")
                    break
        if not is_pruned:
            vprint("initializing lp.", end= '  ')
            tsp_lp.initialize_tsp_lp()
            define_domain(tsp_lp, sub_lp, cutting_plane_binder, global_min_cuts, vprint)
            
            vprint("Solving LP...   ", end='  ')
            objective_value, is_feasible = tsp_lp.solve_lp()
            this_tree_node = {'obj_value':objective_value, 'x_lead':sub_lp}
            if not is_feasible:
                vprint("Not feasible. Pruning and continuing.")
                pruned_tree_node_list.append(this_tree_node)
            else:
                vprint("Feasible with objective value:  "+str(objective_value))
                if objective_value > current_UB:
                    vprint("This is greater than the current upper bound. Pruning and continuing.")
                    pruned_tree_node_list.append(this_tree_node)
                else:
                    vprint("This is less than or equal to the current upper bound.")
                    x_solution = tsp_lp.get_existing_solution()
                    #add to live nodes
                    live_tree_node_list.append(this_tree_node)
                    #is the solution integral
                    is_integral = False
                    is_tour = False
                    vprint("Is the solution integral ?    ", end='')
                    if sum([abs(item - round(item)) for item in x_solution]) < 1e-13:
                        is_integral = True
                        vprint("Yes")
                    else:
                        vprint("No")
                    #check if it's a tour
                    vprint("Is the solution a tour?    ", end='')
                    is_tour = check_if_subgraph_is_tour(tsp_graph, x_solution, is_integral)
                    if is_tour:
                        if is_integral:
                            vprint("Yes, the solution is a vertex tour!")
                            vprint("trying for an improvement with 2-opt... ")
                            integral_x_solution, integral_objective_value, _ = my_2_opt(tsp_graph.graph_inv, x_solution, objective_value, limit=1e4)
                            needs_update, visited_best_x_standin = have_we_visited_the_best_x_solution_as_a_tree_node(x_solution, \
                                                                    best_solution, objective_value, current_UB, integral_objective_value)
                            if needs_update:
                                 visited_best_x = visited_best_x_standin
                        else:
                            vprint("The solution is non-integral, but it provides a vertex tour!")
                            integral_x_solution = np.heaviside(x_solution-1e-13, [1]*len(x_solution))
                            integral_objective_value = tsp_lp.get_objective_value_from_vector(integral_x_solution)
                            if integral_objective_value < current_UB:
                                visited_best_x = False
                        this_sol = Solution(tsp_graph, integral_objective_value, integral_x_solution, name=str(len(visited_tree_node_list))+"B")
                        if integral_objective_value < current_UB:
                            vprint("We have a new UB. Updating x_best and current_UB.")
                            best_solution = copy.deepcopy(this_sol)
                            current_UB = this_sol.get_objective_value()
                            prune_from_new_UB(current_UB, live_tree_node_list, pruned_tree_node_list, vprint)
                            solution_corral.append(this_sol)
                            give_update_nonverbose = True
                    else:
                        vprint("  No")
                    if not (is_integral and is_tour):
                        vprint("Solution is not itself a vertex tour, so we will perform these 3 steps: ")
                        #look for cuts (randomly).
                        vprint("1.   Looking for cutting planes by walking...", end='')
                        cutting_plane_dict = find_random_subgraph_cut(tsp_graph, x_solution, stop_crit_new_planes=20, stop_crit_iters=600)
                        vprint(" done.  "+str(len(cutting_plane_dict['start']))+" were found.")
                        if len(cutting_plane_dict['start']) > 0:
                            cutting_plane_binder.append(cutting_plane_dict)
                        #else:
                        #    print("We're trying a global min cut.")
                        #    coeff_row, current_cut_value = global_min_cut(tsp_graph, x_solution)
                        #    print("value:  "+str(current_cut_value))
                        #    global_min_cuts.append(coeff_row)
                        
                        #update LB
                        vprint("2.   Possibly updating LB...  current is "+str(current_LB)+", new is ", end='')
                        live_node_min = min([item['obj_value'] for item in live_tree_node_list])
                        if live_node_min > current_LB:
                            give_update_nonverbose = True
                            current_LB = live_node_min
                        vprint(current_LB)
                        #branch
                        vprint("3.   Branching.")
                        subprob_queue.append(sub_lp + [0])
                        subprob_queue.append(sub_lp + [1])
        if len(subprob_queue)==0:
            #Done.
            reason_to_quit = "subprob queue is empty."
            tour = final_message(best_solution.get_x_vector(), current_UB, reason_to_quit, visited_best_x, tsp_graph)
            break
        if abs(current_LB-current_UB) < 1-5e-2:
            #Done.
            reason_to_quit = "optimality gap is zero."
            tour = final_message(best_solution.get_x_vector(), current_UB, reason_to_quit, True, tsp_graph)
            break
        #if twin has been visited, pop the pop
        if len(sub_lp) > 0:
            this_tree_node_twin = sub_lp[:-1] + [1 - sub_lp[-1]] 
            if this_tree_node_twin in visited_tree_node_list:
                current_x_leads = [item['x_lead'] for item in live_tree_node_list]
                if sub_lp[:-1] in current_x_leads:
                    current_x_leads.remove(sub_lp[:-1])
                    live_tree_node_list = [item for item in live_tree_node_list if item['x_lead'] in current_x_leads]        
        vprint("***********************************************")
        vprint("LB: "+str(current_LB)+"   UB: "+str(current_UB))
        vprint("subprob_queue:     "+str(len(subprob_queue)))
        vprint("visited nodes:  "+str(len(visited_tree_node_list))+"    live nodes:  "+str(len(live_tree_node_list))+"    pruned nodes:  "+str(len(pruned_tree_node_list)))
        if give_update_nonverbose:
            choice_manager.turn_back_clock(15)
            if not verbose:
                print("LB: "+str(current_LB)+"   UB: "+str(current_UB)+"   visited nodes:  "+str(len(visited_tree_node_list))+"   subprob_queue: "+str(len(subprob_queue)))
        else:
            choice_manager.advance_clock(1)
    return tour, current_UB

def define_domain(tsp_lp, sub_lp, cutting_plane_binder, global_min_cuts, vprint):
    """
    Policy for applying constraints to our subproblems.
    Inputs:
        tsp_lp (TSPLP): Our object defining the state and behavior of this subproblem.
        sub_lp (list): item at index i is either 0 or 1, giving the rhs of an equality constraint for edge x_i.
        cutting_plane_binder (list): each item is a dict of information for a constraint to be applied to the feasible region.
        global_min_cuts (list):  each item is a list of coefficients for a constraint, an edge set found to be the minimum cut for some sub-LP solution.
        vprint (vPrint): wrapper for output statements.
    """
    vprint("applying tree node constraints.",  end= '  ')
    tsp_lp.apply_tree_node_constraints(sub_lp)
    #add your degree 2 constraints.
    vprint("applying degree 2 constraints.")
    tsp_lp.apply_degree_2_constraints()
    #add your global min cut constraints.
    if len(global_min_cuts)>0:
        tsp_lp.apply_global_min_cut_constraints(global_min_cuts)
    #add up to certain number of additional cutting plane constraints.
    if len(cutting_plane_binder) > 0:
        vprint("cutting plane binder has "+str(len(cutting_plane_binder))+" groups")
        constraint_groups = random.sample(range(len(cutting_plane_binder)), min(20, len(cutting_plane_binder)))
        for group_idx in constraint_groups:
            vprint("applying group_idx "+str(group_idx)+" cutting plane constraints.")
            tsp_lp.apply_cutting_plane_constraints(cutting_plane_binder[group_idx])

def prune_from_new_UB(current_UB, live_tree_node_list, pruned_tree_node_list, vprint):
    """ 
    Look at all live tree node LP values; add tree nodes to prune list if > upper bound.
    Inputs:
        current_UB (double):  upper bound for the optimal tour weight (i.e., the lowest weight found for a vertex tour, so far).
        live_tree_node_list (list):  each element represents a node of the subproblem tree that has children with unknown sub-LP solutions.
        pruned_tree_node_list (list):  each element represents a node of the subproblem tree that terminates its branch.
        vprint (vPrint): wrapper for output statements.  
    """
    pop_these_live_nodes = []
    for i_tnd, tree_node_dict in enumerate(live_tree_node_list):
        if tree_node_dict['obj_value'] >= current_UB:
            vprint("Pruning tree node for being >= current_UB:", end='')
            vprint(tree_node_dict)
            pop_these_live_nodes.append(i_tnd)
            pruned_tree_node_list.append(tree_node_dict)
    for i_tnd in sorted(pop_these_live_nodes, reverse=True):
        del live_tree_node_list[i_tnd]

def have_we_visited_the_best_x_solution_as_a_tree_node(x_solution, best_solution, objective_value, current_UB, integral_objective_value):                            
    """
    Some logic taken out of the main loop for cleanliness. Helps determine if the optimum to the TSP has been visited as a tree node.
    Inputs:
        x_solution (list):  TSP edge set vector, a solution to the subproblem of the current tree node.
        best_solution (Solution): object representing the known lowest weight for a vertex tour.
        objective_value (double): optimal value for the subproblem of the current tree node.
        current_UB (double):  upper bound for the optimal tour weight (i.e., the lowest weight found for a vertex tour, so far).
        integral_objective_value (double): optimal value for an auxiliary problem.
    Outputs:
        needs_update (boolean):  True iff a change to the visited_best_x boolean has been determined necessary by this function.
        visited_best_x_standin (boolean): True if this function has determined that the best known vertex tour has been visited as a tree node.
    """
    
    needs_update = False
    visited_best_x_standin = False
    matches_best_x = False
    fresh_obj_val_beats_UB = False
    obj_val_from_2opt_beats_fresh = False
    obj_val_from_2opt_beats_UB = False
    if [a==b for a,b in zip(x_solution, best_solution.get_x_vector())]:
        matches_best_x=True
    if objective_value < current_UB:
        fresh_obj_val_beats_UB = True
    if integral_objective_value < objective_value:
        obj_val_from_2opt_beats_fresh = True
    if integral_objective_value < current_UB:
        obj_val_from_2opt_beats_UB = True
    if matches_best_x or fresh_obj_val_beats_UB:
        needs_update = True
        if obj_val_from_2opt_beats_fresh:
            visited_best_x_standin = False
        else:
            visited_best_x_standin = True
    else:
        if obj_val_from_2opt_beats_fresh and obj_val_from_2opt_beats_UB:
            visited_best_x_standin = False
            needs_update=True
    return needs_update,visited_best_x_standin

def update_solution_corral_and_find_best(solution_corral):
    """
    We keep known vertex tours in the solution corral. In this method, we weed out unpromising tours, and try to improve existing tours with our 2-opt heuristic.
    Input: 
        solution_corral (list):  Contains Solution objects representing known vertex tours.
    Output: 
        best_sol (Solution): A reference to the best Solution object in the corral.
    """
    print("Revisiting existing solutions... ", end='')
    mark_for_del = []
    for i_s, sol in enumerate(solution_corral):
        if sol.is_at_minimum():
            mark_for_del.append(i_s)
    for i_s in mark_for_del:
        del solution_corral[i_s]
    sol_sample = random.sample(solution_corral, min(3, len(solution_corral)))
    for sol in sol_sample:
        sol.try_improvement(limit=5000)
        if sol.get_name():
            print(" "+sol.get_name()+" ", end='')
        else:
            print(" * ", end='')
    print('...done')
    sol_obj_results = [item.get_objective_value() for item in solution_corral]
    best_sol_idx = sol_obj_results.index(min(sol_obj_results))
    best_sol = solution_corral[best_sol_idx]
    if len(solution_corral) > 6:
        worst_sol_idx = sol_obj_results.index(max(sol_obj_results[:-3]))
        del solution_corral[worst_sol_idx]
    return best_sol

def final_message(x_best, current_UB, reason_to_quit, visited_best_x, g):
    """
    At termination, this function provides output describing what happened.
    Inputs:
        x_best (list):  TSP edge set vector associated with the current upper bound.
        current_UB (double):  upper bound for the optimal tour weight (i.e., the lowest weight found for a vertex tour, so far).
        reason_to_quit (string): explains why the process is terminating.
        visited_best_x (boolean): indicates whether our x_best vector, which might have been found by a heuristic method, has been found as the solution to a subproblem at a node of our tree.
        g (TSPGraph): Our graph (see problem_def.py).
    Output:
        tour (list):  list of vertex indices, sorted in order they are visited in our best vertex tour.
    """
    print(reason_to_quit)
    tour = None
    if visited_best_x:
        print("*** Optimal vertex tour has been found. ***")
        tour = makeSortedTourVerticesFromSolutionVector(g.graph_inv, x_best)
        print(tour)
        print("total weight:   "+str(current_UB))
    else:
        print("~~~  FAILURE  ~~~")
        print("Process has ended by tree exhaustion. However, the upper bound, ")
        print("discovered without the subproblem tree, was never found to be ")
        print("represented on the tree itself.")
    return tour

class vPrint():
    """Wrapper deciding whether a printed statement should appear in output."""
    def __init__(self, verbose):
        self.verbose = verbose
    def __call__(self, some_string, end='\n'):
        if self.verbose:
            print(some_string, end=end)

       





