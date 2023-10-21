"""branch and cut logic, plus utility classes."""
import time
from lp_choice import LPChoiceManager
from problem_def import get_single_data, BnBGraph, LPHandler
from search_account import BnBSearchAccount

def bnb_driver(graph_name, search_config, verbosity=1):
    """
    main function of the branch-and-bound solver. handles initialization and flow control.
    Inputs:
        graph_name (string): path to text file containing graph description. See README.md.
        search_config (dict): here we pass multiple functions to the driver to specialize it 
                              for solving a particular graph problem (e.g. the TSP, or the Max Cut).
        verbosity (int): 0: no messages other than errors and result;  
                         1: lower and upper bounds are reported as discovered;
                         2: annotations describing most logical decisions.
    Outputs:
        a solution (Solution): object representing the optimal problem solution.    
    """
    
    print("Solving:   "+search_config['problem name']+"    for graph:  "+search_config['graph name'])
    
    start = time.time()
    minimize = True
    if search_config['optimization type'] == 'maximize':
        minimize = False

    #get input
    input_data = get_single_data(graph_name)

    #initialize search, get initial upper bound
    vprint = VPrint(verbosity)
    bnb_graph = BnBGraph(input_data)
    bnb_lp = LPHandler(bnb_graph, search_config['initial boundary finder'], vprint, minimize)
    search_accountant = BnBSearchAccount(bnb_graph, minimize, vprint)
    search_accountant.initialize_upper_bound(search_config['greedy heuristic'])
    choice_manager = LPChoiceManager(["favor LB", "random"], [0.5, 0.02], 100)

    while True:
        #strategy management: possibly change LP choice strategy
        choice_manager.rotate_strategy()
        #possibly revisit valid solutions to attempt improvement
        if choice_manager.new_cycle():
            search_accountant.revisit_corral_to_improve_and_cut(search_config['greedy heuristic'], search_config['improvement heuristic'])
        #Choose an LP
        sub_lp = search_accountant.get_next_LP(choice_manager)
        vprint("Now handling sub LP:", 2)
        vprint(sub_lp, 2)
        #Check if this tree node is pruned, then proceed
        if search_accountant.tree_node_is_pruned(sub_lp):
            vprint("This tree node is pruned, continuing.", 2)
            break
        else:
            #initialize LP relaxation
            vprint("initializing lp.", 2, end= '  ')
            bnb_lp.initialize_lp()
           
            #add constraints to the LP
            bnb_lp.define_domain(sub_lp)

            #solve LP
            vprint("Solving LP...   ", 2, end='  ')
            objective_value, is_feasible = bnb_lp.solve_lp()
            search_accountant.process_new_tree_node(objective_value, is_feasible, sub_lp)

            #if branch continues here
            if search_accountant.tree_node_is_live(sub_lp):
                vprint("This tree node is live.", 2)
                x_solution = bnb_lp.get_existing_solution()
                #is the solution integral, valid
                is_integral, is_valid = search_accountant.check_validity(x_solution, search_config['edge set validator'])
                if is_valid:
                    if is_integral:
                        vprint("Found a "+search_config['valid edge set type']+"!", 1)
                        search_accountant.add_tree_node_solution_to_solution_corral(objective_value, x_solution)
                    else:
                        vprint("The solution is non-integral, but it provides a "+search_config['valid edge set type']+"!", 1)
                        integral_x, integral_obj_val = bnb_lp.get_proximal_integral_result_from_nonintegral(x_solution)
                        search_accountant.add_off_tree_solution_to_solution_corral(integral_obj_val, integral_x)
                if not (is_integral and is_valid):
                    vprint("Solution is not itself a "+search_config['valid edge set type']+", so: ", 2)
                    #Get cutting plane(s)
                    vprint("1.   Looking for cutting planes...", 2, end='')
                    bnb_lp.find_new_cutting_planes(search_config['cutting plane finder'], x_solution)
                    #branch
                    if len(sub_lp) == bnb_graph.m:
                        search_accountant.remove_live_nodes_by_sub_lp(sub_lp)
                    else:
                        vprint("3.   Branching.", 2)
                        search_accountant.branch(sub_lp)
        search_accountant.check_for_termination()
        if not search_accountant.is_search_running():
            print(search_accountant.get_reason_to_quit())
            if search_accountant.best_x_was_visited():
                end = time.time()
                print("*** Optimal "+search_config['valid edge set type']+" has been found. ***")
                vertex_set = search_config['vertex set getter'](bnb_graph, search_accountant.get_best_solution())
                search_accountant.final_display(end-start, search_config['vertex set description'], vertex_set)
            else:
                print("~~~  FAILURE  ~~~")
                print("Process has ended by tree exhaustion. However, the upper bound, ")
                print("discovered without the subproblem tree, was never found to be ")
                print("represented on the tree itself.")
            return search_accountant.get_best_solution()
        choice_manager.advance_clock(1)
    return

class VPrint():
    """Wrapper deciding whether a printed statement should appear in output."""
    def __init__(self, verbosity):
        self.verbosity = verbosity
    def __call__(self, some_string, code, end='\n'):
        if code <= self.verbosity:
            print(some_string, end=end)



"""
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
"""