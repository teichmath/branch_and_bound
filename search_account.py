"""
Here we do the bookkeeping. All repeated actions in the branch-and-bound algorithm, especially those 
that update our records, are taken out of the driver's flow and placed here.
"""

from collections import deque
import random
import numpy as np
from problem_def import Solution

class BnBSearchAccount():
    """ stores and writes to subproblem tree information: what nodes we're still remembering, 
        what nodes we want to visit next, what branches have been pruned, and more. """
    def __init__(self, bnb_graph, minimize, vprint):
        self.graph = bnb_graph
        self.vprint = vprint
        self.minimize = minimize
        self.display_objective_sign = 1
        self.UB_name = 'UB'
        self.LB_name = 'LB'
        if not self.minimize:
            self.display_objective_sign = -1
            self.UB_name = 'LB'
            self.LB_name = 'UB'
        self.subprob_queue = deque([[]])  #[0],[1],[0,0],[0,1], etc.
        self.solution_corral = []
        self.visited_tree_node_list = []
        # a "live" tree node is a node that is not a terminus, and
        # that gives a greatest known lower bound (or least known
        # upper bound, in maximization) for all its descendant subproblems.
        self.live_tree_node_list = []
        self.pruned_tree_node_list = []
        self.visited_best_x = False
        self.search_is_running = True
        self.reason_to_quit = None
        self.current_UB = np.inf
        self.current_LB = -np.inf
    def initialize_upper_bound(self, greedy_heuristic):
        self.vprint("initializing upper bound. (24)", 2)
        self.add_fresh_solution_to_solution_corral(greedy_heuristic, 'initial')
    def get_next_LP(self, choice_manager):
        self.vprint("getting the next LP.  (27)", 2)
        next_LP = choice_manager.get_current_strategy()(self.subprob_queue, self.live_tree_node_list)
        self.vprint("removing LP from subprob queue:   "+str(next_LP), 2)
        self.subprob_queue.remove(next_LP)
        self.visited_tree_node_list.append(next_LP)
        return next_LP
    def update_UB(self, some_solution, from_tree):
        obj_value = some_solution.get_objective_value()
        self.vprint("old UB is "+str(self.display_objective_sign * self.current_UB)+\
                    ", new proposed UB is "+str(self.display_objective_sign*obj_value)+"   (34)", 2)
        if obj_value < self.current_UB:
            self.vprint("We have a new UB. Updating x_best and current_UB.  (35)", 2)
            self.current_UB = obj_value
            self.vprint(self.LB_name+": "+str(self.display_objective_sign*self.current_LB)+\
                        "       "+self.UB_name+": "+str(self.display_objective_sign*self.current_UB), 1)
            pop_these_live_nodes = []
            for i_tnd, tree_node_dict in enumerate(self.live_tree_node_list):
                if tree_node_dict['obj_value'] >= self.current_UB:
                    ineq_sign = '>='
                    if not self.minimize:
                        ineq_sign = '<='
                    self.vprint("Pruning tree node for being "+ineq_sign+" current "+self.UB_name+":", 2, end='')
                    self.vprint(tree_node_dict, 2)
                    pop_these_live_nodes.append(i_tnd)
                    self.pruned_tree_node_list.append(tree_node_dict)
            if len(pop_these_live_nodes)>0:
                self.remove_live_nodes_by_idx(pop_these_live_nodes)
        if obj_value <= self.current_UB:
            self.vprint("setting visited_best_x_ to :  ", 2, end=' ')
            self.vprint(from_tree, 2, end=' ')
            self.vprint("(48)", 2)
            self.visited_best_x = from_tree
    def update_LB(self):
        self.vprint("possibly updating "+self.LB_name+".  (51)", 2)
        if len(self.live_tree_node_list)>0:
            live_node_min = min([item['obj_value'] for item in self.live_tree_node_list])
            if live_node_min > self.current_LB:
                self.current_LB = live_node_min
                self.vprint(self.LB_name+": "+str(self.display_objective_sign*self.current_LB)+\
                        "       "+self.UB_name+": "+str(self.display_objective_sign*self.current_UB), 1)
    def branch(self, sub_lp):
        self.vprint("adding to subprob queue:  "+str(sub_lp)+" and 0, 1.", 2)
        self.subprob_queue.append(sub_lp + [0])
        self.subprob_queue.append(sub_lp + [1])
    def check_validity(self, x_solution, edge_set_validator):
        is_integral = False
        is_valid = False
        self.vprint("Is the solution integral ?    ", 2, end='')
        if sum([abs(item - round(item)) for item in x_solution]) < 1e-13:
            is_integral = True
            self.vprint("Yes", 2)
        else:
            self.vprint("No", 2)
        #check if it's a tour
        self.vprint("Is the solution valid?    ", 2, end='')
        is_valid = edge_set_validator(self.graph, x_solution, is_integral)
        if is_valid:
            self.vprint("Yes", 2)
        else:
            self.vprint("No", 2)
        return is_integral, is_valid
    def check_for_termination(self):
        self.vprint("checking for termination.  (79)", 2)
        if len(self.subprob_queue)==0:
            self.reason_to_quit = "subprob queue is empty."
            self.search_is_running = False
        if abs(self.current_LB-self.current_UB) < 1-5e-2:
            self.reason_to_quit = "optimality gap is zero."
            self.search_is_running = False
    def is_search_running(self):
        return self.search_is_running
    def get_reason_to_quit(self):
        return self.reason_to_quit
    def revisit_corral_to_improve_and_cut(self, greedy_heuristic, improvement_heuristic):
        self.vprint("Revisiting existing solutions... ", 2, end='')
        sol_sample = random.sample(self.solution_corral, min(3, len(self.solution_corral)))
        for sol in sol_sample:
            sol.try_improvement(improvement_heuristic, limit=5000)
            self.update_UB(sol, from_tree=False)
            if sol.get_name():
                self.vprint(" "+sol.get_name()+" ", 2, end='')
            else:
                self.vprint(" * ", 2, end='')
        self.vprint('...done', 2)
        sol_obj_results = [item.get_objective_value() for item in self.solution_corral]
        if len(self.solution_corral) > 6:
            worst_sol_idx = sol_obj_results.index(max(sol_obj_results[:-3]))
            del self.solution_corral[worst_sol_idx]
        self.add_fresh_solution_to_solution_corral(greedy_heuristic)
    def tree_node_is_live(self, sub_lp):
        is_live = False
        self.vprint("is tree node live?   (110)", 2)
        for tree_node in self.live_tree_node_list:
            if len(sub_lp)==len(tree_node['x_lead']):
                if len(sub_lp)==0:
                    is_live = True
                    break
                if all([a==b for a,b in zip(tree_node['x_lead'], sub_lp)]):
                    is_live = True
                    break
        return is_live
    def tree_node_is_pruned(self, sub_lp):
        is_pruned = False
        for pruned_node in self.pruned_tree_node_list:
            if len(pruned_node) <= len(sub_lp):
                if all(x1 == x2 for x1, x2 in zip(pruned_node, sub_lp[:len(pruned_node)])):
                    is_pruned = True
                    break
        return is_pruned
    def process_new_tree_node(self, objective_value, is_feasible, sub_lp):
        self.vprint("processing a new tree node.  (126)", 2)
        self.vprint("current "+self.UB_name+" is:  "+str(self.display_objective_sign*self.current_UB), 2)
        this_tree_node = {'obj_value':objective_value, 'x_lead':sub_lp}
        if not is_feasible:
            self.vprint("Not feasible. Pruning and continuing.", 2)
            self.pruned_tree_node_list.append(this_tree_node)
        else:
            self.vprint("Feasible with objective value:  "+\
                        str(self.display_objective_sign*objective_value), 2)
            if objective_value > self.current_UB:
                self.vprint("This is greater than the current "+self.UB_name+". Pruning and continuing.", 2)
                self.pruned_tree_node_list.append(this_tree_node)
            else:
                self.vprint("This is less than or equal to the current "+self.UB_name+".", 2)
                self.vprint("adding to live tree node list:  "+str(this_tree_node), 2)
                self.live_tree_node_list.append(this_tree_node)
        if len(sub_lp) > 0:
            this_tree_node_twin = sub_lp[:-1] + [1 - sub_lp[-1]] 
            if this_tree_node_twin in self.visited_tree_node_list:
                current_x_leads = [item['x_lead'] for item in self.live_tree_node_list]
                if sub_lp[:-1] in current_x_leads:
                    current_x_leads.remove(sub_lp[:-1])
                    removal_list = [idx for idx, item in enumerate(self.live_tree_node_list) if item['x_lead'] not in current_x_leads]
                    self.remove_live_nodes_by_idx(removal_list)
    def remove_live_nodes_by_sub_lp(self, sub_lp):
        x_lead_list = [item['x_lead'] for item in self.live_tree_node_list]
        self.remove_live_nodes_by_idx([x_lead_list.index(sub_lp)])
    def remove_live_nodes_by_idx(self, nodes):
        for i_tnd in sorted(nodes, reverse=True):
            self.vprint("removing live node at idx "+str(i_tnd)+", which is "+str(self.live_tree_node_list[i_tnd]), 2)
            del self.live_tree_node_list[i_tnd]
        self.update_LB()
    def add_fresh_solution_to_solution_corral(self, greedy_heuristic, name=None):
        if name is None:
            name = str(len(self.visited_tree_node_list))
        some_solution = self.get_fresh_solution(greedy_heuristic, name)
        self.add_to_solution_corral(some_solution, from_tree=False)
    def add_tree_node_solution_to_solution_corral(self, obj_value, x_solution):
        this_sol = Solution(self.graph, obj_value, x_solution, name=str(len(self.visited_tree_node_list))+"B", minimize=self.minimize)
        if obj_value <= self.current_UB:
            self.add_to_solution_corral(this_sol, from_tree=True)
    def add_off_tree_solution_to_solution_corral(self, obj_value, x_solution):
        this_sol = Solution(self.graph, obj_value, x_solution, name=str(len(self.visited_tree_node_list))+"B", minimize=self.minimize)
        if obj_value < self.current_UB:
            self.add_to_solution_corral(this_sol, from_tree=False)
    def add_to_solution_corral(self, some_solution, from_tree=False):
        self.solution_corral.append(some_solution)
        self.update_UB(some_solution, from_tree)
    def get_fresh_solution(self, greedy_heuristic, name):
        some_solution = Solution(self.graph, name=name, heuristic=greedy_heuristic, minimize=self.minimize)
        return some_solution
    def get_best_solution(self):
        best_obj = np.inf
        best_sol = None
        for sol in self.solution_corral:
            if sol.get_objective_value() < best_obj:
                best_obj = sol.get_objective_value()
                best_sol = sol
        return best_sol
    def best_x_was_visited(self):
        return self.visited_best_x
    def final_display(self, total_time, vertex_set_description, vertex_set):
        print("Optimal weight:       "+str(self.display_objective_sign*self.current_UB))
        print("Time elapsed:         "+str(total_time))
        print("Tree nodes visited:   "+str(len(self.visited_tree_node_list)))
        print(vertex_set_description+"    "+str(vertex_set))
        print("Optimum:   ")
        print(self.get_best_solution().get_x_vector())

