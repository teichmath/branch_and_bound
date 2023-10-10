import random
import numpy as np

class choiceStrategy:
    """
    Strategies for selecting the next subproblem in the queue.
    """
    def __init__(self, method_str):
        self.method_str = method_str
    def __call__(self, LP_deque, live_tree_node_list):    
        """
        Inputs:
            LP_deque (list): the subproblem queue.
            live_tree_node_list (list):  each element represents a node of the subproblem tree that has children with unknown sub-LP solutions.
        Output:
            next_LP (list): the next subproblem to use. item at index i is either 0 or 1, giving the rhs of an equality constraint for edge x_i.
        """
        next_LP = []
        if self.method_str=="straight":
            """choose the next LP in the queue."""
            next_LP = LP_deque.popleft()
        if self.method_str=="random":
            """choose a random LP in the queue."""
            next_LP = random.sample(LP_deque, 1)[0]
            LP_deque.remove(next_LP)
        if self.method_str=="favor LB":
            """find the node representing the known lower bound of the optimal objective value; choose an LP at a child node."""
            if len(live_tree_node_list) > 0:
                live_obj_value_list = [item['obj_value'] for item in live_tree_node_list]
                live_node_idx_for_LB = live_obj_value_list.index(min(live_obj_value_list))
                LB_rep = live_tree_node_list[live_node_idx_for_LB]['x_lead']
                all_LB_rep_LPs = [item for item in LP_deque if len(item)>0 and all([a==b for a,b in zip(LB_rep, item[:len(LB_rep)])])]
                #next_LP = random.choice(all_LB_rep_LPs)
                next_LP = all_LB_rep_LPs[0]
                LP_deque.remove(next_LP)
            else:
                next_LP = LP_deque.popleft()
        if self.method_str=="favor max":
            """find the live node with the greatest objective value among such nodes; choose an LP at a child node."""
            if len(live_tree_node_list) > 0:
                highest_repping_LP_idx = 0
                repped_value = 0
                live_obj_value_list = [item['obj_value'] for item in live_tree_node_list]
                LP_idx_list = np.random.choice(len(LP_deque), min(15, len(LP_deque)), replace=False)
                for i_lp in LP_idx_list:
                    LP_item = LP_deque[i_lp]
                    for live_item in live_tree_node_list:
                        if len(live_item) <= len(LP_item):
                            if all([a==b for a,b in zip(live_item, LP_item[:len(live_item)])]):
                                if live_item['obj_value'] > repped_value:
                                    repped_value = live_item['obj_value']
                                    highest_repping_LP_idx = i_lp
                next_LP = LP_deque[highest_repping_LP_idx]
                LP_deque.remove(next_LP)
            else:
                next_LP = LP_deque.popleft()
        return next_LP

class LPChoiceManager():
    """
    Manages the process of rotating subproblem choice strategies.
    """
    def __init__(self, strategy_list, share_list, cycle_limit):
        self.LP_choice_strategies = [{'strategy':choiceStrategy(a), 'share':b} for a,b in zip(strategy_list, share_list)]
        self.LP_current_strategy_idx = 0
        self.LP_current_strategy = self.LP_choice_strategies[self.LP_current_strategy_idx]['strategy']
        self.cycle_limit = cycle_limit
        self.first_cycle = True
        self.improvement_clock = {'counter':0, 'limit':cycle_limit*self.LP_choice_strategies[self.LP_current_strategy_idx]['share'], 'max limit':self.cycle_limit}
    def get_current_strategy(self):
        return self.LP_current_strategy
    def time_is_up(self):
        if self.improvement_clock['counter'] == self.improvement_clock['limit']:
            return True
        else:
            return False
    def rotate_strategy(self):
        if self.time_is_up():
            self.LP_current_strategy_idx += 1
            if self.LP_current_strategy_idx >= len(self.LP_choice_strategies):
                self.LP_current_strategy_idx = 0
            self.LP_current_strategy = self.LP_choice_strategies[self.LP_current_strategy_idx]['strategy']
            self.improvement_clock['counter'] = 0
            self.improvement_clock['limit'] = int(self.improvement_clock['max limit']*self.LP_choice_strategies[self.LP_current_strategy_idx]['share'])
            self.first_cycle = False
    def new_cycle(self):
        if self.improvement_clock['counter'] == 0 and self.LP_current_strategy_idx == 0 and not self.first_cycle:
            return True
        else:
            return False
    def advance_clock(self, delta):
        self.improvement_clock['counter'] += delta
    def turn_back_clock(self, delta):
        self.improvement_clock['counter'] = max(0, self.improvement_clock['counter']-delta)