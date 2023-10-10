"""Functions handling subgraphs. These are mostly involved in finding cutting planes that limit the feasible region."""
import random
import numpy as np
import pandas as pd

#How to get that neighbor info list:   
def get_nbr_info_list(graph_dict, main_vertex_idx_str, x=None):
    sub_dict = graph_dict[main_vertex_idx_str]
    #neighboring vertex index, solution vector value for connecting edge, edge index
    if x is not None:
        nbr_info_list = [[key, x[value[1]], value[1]] for key, value in sub_dict.items()] 
    else:
        nbr_info_list = [[key, value[1]] for key, value in sub_dict.items()] 
    return nbr_info_list

def check_if_subgraph_is_tour(tsp_graph, x_solution, is_integral):
    dummy_dict, is_tour = subgraph_a_to_b_walk(tsp_graph, x_solution, np.inf, np.inf, tour_check=True)
    if is_tour:
        if not is_integral:
            return True
        else:
            if sum(x_solution) == tsp_graph.n:
                return True
    return False

def find_random_subgraph_cut(tsp_graph, x_solution, stop_crit_new_planes = 6, stop_crit_iters = 50):
    cutting_plane_dict, dummy_bool = subgraph_a_to_b_walk(tsp_graph, x_solution, stop_crit_new_planes, \
                                                     stop_crit_iters, False)
    return cutting_plane_dict

def subgraph_a_random(tsp_graph, x_solution, stop_crit_new_planes = 6, stop_crit_iters = 50):
    stopping_criterion_number_of_new_cutting_planes = stop_crit_new_planes
    stopping_criterion_number_of_iterations = stop_crit_iters
    number_of_iterations = 0
    number_of_new_cutting_planes = 0
    nz_index_list = []
    found_edge_sets = []
    nonzero_counts_as_we_go = [0]
    while number_of_new_cutting_planes < stopping_criterion_number_of_new_cutting_planes \
        and number_of_iterations < stopping_criterion_number_of_iterations:
        cut_value = 0
        edge_indices = []
        v_ratio = .95**(number_of_iterations+12)
        my_combo_r = max(2, int(tsp_graph.n*v_ratio))
        a_vertices = np.random.choice(tsp_graph.n, my_combo_r, replace=False)
        for key, edge in tsp_graph.graph_inv.items():
            if (edge[1][0] in a_vertices) + (edge[1][1] in a_vertices) == 1:
                cut_value += x_solution[key]
                edge_indices.append(key)
        if cut_value < 2-(1e-13) and edge_indices not in found_edge_sets:
            number_of_new_cutting_planes += 1
            nz_index_list += edge_indices
            nonzero_counts_as_we_go.append(nonzero_counts_as_we_go[-1]+len(edge_indices))
            found_edge_sets.append(edge_indices)
        number_of_iterations += 1
    #make the return dictionary
    lower = np.array([2]*number_of_new_cutting_planes, dtype=np.double)  
    upper = np.array([np.inf]*number_of_new_cutting_planes, dtype=np.double) 
    num_nz = len(nz_index_list)  # for each constraint, we should have counted the number of edges 
    start = np.array(nonzero_counts_as_we_go[:-1])
    index = np.array(nz_index_list)
    value = np.array([1]*num_nz, dtype=np.double)
    return {'lower':lower, 'upper':upper, 'num_nz':num_nz, 'start':start, \
            'index':index, 'value':value}

def global_min_cut(tsp_graph, x_solution):
    graph_dict = tsp_graph.graph
    #make that matrix
    my_matrix = []
    for i_sub_dict in range(len(graph_dict)):
        sub_dict = graph_dict[i_sub_dict]
        this_row = [0 for i in range(tsp_graph.n)]
        for key, value in sub_dict.items():
            this_row[key] = x_solution[value[1]]
        my_matrix.append(this_row)
    merge_record = [[i] for i in range(tsp_graph.n)]
    while len(merge_record)>2:
        #find highest valued edge
        max_edge = -np.inf
        max_row = -1
        max_col = -1
        for i_row in range(len(my_matrix)):
            expurg_list = my_matrix[i_row][:i_row] + [-1e20] + my_matrix[i_row][i_row+1:]
            max_idx_this_row = expurg_list.index(max(expurg_list))
            if my_matrix[i_row][max_idx_this_row] > max_edge:
                max_edge = my_matrix[i_row][max_idx_this_row]
                max_row = i_row
                max_col = max_idx_this_row
        merge_record[max_row] += merge_record[max_col]
        merge_record.pop(max_col)
        next_max_iter = len(my_matrix)
        for i_row in range(next_max_iter):
            my_matrix[i_row][max_row] += my_matrix[i_row][max_col]
            my_matrix[i_row].pop(max_col)
        my_matrix[max_row] = [a+b for a,b in zip(my_matrix[max_row],my_matrix[max_col])]
        my_matrix.pop(max_col)
    coeff_row = [0 for i in range(tsp_graph.m)]
    pos_edge_indices = []
    merge_record_lengths = [len(merge_record[0]), len(merge_record[1])]
    v_list_to_use = merge_record[merge_record_lengths.index(min(merge_record_lengths))]
    for v_idx in v_list_to_use:
        pos_edge_indices += [value[1] for key, value in graph_dict[v_idx].items()]
    edge_idx_list_to_use = pd.Series(pos_edge_indices).drop_duplicates(keep=False).tolist()
    vec_total_actual = 0
    for e_idx in edge_idx_list_to_use:
        coeff_row[e_idx]=1
        vec_total_actual += x_solution[e_idx]
    return coeff_row, vec_total_actual
    
def subgraph_a_to_b_walk(tsp_graph, x_solution, \
                             stopping_criterion_number_of_new_cutting_planes = 6, \
                             stopping_criterion_number_of_iterations = 50, \
                             tour_check = False):
    graph_dict = tsp_graph.graph
    is_tour = False
    closed_loop = False
    break_flag = False
    number_of_iterations = 0            #counter
    number_of_new_cutting_planes = 0    #counter
    nonzero_counts_as_we_go = [0]       #for the "start" LP model creation argument
    nz_index_list = [] #in the LP constraint coefficient matrix, these are the variable indices of
    #                   non-zero terms, row by row as a single vector
    
    while True:
        v_0 = random.randint(0, tsp_graph.n-1)     #choose random start
        #print("new random start: "+str(v_0))
        my_a_dictionary = {v_0:get_nbr_info_list(graph_dict, v_0, x_solution)} #start with just one key-value pair
        hello_flag = 0
        while True:
            #print("HELLO   "+str(hello_flag))
            #b_neighbors:  nbr info for all neighbors in B
            #nz_b_neighbors: the above, but with positive solution vector values
            #nz_b_neighbors_per_member_of_a: the number of new B neighbors gained by each A member
            b_neighbors = [elem for key, value in my_a_dictionary.items() for elem in value \
                        if elem[0] not in my_a_dictionary.keys()]
            nz_b_neighbors = [neighbor for neighbor in b_neighbors if neighbor[1] > 1e-13]
            nz_b_neighbors_per_member_of_a = [len([elem for elem in value \
                                                if elem[0] not in my_a_dictionary.keys() \
                                                    and elem[1] > 1e-13]) \
                                                for key, value in my_a_dictionary.items()]
            #print("b_neighbors:")
            #print(b_neighbors)
            #print("nz_b_neighbors:")
            #print(nz_b_neighbors)
            #print("nz_b_neighbors_per_member_of_a:")
            #print(nz_b_neighbors_per_member_of_a)
            if tour_check:
                #How many neighbors are there per member of A? The correct value is 1,
                # unless this is the first iteration, in which case 2.
                if len(nz_b_neighbors) != 2 or \
                    not all([item <= (max(1, 3-len(my_a_dictionary))) \
                             for item in nz_b_neighbors_per_member_of_a]):
                    break_flag = True
                    break
                nz_b_nbr_zero_dict = graph_dict[nz_b_neighbors[0][0]]
                if nz_b_neighbors[0][0] == nz_b_neighbors[1][0]:
                    closed_loop = True
                elif nz_b_neighbors[1][0] in nz_b_nbr_zero_dict.keys():
                    loop_closing_value = nz_b_nbr_zero_dict[nz_b_neighbors[1][0]][1]
                    if x_solution[loop_closing_value] > 1e-13:
                        closed_loop = True
            #At this point, "my_A_dictionary" has, as its keys, the indices of all
            # starting members of A. Its values are lists of neighbor information.
            #Get the solution vector total for the A-B cut. 
            cut_value = sum([item[1] for item in b_neighbors])
            #print("A:  ", end='')
            #print(my_a_dictionary.keys())
            #print("B:  ", end='')
            #print([elem[0] for elem in b_neighbors])
            #print("cut_value:   "+str(cut_value), end='')
            edge_indices_for_this_cut = [item[2] for item in b_neighbors]
            edge_indices_for_this_cut.sort()
            if cut_value < 2-(1e-13):
                number_of_new_cutting_planes += 1
                nz_index_list += edge_indices_for_this_cut
                cumulative_edge_index_count = \
                    len(edge_indices_for_this_cut) + nonzero_counts_as_we_go[-1]
                nonzero_counts_as_we_go.append(cumulative_edge_index_count)
                #print(number_of_new_cutting_planes)
                #print(nz_index_list)
                #print(nonzero_counts_as_we_go)
            #Bring the B neighbors into set A.
            # If we're not doing a tour check, we also bring in
            #  B neighbors with zero x solution values.
            #  We're actually going to try just bringing in the maximum,
            # whatever it may be, as in the min cut algorithm.
            if not tour_check:
                connecting_values = [item[1] for item in b_neighbors]
                list_idx_of_max = connecting_values.index(max(connecting_values))
                vertex_idx_of_max = b_neighbors[list_idx_of_max][0]
                my_a_dictionary.update({vertex_idx_of_max:get_nbr_info_list(graph_dict, \
                                    vertex_idx_of_max, x_solution)})
            else:
                my_a_dictionary.update({vertex_idx_str:get_nbr_info_list(graph_dict, \
                                    vertex_idx_str, x_solution) for vertex_idx_str \
                                        in [item[0] for item in b_neighbors \
                                            if item[1] > 1e-13]})

            #print("A nodes:")
            #print(my_a_dictionary.keys())
            
            number_of_iterations += 1

            if len(my_a_dictionary) == tsp_graph.n: #all vertices used
                #print("A-B subgraph: all vertices used")
                if tour_check:
                    if closed_loop:
                        is_tour = True
                    break_flag = True
                break
            elif tour_check:
                if closed_loop:
                    break_flag = True
                    break
            if number_of_iterations >= stopping_criterion_number_of_iterations:
                #print("A-B subgraph: no more iterations")
                break_flag = True
                break
            if number_of_new_cutting_planes >= stopping_criterion_number_of_new_cutting_planes:
                #print("A-B subgraph: found enough cutting planes")
                break_flag = True
                break
            hello_flag += 1
        if break_flag:
            break
        if number_of_iterations >= stopping_criterion_number_of_iterations:
            break
        if number_of_new_cutting_planes >= stopping_criterion_number_of_new_cutting_planes:
            break    
    if tour_check:
        return {}, is_tour
    
    #if not a tour check, make the return dictionary
    lower = np.array([2]*number_of_new_cutting_planes, dtype=np.double)  
    upper = np.array([np.inf]*number_of_new_cutting_planes, dtype=np.double) 
    num_nz = len(nz_index_list)  # for each constraint, we should have counted the number of edges 
    start = np.array(nonzero_counts_as_we_go[:-1])
    index = np.array(nz_index_list)
    value = np.array([1]*num_nz, dtype=np.double)
    return {'lower':lower, 'upper':upper, 'num_nz':num_nz, 'start':start, \
            'index':index, 'value':value}, False

