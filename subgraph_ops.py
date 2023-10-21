"""Functions handling subgraphs. These are mostly involved in finding cutting planes that limit the feasible region."""
import random
import copy
import math
import numpy as np
import pandas as pd

def get_nbr_info_list(graph_dict, main_vertex_idx_str, x=None):
    sub_dict = graph_dict[main_vertex_idx_str]
    #neighboring vertex index, solution vector value for connecting edge, edge index
    if x is not None:
        nbr_info_list = [[key, x[value[1]], value[1]] for key, value in sub_dict.items()] 
    else:
        nbr_info_list = [[key, value[1]] for key, value in sub_dict.items()] 
    return nbr_info_list

def check_if_subgraph_is_tour(graph, x_solution, is_integral):
    _, is_tour = subgraph_a_to_b_walk(graph, x_solution, np.inf, np.inf, tour_check=True)
    if is_tour:
        if not is_integral:
            return True
        else:
            if sum(x_solution) == graph.n:
                return True
    return False

def get_final_tour_vertices(graph, solution):
    tour = makeSortedTourVerticesFromSolutionVector(graph.graph_inv, solution.get_x_vector())
    return tour

def get_final_vertices_with_optimal_cut_for_delta(graph, solution):
    x_solution = solution.get_x_vector()
    is_integral = False
    if sum([abs(item - round(item)) for item in x_solution]) < 1e-13:
        is_integral = True
    _, vertices = get_subgraph_delta_inverse(graph, x_solution, is_integral)
    return vertices

def check_if_edge_set_is_cut(graph, x_solution, is_integral):
    reply, _ = get_subgraph_delta_inverse(graph, x_solution, is_integral)
    return reply

def get_subgraph_delta_inverse(graph, x_solution, is_integral):
    if not is_integral:
        return False, []
    edge_set = [i for i in range(graph.m) if x_solution[i]==1]
    edge_complement = [i for i in range(graph.m) if x_solution[i]==0]
    vertex_set_A = get_vertices_of_connected_piece(graph, edge_complement, 0)
    delta_A = get_subgraph_delta(graph.graph, vertex_set_A)
    if len(edge_set)==len(delta_A):
        if all([a==b for a,b in zip(sorted(edge_set), sorted(delta_A))]):
            return True, vertex_set_A
    neighborhood = list(set([key for v in vertex_set_A for key in graph.graph[v].keys()]))
    neighborhood = [v for v in neighborhood if v not in vertex_set_A]
    neighborhood_extension_by_edge_complement = []
    for vertex in neighborhood:
        neighborhood_extension_by_edge_complement += get_vertices_of_connected_piece(graph, edge_complement, vertex)
    neighborhood = list(set(neighborhood + neighborhood_extension_by_edge_complement))
    absorbable_vertices = [v for v in range(graph.n) if v not in neighborhood]
    if len(absorbable_vertices)==0:
        return False, []
    true_A = list(set(vertex_set_A+absorbable_vertices))
    delta_true_A = get_subgraph_delta(graph.graph, true_A)
    if len(edge_set)==len(delta_true_A):
        if all([a==b for a,b in zip(sorted(edge_set), sorted(delta_true_A))]):
            return True, true_A
    return False, []

def get_vertices_of_connected_piece(graph, edge_set, starting_vertex):
    vertex_set = [starting_vertex]
    while True:
        neighbor_set = []
        for v in vertex_set:
            for neighbor_idx, edge in graph.graph[v].items():
                if edge[1] in edge_set:
                    neighbor_set.append(neighbor_idx)
        former_size = len(vertex_set)
        vertex_set = list(set(vertex_set+neighbor_set))
        if len(vertex_set) == former_size:
            break
    return vertex_set

def find_cycles(graph, x_solution):
    cycles=[]
    edge_supply = copy.deepcopy(x_solution)
    #edge_dict = copy.deepcopy(graph.graph_inv)
    while sum(edge_supply)>0:
        some_path = makeSortedTourVerticesFromSolutionVector(graph.graph_inv, edge_supply)
        if some_path[0]==some_path[-1]:
            cycles.append(some_path)
        for key, edge in graph.graph_inv.items():
            if edge[1][0] in some_path or edge[1][1] in some_path:
                edge_supply[key]=0
    return cycles

def try_a_cut_increase(graph, v_set, start_vertex):
    vertex_set = copy.deepcopy(v_set)
    full_key_list = list(graph.graph.keys())
    key_list = full_key_list[start_vertex:] + full_key_list[:start_vertex]
    for v_key in key_list:
        v_dict = graph.graph[v_key]
        my_worth = 0
        for key, value in v_dict.items():
            if key in vertex_set:
                my_worth -= value[1]
            else:
                my_worth += value[1]
        if v_key in vertex_set and my_worth < 0:
            vertex_set.remove(v_key)
        if v_key not in vertex_set and my_worth > 0:
            vertex_set.append(v_key)
    return vertex_set

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
        if v_ratio is None:
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

def get_subgraph_delta(graph_dict, v_list):
    pos_edge_indices = []
    for v_idx in v_list:
        pos_edge_indices += [value[1] for key, value in graph_dict[v_idx].items()]
    edge_idx_list_to_use = pd.Series(pos_edge_indices).drop_duplicates(keep=False).tolist()
    return edge_idx_list_to_use

def global_min_cut(tsp_graph, x_solution):
    merge_record, _ = edge_contractions_max_to_min(tsp_graph, x_solution)
    coeff_row = [0 for i in range(tsp_graph.m)]
    merge_record_lengths = [len(merge_record[0]), len(merge_record[1])]
    v_list_to_use = merge_record[merge_record_lengths.index(min(merge_record_lengths))]
    delta_edge_set = get_subgraph_delta(tsp_graph.graph_dict, v_list_to_use)
    vec_total_actual = 0
    for e_idx in delta_edge_set:
        coeff_row[e_idx]=1
        vec_total_actual += x_solution[e_idx]
    return coeff_row, vec_total_actual

def edge_contractions_max_to_min(tsp_graph, x_solution):
    graph_dict = tsp_graph.graph
    my_matrix = []  # like an adjacency matrix; has x_solution values as its entries.
    my_edge_sets_matrix = [] # like an adjacency matrix; has edge indices as its entries.
    cliques = []
    merge_record = [[i] for i in range(tsp_graph.n)]
    #create matrices. for each vertex, ...
    for i_sub_dict in range(len(graph_dict)):
        sub_dict = graph_dict[i_sub_dict]
        this_row = [0 for i in range(tsp_graph.n)]
        this_edge_row = [[] for i in range(tsp_graph.n)]
        #...for each neighbor, record the solution vector value and index of the connecting edge.
        for key, value in sub_dict.items():
            this_row[key] = x_solution[value[1]]
            this_edge_row[key].append(value[1])
        my_matrix.append(this_row)
        my_edge_sets_matrix.append(this_edge_row)
    #contract edges iteratively
    while len(merge_record)>3:
        max_edge = -np.inf
        max_row = -1
        max_col = -1
        #get the maximally valued edge that is not on the diagonal
        for i_row in range(len(my_matrix)):
            expurg_list = my_matrix[i_row][:i_row] + [-1e20] + my_matrix[i_row][i_row+1:]
            max_idx_this_row = expurg_list.index(max(expurg_list))
            if my_matrix[i_row][max_idx_this_row] > max_edge:
                max_edge = my_matrix[i_row][max_idx_this_row]
                max_row = i_row
                max_col = max_idx_this_row
        #we will now contract the edge connecting the vertices represented by the max row and the max col.
        #record the contraction by grouping the affected vertices
        merge_record[max_row] += merge_record[max_col]
        merge_record.pop(max_col)
        next_max_iter = len(my_matrix)
        #in both matrices (carrying the solution's edge values and the edge indices, respectively), group together entries to 
        # record the contraction
        for each_matrix in [my_matrix, my_edge_sets_matrix]:
            for i_row in range(next_max_iter):
                each_matrix[i_row][max_row] += each_matrix[i_row][max_col]
                each_matrix[i_row].pop(max_col)
            each_matrix[max_row] = [a+b for a,b in zip(each_matrix[max_row],each_matrix[max_col])]
            each_matrix.pop(max_col)
        my_edges = sorted(list(set(my_edge_sets_matrix[max_row][max_row])))
        clique_size = len(merge_record[max_row])
        if len(my_edges) == clique_size*(clique_size-1)/2:
            clique_weight_limit = math.floor((clique_size**2)/4) + 1e-14
            cliques.append({'size':clique_size, 'weight':my_matrix[max_row][max_row]/2, \
                        'limit':clique_weight_limit, 'edges':my_edges})
    return merge_record, cliques


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
        #print("@@@@@@@GOING@@@@@@@@")
        v_0 = random.randint(0, tsp_graph.n-1)     #choose random start
        #print("new random start: "+str(v_0))
        my_a_dictionary = {v_0:get_nbr_info_list(graph_dict, v_0, x_solution)} #start with just one key-value pair
        hello_flag = 0
        while True:
            #print("#######GOING########")
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
            #print("222  b_neighbors:")
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
            #print("247   A:  ", end='')
            #print(my_a_dictionary.keys())
            #print("      B:  ", end='')
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
                #print("line 260")
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
            #print("281   A nodes:")
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

def makeSolutionVectorFromSortedTourEdges(graph_inv, tour_edges):
    """
    Inputs:
        graph_inv (dict): edge-wise description of the graph. See problem_def.py.
        tour_edges (list): indices of edges, listed in order visited to perform a vertex tour.
    Output:
        x_solution (list): edge set indicator vector.
        """
    x_solution = [0 for i in range(len(graph_inv))]
    for key, edge in graph_inv.items():
        if edge in tour_edges:
            x_solution[key]=1
    return x_solution

def makeSortedTourVerticesFromSolutionVector(graph_inv, x_solution):
    """
    Inputs:
        graph_inv (dict): edge-wise description of the graph. See problem_def.py.
        x_solution (list): edge set indicator vector for a vertex tour.
    Output:
        tour (list): indices of graph vertices, listed in order visited to perform a vertex tour. Final value equals first value.
        """
    # list of all edge indices in binary x solution
    edge_indices = [i for i in range(len(x_solution)) if x_solution[i]==1]
    # initialized tour has the first and second vertices belonging to the first edge
    tour = [graph_inv[edge_indices[0]][1][0], graph_inv[edge_indices[0]][1][1]]
    # performing an action edge-many times
    for i_edge_idx in range(1, len(edge_indices)):
        # action involves each edge...
        for edge_idx_prime in edge_indices:
            #if this edge has the end of the tour, and isn't the same edge we alread used...
            if tour[-1] in graph_inv[edge_idx_prime][1] and tour[-2] not in graph_inv[edge_idx_prime][1]:
                #...then add the vertex that hasn't been used yet.
                tour += [item for item in graph_inv[edge_idx_prime][1] if item != tour[-1]]
                break
        if tour[0] == tour[-1]:
            break
    return tour

def makeSortedTourEdgesFromSortedTourVertices(graph_inv, tour):
    """
    Inputs:
        graph_inv (dict): edge-wise description of the graph. See problem_def.py.
        tour (list): indices of graph vertices, listed in order visited to perform a vertex tour. Final value equals first value.
    Output:
        tour_edges (list): indices of edges, listed in order visited to perform a vertex tour.
        """
    tour_edges = []
    for g in range(len(graph_inv) - 1, -1, -1):
        if abs(tour.index(graph_inv[g][1][0]) - tour.index(graph_inv[g][1][1])) in (1, len(tour)-2):
            tour_edges.append(graph_inv[g])
    return tour_edges

def makeSortedTourVerticesFromSortedTourEdges(tour_edges):
    """
    Input:
        tour_edges (list): indices of edges, listed in order visited to perform a vertex tour.
    Output:
        tour (list): indices of graph vertices, listed in order visited to perform a vertex tour. Final value equals first value.
        """
    tour = [tour_edges[0][1][0], tour_edges[0][1][1]]
    for i in range(len(tour_edges)):
        for edge in tour_edges:
            if tour[-1] in edge[1] and tour[-2] not in edge[1]:
                tour += [item for item in edge[1] if item != tour[-1]]
                break
    return tour

def makeVertexSetFromEdgeVector(graph, x_edge_vector):
    vertex_set = []
    for edge_idx, x_value in enumerate(x_edge_vector):
        if x_value==1:
            vertex_set += [graph.graph_inv[edge_idx][1][0], graph.graph_inv[edge_idx][1][1]]
    vertex_set = list(set(vertex_set))
    return vertex_set
