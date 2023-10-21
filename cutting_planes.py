"""Methods for bounding the feasible region."""
import math
import numpy as np
from subgraph_ops import subgraph_a_to_b_walk, edge_contractions_max_to_min

def compute_degree_2_constraints(graph):
    """
    bounds to require that every vertex is degree 2.
    Inputs:
        graph: object representing the graph of the problem. see problem_def.py.
    Outputs:
        a dictionary in the original HiGHs format for initializing LPs.        
    """
    full_index_list_of_nonzero_coefficients = []
    nonzero_counts_as_we_go = [0]
    for u in range(0, graph.n):
        #u is some vertex in the graph. We want all adjacent edges to have combined indicator value of 2.
        adj_edge_indices = [value[1] for key, value in graph.graph[u].items()]
        full_index_list_of_nonzero_coefficients += adj_edge_indices
        cumulative_edge_index_count = len(adj_edge_indices) + nonzero_counts_as_we_go[-1]
        nonzero_counts_as_we_go.append(cumulative_edge_index_count)
    lower = np.array([2]*graph.n, dtype=np.double)
    upper = np.array([2.01]*graph.n, dtype=np.double)
    num_nz = len(full_index_list_of_nonzero_coefficients)
    start = np.array(nonzero_counts_as_we_go[:-1])
    index = np.array(full_index_list_of_nonzero_coefficients)
    value = np.array([1]*num_nz, dtype=np.double)
    return {'lower':lower, 'upper':upper, 'num_nz':num_nz, 'start':start, 'index':index, 'value':value}

def find_cutting_planes_from_low_graph_cuts(graph, x_solution, stop_crit_new_planes = 6, stop_crit_iters = 50):
    """
    retrieves lower bound constraints on graph cuts, found greedily.
    Inputs:
        graph: object representing the graph of the problem. see problem_def.py.
        x_solution (list): an edge set indicator vector.
        stop_crit_new_planes (int): a stopping criterion.
        stop_crit_iters (int): a stopping criterion.
    Outputs:
        cutting_plane_dict (dict): a dictionary in the original HiGHs format for initializing LPs. 
            """
    cutting_plane_dict, _ = subgraph_a_to_b_walk(graph, x_solution, stop_crit_new_planes, \
                                                     stop_crit_iters, False)
    return cutting_plane_dict

def find_cutting_planes_from_cliques_violating_cut_edge_set(graph, x_solution):
    """
    retrieves upper bound constraints on clique edges, found greedily.
    Inputs:
        graph: object representing the graph of the problem. see problem_def.py.
        x_solution (list): an edge set indicator vector.
    Outputs:
        cutting_plane_dict (dict): a simple representation of the constraints found. 
    """
    _, cliques = edge_contractions_max_to_min(graph, x_solution)
    coeff_rows=[]
    rhs_values=[]
    for clique in cliques:
        if clique['weight'] > clique['limit']:
            coeff_row = [0 for i in range(graph.m)]
            for e_idx in clique['edges']:
                coeff_row[e_idx]=1
            coeff_rows.append(coeff_row)
            rhs_values.append(clique['limit'])
    odd_cycles = get_heavy_k_threes(graph, x_solution)
    for k_3 in odd_cycles:
        k_3_edges = [graph.graph[k_3[0]][k_3[1]][1], graph.graph[k_3[0]][k_3[2]][1], graph.graph[k_3[1]][k_3[2]][1]] 
        coeff_row = [0 for i in range(graph.m)]
        for e_idx in k_3_edges:
            coeff_row[e_idx]=1
        coeff_rows.append(coeff_row)
        rhs_values.append(np.double(2))
    #'equality' is the orientation of the bound. '1' means upper bound, '-1' would mean lower bound.
    cutting_plane_dict = {'rows':coeff_rows, 'equality':[1]*len(coeff_rows), 'rhs':rhs_values}
    return cutting_plane_dict

def get_heavy_k_threes(graph, x_solution, individual_edge_limit=0.7):
    """
    finds K3 cliques that are heavier than some given amount.
    Inputs:
        graph: object representing the graph of the problem. see problem_def.py.
        x_solution (list): an edge set indicator vector.
        individual_edge_limit (double): we will consider only edges valued at or greater than this number.
    Outputs:
        odd_cycles (list of lists):  each nested list has the indices of 3 vertices that form a heavy clique. 
    """
    odd_cycles = []
    heavy_following_neighbors = {}
    for v_1 in range(graph.n):
        if v_1 not in heavy_following_neighbors.keys():
            heavy_following_neighbors[v_1] = get_heavy_following_neighbors(graph, x_solution, v_1, individual_edge_limit)
        for v_2 in heavy_following_neighbors[v_1]:
            if v_2 not in heavy_following_neighbors.keys():
                heavy_following_neighbors[v_2] = get_heavy_following_neighbors(graph, x_solution, v_2, individual_edge_limit)
            for v_3 in heavy_following_neighbors[v_2]:
                 if v_3 in heavy_following_neighbors[v_1]:
                      odd_cycles.append([v_1, v_2, v_3])
    return odd_cycles

def get_heavy_following_neighbors(graph, x_solution, i_v, individual_edge_limit):
    hfn=[]
    for key in graph.graph[i_v].keys():
         if key > i_v and x_solution[key] >= individual_edge_limit:
              hfn.append(key)
    return hfn

def get_global_graph_cut_upper_bound(graph):
    """returns a global upper bound for a graph cut edge set."""
    if graph.n % 2 == 0:
        rhs = np.double((graph.n**2)/4)
    else:
        half_n = math.floor(graph.n/2)
        rhs = half_n*(half_n+1)
    coeff_row = [1 for i in range(graph.m)]
    cutting_plane_dict = {'rows':[coeff_row], 'equality':[1], 'rhs':[rhs]}
    return cutting_plane_dict


