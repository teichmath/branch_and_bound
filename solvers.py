from bnb_driver import bnb_driver
from cutting_planes import find_cutting_planes_from_cliques_violating_cut_edge_set, \
                            find_cutting_planes_from_low_graph_cuts, \
                                compute_degree_2_constraints, \
                                    get_global_graph_cut_upper_bound
from heuristics import nearest_neighbor_vertex_tour, vertex_tour_2_opt, \
                                                    get_random_graph_cut, \
                                                    find_a_heavier_graph_cut
from subgraph_ops import check_if_subgraph_is_tour, check_if_edge_set_is_cut, \
                get_final_tour_vertices, get_final_vertices_with_optimal_cut_for_delta
                            
def tsp_solver(graph_name):
    verbosity = 1
    search_config = { 
        'problem name': 'Traveling Salesman Problem',
        'graph name': graph_name,
        'optimization type': 'minimize',
        'valid edge set type': 'vertex tour',
        'edge set validator': check_if_subgraph_is_tour, 
        'greedy heuristic': nearest_neighbor_vertex_tour,
        'improvement heuristic': vertex_tour_2_opt,
        'initial boundary finder': compute_degree_2_constraints,
        'cutting plane finder': find_cutting_planes_from_low_graph_cuts,
        'vertex set description': "Optimal tour vertex sequence: ",
        'vertex set getter': get_final_tour_vertices
    }
    solution = bnb_driver(graph_name, search_config, verbosity)

def max_cut_solver(graph_name):
    verbosity = 1
    search_config = {
        'problem name': 'Max Cut Problem',
        'graph name': graph_name,
        'optimization type': 'maximize',
        'valid edge set type': 'graph cut',
        'edge set validator': check_if_edge_set_is_cut, 
        'greedy heuristic': get_random_graph_cut,
        'improvement heuristic': find_a_heavier_graph_cut,
        'initial boundary finder': get_global_graph_cut_upper_bound,
        'cutting plane finder': find_cutting_planes_from_cliques_violating_cut_edge_set,
        'vertex set description': "Optimal cut is delta of vertex set:  ",
        'vertex set getter': get_final_vertices_with_optimal_cut_for_delta
    }
    solution = bnb_driver(graph_name, search_config, verbosity)





