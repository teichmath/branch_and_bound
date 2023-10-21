"""Greedy methods for finding special edge sets."""
import math
import random
import numpy as np
from subgraph_ops import makeSortedTourVerticesFromSolutionVector, \
                        makeSolutionVectorFromSortedTourEdges, \
                        makeSortedTourEdgesFromSortedTourVertices, \
                        try_a_cut_increase, get_subgraph_delta, \
                        makeVertexSetFromEdgeVector

def get_random_graph_cut(graph, k=None):
    """
    retrieves a graph cut found randomly.
    Inputs:
        graph: object representing the graph of the problem. see problem_def.py.
        k (int): size of vertex set to sample. 
    Outputs:
        edge_set (list): indices of edges in the cut.
        total_weight (double): total of edge weights in edge set.
    """
    if k is None:
        k = math.floor(graph.n/2)
    edge_set = get_random_cut_edge_set(graph, k)
    weights = [x[0] for x in graph.graph_inv.values()]
    total_weight = sum([e*w for e,w in zip(edge_set,weights)])
    return edge_set, total_weight

def get_random_cut_edge_set(graph, k):
    """
    gets random graph cut by randomly sampling vertices, then finding the delta of the vertex set.
    Inputs:
        graph: object representing the graph of the problem. see problem_def.py.
        k (int): size of vertex set to sample. 
    Outputs:
        edge_set (list): indices of edges in the cut.
    """
    vertex_A_set = random.sample([i for i in range(graph.n)], k)
    edges_in_cut = get_subgraph_delta(graph.graph, vertex_A_set)
    edge_set = [0 for i in range(graph.m)]
    for edge_idx in edges_in_cut:
        edge_set[edge_idx]=1
    return edge_set

def find_a_heavier_graph_cut(graph, x_vector, limit=100, start_vertex=0):
    """
    operates on a set of vertices, the delta of which is a given graph cut, to increase 
    the cut's weight.
    Inputs:
        graph: object representing the graph of the problem. see problem_def.py.
        x_vector (list):  edge set indicator vector for the graph cut.
        limit (int): a stopping criterion.
        start_vertex (int): the index of the vertex on which to begin.
    Outputs:
        x_edge_set (list): edge set indicator vector for the possibly improved graph cut.
        best_cut_weight (double): total of edge weights in the found graph cut.
        next_vertex (int): the index of the vertex on which to begin, if and when this 
        function is called again.
    """
    next_vertex = start_vertex
    x_edge_set = x_vector
    vertex_set = makeVertexSetFromEdgeVector(graph, x_edge_set)
    weights = [x[0] for x in graph.graph_inv.values()]
    best_cut_weight = sum([e*w for e, w in zip(x_edge_set, weights)])
    cut_search_iter = 0
    while True:
        vertex_set = try_a_cut_increase(graph, vertex_set, next_vertex)    
        edge_set = get_subgraph_delta(graph.graph, vertex_set)
        this_x_edge_set = [0 for i in range(graph.m)]
        for edge_idx in edge_set:
            this_x_edge_set[edge_idx] = 1
        this_cut_weight = sum([e*w for e, w in zip(this_x_edge_set, weights)])
        next_vertex+=1
        if next_vertex >= graph.n:
            next_vertex = 0
        if this_cut_weight <= best_cut_weight:
            break
        else:
            best_cut_weight = this_cut_weight
            x_edge_set = this_x_edge_set
        cut_search_iter += 1
        if cut_search_iter > limit:
            break
    return x_edge_set, best_cut_weight, next_vertex

def vertex_tour_2_opt(graph, x_solution, limit=100, start_edge=0):
    """
    2-Opt algorithm for searching for improvements for an existing vertex tour.
    Inputs:
        graph_inv (dict): edge-wise description of the graph. See problem_def.py.
        x_solution (list): edge set vector we begin the process with (should 
                            already be a vertex tour).
        limit (int): maximum iterations.
        start_edge (int): edge in vector at which to begin the process.
    Outputs:
        x_solution (list): edge set vector (a vertex tour) at the end of the process.
        weight (double): vertex tour weight at the end of the process.
        next_edge (int): if method were used again for this graph and with this 
                            solution, this would be the start edge to use to 
                            continue the process.    
    """
    graph_inv = graph.graph_inv
    tour = makeSortedTourVerticesFromSolutionVector(graph_inv, x_solution)
    tour_edges = makeSortedTourEdgesFromSortedTourVertices(graph_inv, tour)
    trials = 0
    reset_after_improvement = True
    end_edge = start_edge
    while trials <= limit and reset_after_improvement:
        reset_after_improvement = False
        #print('2 Opt tries again')
        for i in range(start_edge, len(tour_edges)):
            if trials>limit:
                break
            end_edge = i
            for j in range(i+1, len(tour_edges)):
                trials += 1
                #must butterfly, can't have any in common
                abcd = []
                abcd.extend([tour_edges[i][1][0], tour_edges[i][1][1], \
                             tour_edges[j][1][0], tour_edges[j][1][1]])
                if len(abcd) != len(set(abcd)):
                    continue
                abcd_sorted = False
                while not abcd_sorted:
                    abcd_sorted = True
                    for k in range(3):
                        if tour.index(abcd[k]) > tour.index(abcd[k+1]):
                            abcd_sorted = False
                            abcd[k], abcd[k+1] = abcd[k+1], abcd[k]
                butterfly1 = [item for item in graph_inv.values() if (abcd[0], abcd[2]) \
                              in item or (abcd[2], abcd[0]) in item]
                butterfly2 = [item for item in graph_inv.values() if (abcd[1], abcd[3]) \
                              in item or (abcd[1], abcd[3]) in item]
                if butterfly1 and butterfly2:
                    #must be an improvement.
                    if butterfly1[0][0] + butterfly2[0][0] < tour_edges[i][0] + tour_edges[j][0]:
                        #print('remove: ', tour_edges[i], tour_edges[j])
                        #print('add:    ', butterfly1[0], butterfly2[0])
                        if j < i:
                            del tour_edges[i]
                            del tour_edges[j]
                        else:
                            del tour_edges[j]
                            del tour_edges[i]
                        tour_edges.extend([butterfly1[0], butterfly2[0]])
                        sol_vector = makeSolutionVectorFromSortedTourEdges(graph_inv, tour_edges)
                        tour = makeSortedTourVerticesFromSolutionVector(graph_inv, sol_vector)
                        reset_after_improvement = True
    weight = sum([item[0] for item in tour_edges])
    x_solution = makeSolutionVectorFromSortedTourEdges(graph_inv, tour_edges)
    next_edge = end_edge+1
    if end_edge == len(tour_edges):
        next_edge = 0
    return x_solution, weight, next_edge

def nearest_neighbor_vertex_tour(graph):
    """
    nearest neighbor algorithm for greedily finding a vertex tour.
    Input:
        graph (dict): vertex-wise description of the graph. See problem_def.py.
    Outputs:
        x_final (list): edge set indicator vector for tour.
        total_weight (double): sum of the weights of the tour edges.
    """
    graph_dict = graph.graph
    n = graph.n # number of nodes
    m = graph.m
    nodes = range(n)
    while True:
        visited = []
        # randomly select a start node
        start_node = np.random.randint(n)
        total_weight = 0
        # start visiting. Only stop when all nodes are visited
        current_node = start_node
        visited.append(current_node)
        visit_counter = 0
        while set(visited) != set(nodes) and visit_counter < n+1:
            neighbors = graph_dict[current_node]
            # search for the nearest unvisited neighbor nn
            w = float("inf")
            for v in neighbors:
                if v not in visited:
                    if neighbors[v][0] < w:
                        nn = v
                        w = neighbors[nn][0]
            # update visited and total weight of the tour
            visited.append(nn)
            total_weight += w
            current_node = nn
            visit_counter+=1
        if set(visited) == set(nodes):
            # return to start node
            if start_node in graph_dict[visited[-1]]:
                visited.append(start_node)
                total_weight += graph_dict[nn][start_node][0]
                break
    x_final = np.zeros(m)
    for i in range(len(visited) - 1):
        index = graph_dict[visited[i]][visited[i+1]][1]
        x_final[index] = 1
    return x_final, total_weight