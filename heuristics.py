"""Greedy methods for finding vertex tours."""
import numpy as np

def my_2_opt(graph_inv, x_solution, limit=100, start_edge=0):
    """
    2-Opt algorithm for searching for improvements for an existing vertex tour.
    Inputs:
        graph_inv (dict): edge-wise description of the graph. See problem_def.py.
        x_solution (list): edge set vector we begin the process with (should already be a vertex tour).
        limit (int): maximum iterations.
        start_edge (int): edge in vector at which to begin the process.
    Outputs:
        x_solution (list): edge set vector (a vertex tour) at the end of the process.
        weight (double): vertex tour weight at the end of the process.
        next_edge (int): if method were used again for this graph and with this solution, this would be the start edge to use to continue the process.    
    """
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
                abcd.extend([tour_edges[i][1][0], tour_edges[i][1][1], tour_edges[j][1][0], tour_edges[j][1][1]])

                if len(abcd) != len(set(abcd)):
                    continue

                abcd_sorted = False
                while not abcd_sorted:
                    abcd_sorted = True
                    for k in range(3):
                        if tour.index(abcd[k]) > tour.index(abcd[k+1]):
                            abcd_sorted = False
                            abcd[k], abcd[k+1] = abcd[k+1], abcd[k]
                butterfly1 = [item for item in graph_inv.values() if (abcd[0], abcd[2]) in item or (abcd[2], abcd[0]) in item]
                butterfly2 = [item for item in graph_inv.values() if (abcd[1], abcd[3]) in item or (abcd[1], abcd[3]) in item]

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
                        tour = makeSortedTourVerticesFromSortedTourEdges(tour_edges)
                        reset_after_improvement = True
    weight = sum([item[0] for item in tour_edges])
    x_solution = makeSolutionVectorFromSortedTourEdges(graph_inv, tour_edges)
    next_edge = end_edge+1
    if end_edge == len(tour_edges):
        next_edge = 0
    return x_solution, weight, next_edge

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
    edge_indices = [i for i in range(len(x_solution)) if x_solution[i]==1]
    tour = [graph_inv[edge_indices[0]][1][0], graph_inv[edge_indices[0]][1][1]]
    for i_edge_idx in range(1, len(edge_indices)):
        for edge_idx_prime in edge_indices:
            if tour[-1] in graph_inv[edge_idx_prime][1] and tour[-2] not in graph_inv[edge_idx_prime][1]:
                tour += [item for item in graph_inv[edge_idx_prime][1] if item != tour[-1]]
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

def my_nn(graph):
    """
    Nearest Neighbor algorithm for greedily finding a vertex tour.
    Input:
        graph (dict): vertex-wise description of the graph. See problem_def.py.
    Outputs:
        visited (list): 
    """
    n = len(graph.keys()) # number of nodes
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
            neighbors = graph[current_node]
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
            if start_node in graph[visited[-1]]:
                visited.append(start_node)
                total_weight += graph[nn][start_node][0]
                break

    return visited, total_weight