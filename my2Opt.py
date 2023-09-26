import my_utils
import numpy as np

def my2Opt(graph, tour, UB, limit):

    tour_edges = []

    trials = 0

    for g in range(len(graph) - 1, -1, -1):
        if abs(tour.index(graph[g][1][0]) - tour.index(graph[g][1][1])) in (1, len(tour)-2):
            tour_edges.append(graph[g])

    reset_after_improvement = True
    while trials <= limit and reset_after_improvement:
        reset_after_improvement = False
        print('2 Opt tries again')
        for i in range(len(tour_edges)):
            for j in range(i+1, len(tour_edges)):
                trials = trials + 1
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
                butterfly1 = [item for item in graph if (abcd[0], abcd[2]) in item or (abcd[2], abcd[0]) in item]
                butterfly2 = [item for item in graph if (abcd[1], abcd[3]) in item or (abcd[1], abcd[3]) in item]

                if butterfly1 and butterfly2:
                    #must be an improvement.
                    if butterfly1[0][0] + butterfly2[0][0] < tour_edges[i][0] + tour_edges[j][0]:
                        print('remove: ', tour_edges[i], tour_edges[j])
                        print('add:    ', butterfly1[0], butterfly2[0])

                        if j < i:
                            del tour_edges[i]
                            del tour_edges[j]
                        else:
                            del tour_edges[j]
                            del tour_edges[i]

                        tour_edges.extend([butterfly1[0], butterfly2[0]])
                        tour = makeTour(tour_edges)
                        reset_after_improvement = True

    weight = sum([item[0] for item in tour_edges])
    return (tour, weight)

def makeTour(tour_edges):
    tour = []

    tour.extend([tour_edges[0][1][0], tour_edges[0][1][1]])

    uhoh_counter = 0

    while len(tour) < len(tour_edges)+1:
        for i in range(len(tour_edges)):
            if tour_edges[i][1][0] == tour[len(tour)-1]:
                if tour_edges[i][1][1] not in tour:
                    tour.append(tour_edges[i][1][1])
                if len(tour) == len(tour_edges) and tour_edges[i][1][1] == tour[0]:
                    tour.append(tour_edges[i][1][1])
            if tour_edges[i][1][1] == tour[len(tour) - 1]:
                if tour_edges[i][1][0] not in tour:
                    tour.append(tour_edges[i][1][0])
                if len(tour) == len(tour_edges) and tour_edges[i][1][0] == tour[0]:
                    tour.append(tour_edges[i][1][0])

        uhoh_counter+=1
        if uhoh_counter == 800:
            print 'Trouble in while loop of my2opt.maketour'
            print 'tour edges: ',tour_edges
            print 'tour so far: ',tour


    return tour