import my_utils
import numpy as np
import random


def ezComb(g_inv, x, number_of_nodes):

    halfset = []  # ((edge_wt, (edge_node_1, edge_node_2)), index of edge)

    for i in range(1,len(x)):
        if abs(x[i] - 0.5) < 1e-6:
            halfset.append((g_inv[i], i))

    comb_attempts = 0

    if halfset:

        while comb_attempts < 5:

            comb_attempts+=1

            treelist = []  # (node, parent, rank)
            fertility = np.full(number_of_nodes + 1, True, dtype=bool)
            handle_edge_indices = []
            handle_nodes = []
            handle_weight = 0
            handle_achieved = False
            start_node = random.choice(halfset)[0][1][0]
            treelist.append((start_node, -1, 0))

            while True:

                #select hot node; for first run through, it will be the start node
                hot_node = [(-1)]
                arms_of_cycle = [[],[]]
                for node in treelist:
                    if fertility[node[0]]:
                        hot_node = node
                        break

                if hot_node[0] == -1:
                    break

                for edge in halfset:
                    #find 1/2 edge with hot node.
                    if hot_node[0] in edge[0][1]:
                        new_child_node = edge[0][1][0]
                        if hot_node[0] == edge[0][1][0]:
                            new_child_node = edge[0][1][1]
                        if new_child_node == hot_node[1]:
                            continue
                        # check for odd/even overlap
                        existing_mentions = [item for item in treelist if item[0] == new_child_node]
                        treelist.append((new_child_node, hot_node[0], hot_node[2]+1))
                        if existing_mentions:
                            oddeven = ' '
                            if (hot_node[2]+1)%2==0:
                                oddeven += 'even'
                            else:
                                oddeven += 'odd'
                            for mention in existing_mentions:
                                if (mention[2]%2 == 0 and 'odd' in oddeven) or (mention[2]%2 != 0 and 'even' in oddeven):
                                    arms_of_cycle[0].append((new_child_node, hot_node[0], hot_node[2]+1))
                                    arms_of_cycle[0].append(hot_node)
                                    arms_of_cycle[1].append(((new_child_node, mention[1], mention[2])))
                                    break
                            if arms_of_cycle[0]:
                                break


                fertility[hot_node[0]] = False

                if arms_of_cycle[0]:
                    rank_to_find = 1
                    while rank_to_find > 0:
                        parent_to_find = arms_of_cycle[0][-1][1]
                        rank_to_find = arms_of_cycle[0][-1][2] - 1
                        parent_set = [item for item in treelist if item[0] == parent_to_find and item[2] == rank_to_find]
                        if parent_set:
                            arms_of_cycle[0].append(parent_set[0])
                    rank_to_find = 1
                    connection_set = []
                    while rank_to_find > 0:
                        parent_to_find = arms_of_cycle[1][-1][1]
                        rank_to_find = arms_of_cycle[1][-1][2] - 1
                        parent_set = [item for item in treelist if item[0] == parent_to_find and item[2] == rank_to_find]
                        if parent_set:
                            arms_of_cycle[1].append(parent_set[0])
                        connection_set = [item for item in arms_of_cycle[0] if item[0] == arms_of_cycle[1][-1][0]]
                        if connection_set:
                            break
                    if connection_set:
                        for node_info in arms_of_cycle[1]:
                            if node_info[0] == connection_set[0][0]:
                                break
                            else:
                                if node_info[0] not in handle_nodes:
                                    handle_nodes.append(node_info[0])
                                if node_info[1] >= 0 and node_info[1] not in handle_nodes:
                                    handle_nodes.append(node_info[1])
                                edge_lookup = [item for item in halfset if item[0][1][0] in [node_info[0], node_info[1]] and item[0][1][1] in [node_info[0], node_info[1]] ]
                                if edge_lookup:
                                    handle_edge_indices.append(edge_lookup[0][1])
                                    handle_weight += edge_lookup[0][0][0]
                        for i in range(0,len(arms_of_cycle[0])):
                            if arms_of_cycle[0][i][0] == connection_set[0][0]:
                                break
                            else:
                                if arms_of_cycle[0][i][0] not in handle_nodes:
                                    handle_nodes.append(arms_of_cycle[0][i][0])
                                if arms_of_cycle[0][i][1] >= 0 and arms_of_cycle[0][i][1] not in handle_nodes:
                                    handle_nodes.append(arms_of_cycle[0][i][1])

                                edge_lookup = [item for item in halfset if item[0][1][0] in [arms_of_cycle[0][i][0], arms_of_cycle[0][i][1]] and item[0][1][1] in [arms_of_cycle[0][i][0], arms_of_cycle[0][i][1]] ]
                                if edge_lookup:
                                    handle_edge_indices.append(edge_lookup[0][1])
                                    handle_weight += edge_lookup[0][0][0]

                        print 'got handle.'
                        handle_achieved = True
                        break

            #We have a loop. Now we need teeth!
            teeth = []
            if handle_achieved:
                one_wts = [i for i in range(len(x)) if x[i] >= 0.8]
                one_edges = [(g_inv[i], i, x[i]) for i in one_wts]
                potential_teeth = [item for item in one_edges if (item[0][1][0] in handle_nodes or item[0][1][1] in handle_nodes) and (item[0][1][0] not in handle_nodes or item[0][1][1] not in handle_nodes)]
                tooth_ends = []
                for node in handle_nodes:
                    teeth_for_this_node = [item for item in potential_teeth if node in item[0][1]]
                    if teeth_for_this_node:
                        print 'at 133, hello'
                        my_tooth = max(teeth_for_this_node, key=lambda x:x[2])
                        print my_tooth,' is my tooth'
                        if my_tooth[0][1][0] == node:
                            if my_tooth[0][1][1] not in tooth_ends:
                                tooth_ends.append(my_tooth[0][1][1])
                                teeth.append(my_tooth)
                        if my_tooth[0][1][1] == node:
                            if my_tooth[0][1][0] not in tooth_ends:
                                tooth_ends.append(my_tooth[0][1][0])
                                teeth.append(my_tooth)

                if len(teeth) >= 3:
                    if len(teeth)%2 == 0:
                        teeth.remove(min(teeth, key=lambda x:x[2]))
                else:
                    teeth = []

            if teeth:
                # constraint is x(handle) + x(teeth) <= |handle nodes| + |teeth nodes| - |teeth| - (k+1)
                tooth_edge_indices = [item[1] for item in teeth]
                k = (len(teeth) - 1)/2
                rhs = len(handle_nodes) + len(teeth) - k - 1
                if len(handle_edge_indices)*0.5 + len(tooth_edge_indices) > rhs:
                    edge_set_for_comb_inequality = list(set(handle_edge_indices) | set(tooth_edge_indices))
                    print 'comb found. handle edges, with values: ',handle_edge_indices,' | ',[x[i] for i in handle_edge_indices]
                    print '            tooth edges, with values:  ',tooth_edge_indices,'  | ',[x[i] for i in tooth_edge_indices]
                    print '            rhs:                       ',rhs
                    return edge_set_for_comb_inequality, rhs
    else:
        print 'no 0.5 values found in comb handle search!'


    return [], 0
