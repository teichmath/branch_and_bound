import my_utils
import numpy as np
import copy

# ww17
# Last updated 04/15/2018

def myCF(graph):
	"""
	Implementation of Christofides' Heuristic for Upper Bound.
	Input:
		graoh: a dictionary of dictionary. 
			graph[u][v][0] = weight of (u,v)
			graph[u][v][1] = edge index of (u,v)
	Output:
		visited: a list.
			A list of visited nodes that form the tour.
		total_weight: intger. 
			Total weight of the tour.
	"""
	#graph = my_utils.make_graph(data)[0]
	n = len(graph.keys()) # number of nodes
	#m = data[0][1] # number of edges

	#nodes = range(n)
	#visited = []

	# get a min cost spanning tree with Kruskal's.
	MSTedges = kruskalsMST(n, g_inv)

	#Find W, the set of nodes which have an odd degree in the graph induced by MSTedges.
	W = []
	degrees = [0 for i in range(n+1)]
	for i in range(n+1):
		for j in (i+1,n+1):
			if graph[i][j]:
				degrees[i]+=1
				degrees[j]+=1
	#Get the subgraph of G induced by W.
	G_w = []
	for edge in g_inv:
		if degrees[edge[1][0]] % 2 != 0 and degrees[edge[1][1]] % 2 != 0:
			if(edge[1][0] not in W):
				W.append(edge[1][0])
			if(edge[1][1] not in W):
				W.append(edge[1][1])
			G_w.append(edge)
	# Find a perfect matching of G_w with the Blossom Algorithm.


def kruskalsMST(n, g_inv):
	# list block numbers by node index
	node2block = [i for i in range(n)]
	# list block cardinalities with block number as index
	block2card = [1 for i in range(n)]
	# initialize edge set of MST
	MSTedges = []
	#go through edges
	for edge in g_inv:
		#are the nodes in separate blocks
		node1 = edge[1][0]
		node2 = edge[1][1]
		if node2block[node1] != node2block[node2]:
			MSTedges.append(edge)
			node_host, node_guest = node1, node2
			if block2card[node2block[node1]]  < block2card[node2block[node2]]:
				node_host, node_guest = node2, node1
			block2card[node2block[node_host]] += block2card[node2block[node_guest]]
			block2card[node2block[node_guest]] = 0
			guestblock = node2block[node_guest]
			hostblock = node2block[node_host]
			for j in range(n):
				if node2block[j] == guestblock:
					node2block[j] = hostblock
	return MSTedges


def blossomForPerfectMatching(G_w, n):
	#initialize
	M = []  # matching on the main graph (G_w); elements are (weight, (u,v))
	M_prime = []  # matching on the graph with blossoms shrunk; elements are (weight, (u,v))
	#G_prime_all = []  # edges of G; elements are (u,v)
	nodes_in_G_prime = [] #list of node names
	for edge in G_w:
	#	G_prime_all.append(edge[1])
		if edge[1][0] not in nodes_in_G_prime:
			nodes_in_G_prime.append(edge[1][0])
		if edge[1][1] not in nodes_in_G_prime:
			nodes_in_G_prime.append(edge[1][1])
#	G_prime = copy.deepcopy(G_prime_all) #edges of mutable auxiliary graph; elements are (u,v)
	G_prime = copy.deepcopy(G_w) #now full edges!
	r = 0 # root node
	new_high_index = n + 1
	blossoms = [] # we keep track of contracted odd cycles here
	any_nodes_exposed = True

	status = 'begin tree'

	#Get an M prime exposed node. if none exists, we have a perfect matching.
	while any_nodes_exposed:

		any_nodes_exposed = False

		for node in nodes_in_G_prime:
			# if node is not mentioned in M_prime, it's M prime exposed
			if not [item for item in M_prime if node in item[1]]:
				r = node
				any_nodes_exposed = True
				break

		if not any_nodes_exposed:
			return M_prime

		#Initialize tree

		print 'r :', r

		T = [] # collection of tree edges; elements are (u,v)
		AB_T = [(r,)] # list of tuples representing nodes at index distance from root node. elements are integers. tuples at
					# even indices are class B, at odd indices are class A.

		found_an_edge = True
		status = 'no edge'

		while found_an_edge:

			print 'T: ',T,' AB_T: ',AB_T,' M_prime: ',M_prime

			#Get an edge vw with v in B_T and w not in A_T. If none exists, G_w has no perfect matching (because we still have
			# an exposed node).
			found_an_edge = False
			v = 0
			w = 0
			#find an available v
			for i in range(0, len(AB_T)):
				if i % 2 == 0:
					for node in AB_T[i]:
						#make set of edges with v, drawing from G_prime
						these_edges = [item for item in G_prime if node in item[1]]
						print 'these edges: ',these_edges
						for edge in these_edges:
							status = 'no edge'
							node_to_check = edge[1][0]
							if edge[1][0] == node:
								node_to_check = edge[1][1]
							# for each, check if other node is in AB_T (in a B group)
							print 'node to check: ',node_to_check
							for j in range(len(AB_T)):
								if node_to_check in AB_T[j]:
									status = 'w in A'
									if j%2==0:
										#if so, update status and break
										status = 'w in B'
										v = node
										w = node_to_check
										found_an_edge = True
										case_w_in_B(M_prime, T, G_prime, v, w, AB_T, new_high_index)
									break
							if status == 'no edge':
								#if not, see if it's mentioned in M_prime
								v = node
								w = node_to_check
								if [item for item in M_prime if node_to_check in item[1]]:
									status = 'w not in T, is M_prime covered'
									case_w_Covered(T, v, w, AB_T, M_prime)
								else:
									status = 'w not in T, is M_prime exposed'
									case_w_Exposed(G_w, G_prime, M_prime, T, r, v, w, AB_T, blossoms)
								found_an_edge = True
							if status in ('w in B', 'w not in T, is M_prime exposed'):
								break
						if status in ('w in B', 'w not in T, is M_prime exposed'):
							break
				if status in ('w in B', 'w not in T, is M_prime exposed'):
					break
			if status == 'w not in T, is M_prime exposed':
				break

		if found_an_edge == False:
			print 'No perfect matching!'
			return


def case_w_in_B(M_prime, T, G_prime, v, w, AB_T, new_high_index):
	shrink(M_prime, T, G_prime, v, w, AB_T, new_high_index)
	new_high_index += 1
	if r == v or r == w:
		r = new_high_index

def case_w_Covered(T, v, w, AB_T, M_prime):
	# There must exist an edge wz in M_prime such that z is not in T either.
	full_edge_with_z = next(item for item in M_prime if w in item[1])
	z = full_edge_with_z[1][0]
	if (full_edge_with_z[1][0] == w):
		z = full_edge_with_z[1][1]

	extendT(T, v, w, AB_T)
	extendT(T, w, z, AB_T)
	print w,'covered;  T: ', T, ' AB_T: ', AB_T


def case_w_Exposed(G_w, G_prime, M_prime, T, r, v, w, AB_T, blossoms):

	augmentMPrime(G_prime, M_prime, T, r, v, w, AB_T)

	need_a_G_copy = False

	if blossoms:
		need_a_G_copy = True

	# Extend M_prime to a matching of G
	# iterate through contracted blossoms; expand each
	for blossom in blossoms:
		# which edge and node are connected to contracted blossom by M_prime? (there exists only one)
		edge_to_contracted_blossom = next(item for item in M_prime if blossom[0] in item[1])
		node_to_contracted_blossom = edge_to_contracted_blossom[0][1][0]
		if edge_to_contracted_blossom[0][1][0] == blossom[0]:
			node_to_contracted_blossom = edge_to_contracted_blossom[0][1][1]
		# which node in the blossom is adjacent to that? (there might be more than one, we choose one)
		edges_we_could_include = [item for item in G_w if node_to_contracted_blossom in item[1]]
		node_we_cover_first_in_blossom = edges_we_could_include[0][1][0]
		if edges_we_could_include[0][1][0] == node_to_contracted_blossom:
			node_we_cover_first_in_blossom = edges_we_could_include[0][1][1]
		# Update M_prime
		M_prime.remove(edge_to_contracted_blossom[0])
		M_prime.append(edges_we_could_include[0])
		# get blossom edges
		edges_in_blossom = [item for item in G_w if item[1][0] in blossom and item[1][1] in blossom]
		nodes_used_in_blossom = (node_we_cover_first_in_blossom)
		accept_edges = false
		# find next edges to consider
		blossom_edges_remain = true;
		while blossom_edges_remain:
			my_next_blossom_edges = [item for item in edges_in_blossom if
									 len(intersection(item[1], nodes_used_in_blossom)) == 1]
			if my_next_blossom_edges == false:
				blossom_edges_remain = false
			if (accept_edges):
				for edge in my_next_blossom_edges:
					M_prime.append(edge)
					nodes_used_in_blossom.append(edge[1][0])
					nodes_used_in_blossom.append(edge[1][1])
			accept_edges ^= true

	# Reset G_prime
	if need_a_G_copy:
		G_prime = copy.deepcopy(G_w)



def extendT(T, v, w, AB_T):
	T.append((v, w))
	for j in range(len(AB_T)-1, -1, -1):
		if v in AB_T[j]:
			if j < len(AB_T)-1:
				AB_T[j+1] += (w,)
			else:
				AB_T.append((w,))

def removeEdgeFromT(T, edge, AB_T):
	T.remove(edge)
	#get indices of edge[0] and edge[1]; if 0<1, find 1 in AB_T and remove it; else find 0 and remove.

def augmentMPrime(g_inv, M_prime, T, r, v, w, AB_T):

	extendT(T, v, w, AB_T)

	P = getPath(r, w, T, AB_T)

	print 'in augment: P is ',P

	for edge in P:
		contribute = True
		existing_edge = [item for item in M_prime if edge in item]
		if(existing_edge):
			print 'existing_edge: ',existing_edge
			contribute = False
			M_prime.remove(existing_edge[0])
		if(contribute):
			print 'the edge is ',edge
			edge_to_add = [item for item in g_inv if edge in item or (edge[1], edge[0]) in item]
			print 'edge_to_add: ',edge_to_add
			M_prime.append(edge_to_add[0])

def getPath(back_node, forward_node, T, AB_T):
	P = []
	path_started = False
	current_node = forward_node
	for node_group in reversed(AB_T):
		if current_node == back_node:
			break
		for node in node_group:
			if path_started:
				if (node, current_node) in T or (current_node, node) in T:
					if(node, current_node) in T:
						P.append((node, current_node))
					else:
						P.append((current_node, node))
					current_node = node
					break
			else:
				if node == forward_node:
					path_started = True
					break
	return P

def shrink(M_prime, T, G_prime, v, w, AB_T, new_high_index):
	P = getPath(v, w, T, AB_T)
	blossom = [(new_high_index, ())]
	for edge in P:
		for indx in (0, 1):
			if edge[indx] not in blossom:
				blossom[0][1].append(edge[indx])
				for node_group in AB_T:
					if edge[indx] in node_group:
						node_group.remove(edge[indx])
						break

	for edge in G_prime:
		edge_status = 'keep'
		if edge[1][0] in blossom[1] and edge[1][1] in blossom[1]:
			edge_status = 'remove'
		elif edge[1][0] in blossom[1]:
			edge_status = 'rename 0'
		elif edge[1][1] in blossom[1]:
			edge_status = 'rename 1'

		if edge_status == 'rename 0':
			edge[1][0] = new_high_index
		if edge_status == 'rename 1':
			edge[1][1] = new_high_index
		if edge_status == 'remove':
			if edge in M_prime:
				M_prime.remove(edge)
			if edge[1] in T:
				removeEdgeFromT(edge[1])
			G_prime.remove(edge)

	return blossom
