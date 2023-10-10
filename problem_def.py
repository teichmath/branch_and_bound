"""Here we define the problem space. These are classes to give state and behavior to our graph, our LPs, and our solutions."""
import os
import copy
import numpy as np	
import scipy
from heuristics import my_nn, my_2_opt

def get_single_data(fname):
	"""
	This function gets the entries in the file fname
	Input: 
		fname (string): path of the file, e.g. 'TSP/st70.txt'
	Output: 
		data (list): data[0] gives [#nodes #edges];  data[1:] gives [node1 node2 weight].
	"""
	print("Reading file: "+fname)
	data = []
	data_dir = os.path.dirname(__file__)
	f = open(os.path.join(data_dir, fname),"r")
	f1 = f.readlines()
	for line in f1:
		line = line.split()
		if line:
			line = [int(i) for i in line]
			data.append(line)
	return data

class Solution:
	"""
	A vertex tour for a particular graph. The tour can be considered a work in process, since it can be improved within the Solution instance. 
	"""
	def __init__(self, tsp_graph, obj_value=None, x_vector=None, name=None):
		"""We can initialize a Solution either with an existing x vector, or with simply a graph, in which case a fresh vertex tour is found."""
		self.tsp_graph = tsp_graph
		self.name = name
		self.next_edge = 0
		self.improved_this_round = False
		self.stuck_at_local_minimum = False
		if not obj_value or not x_vector:
			visited, upper_bound = my_nn(self.tsp_graph.graph)
			self.objective_value = upper_bound
			self.x_vector = [0 for i in range(self.tsp_graph.m)]
			for i in range(len(visited) - 1):
				index = self.tsp_graph.graph[visited[i]][visited[i+1]][1]
				self.x_vector[index] = 1
		else:
			self.objective_value = obj_value
			self.x_vector = x_vector
	def get_objective_value(self):
		return self.objective_value
	def get_x_vector(self):
		return self.x_vector
	def try_improvement(self, limit=1e3):
		"""Attempt an improvement to the instance's vertex tour by way of the 2-opt method."""
		x_new, objv_new, next_edge = my_2_opt(self.tsp_graph.graph_inv, self.x_vector, self.objective_value, limit, self.next_edge)
		objv_improved = False
		if objv_new < self.objective_value:
			objv_improved = True
			self.improved_this_round = True
			self.objective_value = objv_new
			self.x_vector = x_new
		if next_edge==0:
			if not self.improved_this_round:
				self.stuck_at_local_minimum=True
			self.improved_this_round = False
		self.next_edge = next_edge
		return objv_improved
	def get_name(self):
		return self.name
	def is_at_minimum(self):
		"""This will be True if we attempted an improvement and exhausted the 2-opt without finding one."""
		return self.stuck_at_local_minimum

class TSPGraph:
	"""
	Stores information we need about the graph.
	"""
	def __init__(self, data):
		self.graph, self.graph_inv, self.graph_inv_sb, self.optimal, self.n, self.m = self.make_graph(data)
		self.deg_2_constraints_dict = self.compute_degree_2_constraints()
 	
	def make_graph(self, data):
		"""
		This function makes a graph out of the input data. It gives two different representations of the graph.
		Input: 
			data: a list of lists. 
				data[0] gives [#nodes #edges]
				data[1:] gives [node1 node2 weight]
				data[-1] potentially gives the minimum weight (if provided)
		Output: 
			graph: dictionary (vertex-wise) of dictionaries (edge-wise):   
				{ vert_idx: {other_vert_idx:(edg_wt, edg_idx), ... , other_vert_idx:(edg_wt, edg_idx)},   
				... ,
				vert_idx: {other_vert_idx:(edg_wt, edg_idx), ... , other_vert_idx:(edg_wt, edg_idx)} }
			graph_inv: dictionary (edge-wise):
				{edg_idx: (edg_wt, (vert_idx, other_vert_idx)),  ... , edg_idx: (edg_wt, (vert_idx, other_vert_idx)) }
				edges are ordered first by initial vertex index, and then by weight, in descending order.
				The edge indices in graph_inv match the edge indices in graph. 
			graph_inv_sb: dictionary (edge-wise):
				{edg_idx: (edg_wt, (vert_idx, other_vert_idx)),  ... , edg_idx: (edg_wt, (vert_idx, other_vert_idx)) }
				edges are ordered only by weight, in ascending order.
				The edge indices here do not match those in graph or graph_inv. This is a different indexing of the edges. 
			optimal: int. the minimum weight provided by the input file.
			n: int. number of vertices.
			m: int. number of edges.
		"""
		graph = {}
		graph_inv = {}
		optimal = None
		edge_index = -1
		for item in data:
			if len(item) ==1:
				optimal = item[0] # minimum weight
			if len(item) == 2:
				n = item[0]  # number of nodes
				m = item[1]  # number of edges
			if len(item) ==3:
				u = item[0]  # first node
				if u not in graph:
					graph[u] = {}
				v = item[1] # second node
				if v not in graph:
					graph[v] = {}
				# symmetric TSP graph
				weight = item[2]
				graph[u][v] = weight, edge_index
				graph[v][u] = weight, edge_index
				graph_inv[edge_index] = weight, (u,v)
			edge_index += 1
		graph_inv_sb = copy.deepcopy(graph_inv)
		self.quicksort(graph_inv_sb)
		return graph, graph_inv, graph_inv_sb, optimal, n, m

	def quicksort(self, edgewts):
		self.quicksort_drive(edgewts, 0, len(edgewts) - 1)

	def quicksort_drive(self, somewts, start, end):
		if start >= end:
			return
		p = self.partition(somewts, start, end)
		self.quicksort_drive(somewts, start, p - 1)
		self.quicksort_drive(somewts, p + 1, end)

	def partition(self, somewts, start, end):
		b = a = start
		while a < end:
			if somewts[a][0] <= somewts[end][0]:
				#swap a and b
				somewts[b], somewts[a] = somewts[a], somewts[b]
				b += 1
			a += 1
		somewts[b], somewts[end] = somewts[end], somewts[b]
		return b

	def vector2graph(self, g_inv,x):
		"""
		Convert a vector of weights into a new graph.
		Input: 
			g_inv: a dictionary of dictionary. 
				The global graph given in the inverse form. For each edge index i, 
				graph[i] = (weight, (u,v)), where (u,v) corresponds to the edge(i).
			x: a list. 
				The indices i of x correspond to the indices of the edges in g, and x[i] gives the weight that
				we want to set up for edge(i) in the new graph.
		Output:
			new_graph: a dictionary of dictionary. The new graph according to x.
		"""
		new_graph ={}

		for i in range(len(x)):
			node1 = g_inv[i][1][0]
			if node1 not in new_graph:
				new_graph[node1] = {}
			node2 = g_inv[i][1][1]
			if node2 not in new_graph:
				new_graph[node2] = {}
			new_graph[node1][node2] = x[i], i
			new_graph[node2][node1] = x[i], i

		return new_graph
	
	def compute_degree_2_constraints(self):
		full_index_list_of_nonzero_coefficients = []
		nonzero_counts_as_we_go = [0]
		for u in range(0, self.n):
			#u is some vertex in the graph. We want all adjacent edges to have combined indicator value of 2.
			adj_edge_indices = [value[1] for key, value in self.graph[u].items()]
			#print(str(u)+" edges are "+str(adj_edge_indices))
			full_index_list_of_nonzero_coefficients += adj_edge_indices
			cumulative_edge_index_count = len(adj_edge_indices) + nonzero_counts_as_we_go[-1]
			nonzero_counts_as_we_go.append(cumulative_edge_index_count)
		lower = np.array([2]*self.n, dtype=np.double)
		upper = np.array([2.01]*self.n, dtype=np.double)
		num_nz = len(full_index_list_of_nonzero_coefficients)
		start = np.array(nonzero_counts_as_we_go[:-1])
		index = np.array(full_index_list_of_nonzero_coefficients)
		value = np.array([1]*num_nz, dtype=np.double)
		return {'lower':lower, 'upper':upper, 'num_nz':num_nz, 'start':start, 'index':index, 'value':value}

class TSPLP:
	def __init__(self, tsp_graph):
		self.tsp_graph = tsp_graph
		self.weights = [x[0] for x in self.tsp_graph.graph_inv.values()]
		self.initialize_tsp_lp()

	def initialize_tsp_lp(self):
		self.cost = np.array(self.weights, dtype=np.double)
		self.var_bounds = [(0,1) for var in range(self.tsp_graph.m)]
		self.domain_coeff_matrix = []
		self.constraint_ubs = []
		self.domain_coeff_matrix_equality = []
		self.constraint_equality_values = []
		self.res = None

	def solve_lp(self):
		if len(self.domain_coeff_matrix)==0:
			if len(self.domain_coeff_matrix_equality)==0:
				self.res = scipy.optimize.linprog(self.cost, bounds=self.var_bounds)
			else:
				self.res = scipy.optimize.linprog(self.cost, A_eq=np.array(self.domain_coeff_matrix_equality), \
									  b_eq=np.array(self.constraint_equality_values), bounds=self.var_bounds)
		else:
			if len(self.domain_coeff_matrix_equality)==0:
				self.res = scipy.optimize.linprog(self.cost, A_ub=np.array(self.domain_coeff_matrix), \
									  b_ub=np.array(self.constraint_ubs), bounds=self.var_bounds)
			else:
				self.res = scipy.optimize.linprog(self.cost, A_ub=np.array(self.domain_coeff_matrix), \
									  b_ub=np.array(self.constraint_ubs), \
										A_eq=np.array(self.domain_coeff_matrix_equality), \
											b_eq=np.array(self.constraint_equality_values), \
												bounds=self.var_bounds)
		self.feasibility_message()
		return self.res.fun, self.res.success
		
	def is_feasible(self):
		if not self.res:
			print("no solving attempt made, hence no result.")
			return False
		else:
			self.feasibility_message()
			return self.res.success

	def feasibility_message(self):
		if not self.res.success:
			if self.res.status != 2:
				print("solve_lp not successful; code "+str(self.res.status))
				print("0 : Optimization terminated successfully. 1 : Iteration or time limit reached. "+\
				"2 : Problem appears to be infeasible. 3 : Problem appears to be unbounded."+\
				"4 : The HiGHS solver ran into a problem.")

	def get_existing_solution(self):
		if not self.res:
			print("no solving attempt made, hence no result.")
			return None
		else:		
			return self.res.x

	def get_objective_value_from_vector(self, x_vector):
		if len(self.cost) == len(x_vector):
			return sum([a*b for a,b in zip(self.cost, x_vector)])
		else:
			return np.inf

	def apply_tree_node_constraints(self, edge_binary_values):
		for i_v, val in enumerate(edge_binary_values):
			self.var_bounds[i_v] = (val, val)
	
	def apply_comb_inequality_constraint(self, handle_edges, tooth_edges, rhs):
		coeff_row = [0 for i in range(self.tsp_graph.m)]
		for h_idx in handle_edges:
			coeff_row[int(h_idx)] = 0.5
		for t_idx in tooth_edges:
			coeff_row[int(t_idx)] = 1
		self.domain_coeff_matrix.append(coeff_row)
		self.constraint_ubs.append(np.double(rhs))

	def apply_global_min_cut_constraints(self, min_cut_rows):
		for row in min_cut_rows:
			self.domain_coeff_matrix.append([-1*item for item in row])
		self.constraint_ubs += [np.double(-2) for i in range(len(min_cut_rows))]

	def apply_degree_2_constraints(self, equality=False):
		self.apply_constraints_from_highs_dict_format(self.tsp_graph.deg_2_constraints_dict, equality)	
	
	def apply_cutting_plane_constraints(self, cutting_plane_dict, equality=False):
		self.apply_constraints_from_highs_dict_format(cutting_plane_dict, equality)

	def apply_constraints_from_highs_dict_format(self, c, equality):
		index_position = 0
		start_plus_inf = [int(item) for item in c['start']] + [np.inf]
		a_base_array = []
		a_full_array = []
		b_full_vector = []
		for i_nz in range(len(c['start'])):
			this_row = [0 for i in range(self.tsp_graph.m)]
			p = start_plus_inf[i_nz]
			if i_nz < len(c['start']) - 1:
				q = start_plus_inf[i_nz+1]
				for idx, val in zip(c['index'][p:q], c['value'][p:q]):
					this_row[idx]=val
			else:
				for idx, val in zip(c['index'][p:], c['value'][p:]):
					this_row[idx]=val
			a_base_array.append(this_row)
		for base_row, low_bound, hi_bound in zip(a_base_array, c['lower'], c['upper']):
			if not equality:
				if low_bound > -np.inf:
					a_full_array.append([-1*coeff for coeff in base_row])
					b_full_vector.append(-1*low_bound)
				if hi_bound < np.inf:
					a_full_array.append(base_row)
					b_full_vector.append(hi_bound)
			else:
				a_full_array.append(base_row)
				b_full_vector.append(low_bound)
		if not equality:
			self.domain_coeff_matrix += a_full_array
			self.constraint_ubs += b_full_vector
		else:
			self.domain_coeff_matrix_equality += a_full_array
			self.constraint_equality_values += b_full_vector