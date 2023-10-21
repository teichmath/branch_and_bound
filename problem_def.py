"""Here we define the problem space. These are classes to give state and behavior to our graph,
	our LPs, and our solutions."""
import os
import copy
import numpy as np
import scipy

def get_single_data(fname):
	"""
	gets the entries in the file fname.
	Input: 
		fname (string): path of the file, e.g. 'TSP/st70.txt'
	Output: 
		data (list): data[0] gives [#nodes #edges];  data[1:] gives [node1 node2 weight].
	"""
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
	a vertex tour for a particular graph. the tour can be considered a work in progress, 
	since it can be improved within the Solution instance. 
	"""
	def __init__(self, graph, obj_value=None, x_vector=None, name=None, heuristic=None, minimize=True):
		"""
		we can initialize a Solution either with an existing x vector, or with simply 
		a graph, in which case a fresh valid solution is generated here.
		Inputs:
        	graph: object representing the graph of the problem. see below.
			obj_value (double): objective value for an existing solution.
			x_vector (list): edge set indicator vector for an existing solution.
			name (string): accessed primarily for debugging.
			heuristic (method): for finding a fresh valid solution.
				Input:  a graph object.
				Outputs:  (list) an edge set indicator vector;
						  (double) the objective value associated with the list.
			minimize (boolean): indicates whether the problem is to minimize or 
								maximize the objective. 
		"""
		self.graph = graph
		self.name = name
		self.improvement_bookmark = 0
		self.improved_this_round = False
		self.stuck_at_local_minimum = False
		self.minimize = minimize
		if (obj_value is None or x_vector is None) and heuristic:
			self.x_vector, self.objective_value = heuristic(self.graph)
			if not self.minimize:
				self.objective_value *= -1
		else:
			self.objective_value = obj_value
			self.x_vector = x_vector
			if obj_value is None:
				self.objective_value = np.inf
	def get_objective_value(self):
		return self.objective_value
	def get_x_vector(self):
		return self.x_vector
	def try_improvement(self, improvement_heuristic, limit=1e3):
		"""
		attempt an improvement to the instance's vertex tour by way of the 2-opt method.
		Inputs:
			improvement_heuristic (method): for finding an improved valid solution.
				Input:  a graph object;
						(list) an edge set indicator vector;
						(int) a stopping criterion;
						(int) a bookmark index, for either a vertex or an edge, depending 
							on the formulation of the heuristic. if the process is stopped 
							before every vertex or edge is considered, this is a way to 
							keep track of where it left off, in case it makes sense to 
							pick up in the same place when revisiting the solution.
				Outputs:  (list) an edge set indicator vector;
						  (double) the objective value associated with the list;
						  (int) a bookmark index, as described above.
			limit (int): a stopping criterion.
		Outputs:
			objv_improved (boolean): reports whether an improvement was found.
		"""
		x_new, objv_new, improvement_bookmark = improvement_heuristic(self.graph, self.x_vector, \
																limit, self.improvement_bookmark)
		if not self.minimize:
			objv_new *= -1
		objv_improved = False
		if objv_new < self.objective_value:
			objv_improved = True
			self.improved_this_round = True
			self.objective_value = objv_new
			self.x_vector = x_new
		#if the bookmark is 0, and we didn't just see an improvement, we
		# suppose that this solution is stuck at a local minimum. This
		# is not generally true. In the future it would make more sense
		# to give up on this solution if the bookmark has made any full
		# rotation without an improvement.
		if improvement_bookmark==0:
			if not self.improved_this_round:
				self.stuck_at_local_minimum=True
			self.improved_this_round = False
		self.improvement_bookmark = improvement_bookmark
		return objv_improved
	def get_name(self):
		return self.name
	def is_at_minimum(self):
		"""this will be True if we attempted an improvement and exhausted the 2-opt without \
			finding one."""
		return self.stuck_at_local_minimum

class BnBGraph:
	"""  stores information we need about the graph. """
	def __init__(self, data):
		self.graph, self.graph_inv, self.graph_inv_sb, self.optimal, self.n, self.m = \
																	self.make_graph(data)
	def make_graph(self, data):
		"""
		makes a graph out of the input data. it gives two different 
		representations of the graph.
		Input: 
			data: a list of lists. 
				data[0] gives [#nodes #edges]
				data[1:] gives [node1 node2 weight]
				data[-1] potentially gives the minimum weight (if provided)
		Output: 
			graph: dictionary (vertex-wise) of dictionaries (edge-wise):   
				{ vert_idx: {other_vert_idx:(edg_wt, edg_idx), ... },   
				... ,
				vert_idx: {other_vert_idx:(edg_wt, edg_idx), ... )} }
			graph_inv: dictionary (edge-wise):
				{edg_idx: (edg_wt, (vert_idx, other_vert_idx)),  
											... , 
												edg_idx: (edg_wt, (vert_idx, other_vert_idx)) }
				edges are ordered first by initial vertex index, and then by weight, 
				in descending order. The edge indices in graph_inv match the edge indices 
				in graph. 
			graph_inv_sb: dictionary (edge-wise):
				{edg_idx: (edg_wt, (vert_idx, other_vert_idx)),  
											... , 
												edg_idx: (edg_wt, (vert_idx, other_vert_idx)) }
				edges are ordered only by weight, in ascending order.
				The edge indices here do not match those in graph or graph_inv. 
				This is a different indexing of the edges. 
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

class LPHandler:
	""" handles our LP solver and stores cutting planes. """
	def __init__(self, graph, initial_bound_finder, vprint, minimize=True):
		self.graph = graph
		self.weights = [x[0] for x in self.graph.graph_inv.values()]
		if not minimize:
			self.weights = [-1*x[0] for x in self.graph.graph_inv.values()]
		self.cutting_plane_highs_dict_binder = []
		self.cutting_plane_row_rhs_dict_binder = []
		self.initial_cutting_planes = initial_bound_finder(graph)
		self.sort_new_planes_dict(self.initial_cutting_planes)
		self.initialize_lp()
		self.vprint = vprint

	def initialize_lp(self):
		self.cost = np.array(self.weights, dtype=np.double)
		self.var_bounds = [(0,1) for var in range(self.graph.m)]
		self.domain_coeff_matrix = []
		self.constraint_ubs = []
		self.domain_coeff_matrix_equality = []
		self.constraint_equality_values = []
		self.res = None

	def solve_lp(self):
		""" 
		uses scipy's interior point LP solver, which is actually the HiGHs solver.
		a commercial solver (such as Gurobi's) used here could speed things up.
		we store the LP result as self.res.
		note that the LP solved here is only a subproblem of our main problem, and 
		it is a relaxation at that. 
		Outputs:
			self.res.fun (double): objective value of solution.
			self.res.success (boolean): indicates whether the problem was feasible.
						note that our objective function is simply an inner product, 
						and this is an LP relaxation; therefore "not feasible" 
						would mean an empty feasible set.
		"""
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
			self.vprint("no solving attempt made, hence no result.", 0)
			return False
		else:
			self.feasibility_message()
			return self.res.success

	def feasibility_message(self):
		if not self.res.success:
			if self.res.status != 2:
				self.vprint("solve_lp not successful; code "+str(self.res.status), 0)
				self.vprint("0 : Optimization terminated successfully. 1 : Iteration or time limit reached. "+\
				"2 : Problem appears to be infeasible. 3 : Problem appears to be unbounded."+\
				"4 : The HiGHS solver ran into a problem.", 0)

	def get_existing_solution(self):
		if not self.res:
			self.vprint("no solving attempt made, hence no result.", 0)
			return None
		else:		
			return self.res.x

	def get_proximal_integral_result_from_nonintegral(self, x_solution):
		integral_x_solution = np.heaviside(x_solution-1e-13, [1]*len(x_solution))
		integral_objective_value = self.get_objective_value_from_vector(integral_x_solution)
		return integral_x_solution, integral_objective_value

	def get_objective_value_from_vector(self, x_vector):
		if len(self.cost) == len(x_vector):
			return sum([a*b for a,b in zip(self.cost, x_vector)])
		else:
			return np.inf

	def define_domain(self, sub_lp):
		"""
		policy for applying constraints to our subproblems.
		Inputs:
			sub_lp (list): item at index i is either 0 or 1, giving the rhs of an equality constraint for edge x_i.
		"""
		self.vprint("applying tree node constraints.", 2, end= '  ')
		self.apply_tree_node_constraints(sub_lp)
		self.apply_cutting_plane_constraints()
	
	def find_new_cutting_planes(self, cutting_plane_finder, x_solution):
		cutting_plane_dict = cutting_plane_finder(self.graph, x_solution)
		self.sort_new_planes_dict(cutting_plane_dict)

	def sort_new_planes_dict(self, new_dict):
		if 'rhs' in new_dict.keys():
			self.cutting_plane_row_rhs_dict_binder.append(new_dict)
		elif 'lower' in new_dict.keys():
			self.cutting_plane_highs_dict_binder.append(new_dict)
	
	def apply_tree_node_constraints(self, edge_binary_values):
		for i_v, val in enumerate(edge_binary_values):
			self.var_bounds[i_v] = (val, val)
	
	def apply_cutting_plane_constraints(self):
		for highs_constraint_group in self.cutting_plane_highs_dict_binder:
			self.apply_constraints_from_highs_dict_format(highs_constraint_group)
		for row_rhs_constraint_group in self.cutting_plane_row_rhs_dict_binder:
			self.apply_constraints_from_row_rhs_format(row_rhs_constraint_group)

	def apply_constraints_from_row_rhs_format(self, c):
		"""
		Inputs:
			c (dict): a dictionary with these keys and values:
				'rows': a list of lists, with each list a vector of all coefficients for a constraint.
				'equality': a list of integers, and each is either -1,0,or 1.
							-1 indicates the constraint is "greater than"
							0 indicates the constraint is at equality
							1 indicates the constraint is "less than."
				'rhs': a list of doubles giving the right hand side of each constraint.
		"""
		for row, sign, rhs in zip(c['rows'], c['equality'], c['rhs']):
			if sign != 0:
				self.domain_coeff_matrix.append([sign*item for item in row])
				self.constraint_ubs.append(sign*rhs)
			else:
				self.domain_coeff_matrix_equality.append(row)
				self.constraint_equality_values.append(rhs)

	def apply_constraints_from_highs_dict_format(self, c, equality=False):
		"""
		Inputs:
			c (dict): a dictionary with these keys:
			'lower', 'upper', 'num_nz', 'start', 'index', 'value'
			see HiGHs documentation for their meanings.
			equality (boolean): indicates if these are equality constraints.
		
		"""
		start_plus_inf = [int(item) for item in c['start']] + [np.inf]
		a_base_array = []
		a_full_array = []
		b_full_vector = []
		for i_nz in range(len(c['start'])):
			this_row = [0 for i in range(self.graph.m)]
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

	def apply_comb_inequality_constraint(self, handle_edges, tooth_edges, rhs):
		coeff_row = [0 for i in range(self.graph.m)]
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
		self.apply_constraints_from_highs_dict_format(self.graph.deg_2_constraints_dict, equality)	
