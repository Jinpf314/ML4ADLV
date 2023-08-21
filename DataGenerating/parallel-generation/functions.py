from sage.all import *
import collections
import random

identity=lambda x:x
sigma_af = None
sigma_classical = None
sigma_translation = None
sigma_permutation = None
is_quasi_split = None
is_split = None
def setRandomSeed(seed):
	random.seed(seed)
	set_random_seed(seed)
setRandomSeed(1)
def InitializeRootSystem(D, sigma=None):
	global R,W,We,We_W0P,We_FW,theta,root_poset,root_lattice,root_space,coroot_lattice,coroot_space,coweight_lattice,coweight_space
	global Rd,Wd
	global cartan_matrix,cartan_inverse,roots_to_ambient,roots_to_ambient_inverse,ambient_space
	global simple_roots,simple_coroots,simple_reflections,simple_affine_reflections,fundamental_coweights
	global sigma_af,sigma_classical,sigma_translation,sigma_permutation,is_quasi_split,is_split,sigma_classical_orbits,sigma_classical_ambient,sigma_classical_ambient_inv
	R = RootSystem(D.cartan_type())
	cartan_matrix = R.cartan_matrix()
	W = WeylGroup(cartan_matrix)
	We = ExtendedAffineWeylGroup(W)
	We_W0P = We.W0P()
	We_FW = We.FW()
	root_poset = R.root_poset()
	theta = root_poset.maximal_elements()[0].element
	root_lattice = R.root_lattice()
	root_space = R.root_space()
	coroot_lattice = R.coroot_lattice()
	coroot_space = R.coroot_space()
	coweight_lattice = R.coweight_lattice()
	coweight_space = R.coweight_space()
	cartan_inverse = cartan_matrix.inverse()
	Rd = R.dual
	Wd = WeylGroup(Rd)
	simple_roots = root_space.alpha()
	simple_coroots = coroot_space.alpha()
	simple_reflections = W.simple_reflections()
	simple_affine_reflections = We_FW.simple_reflections()
	fundamental_coweights = coweight_lattice.basis()
	if sigma is None or sigma in ['s', 'split']:
		sigma = PermutationGroupElement(())
	elif sigma in ['qs', 'quasi-split']:
		automorphisms = D.automorphism_group()
		assert(automorphisms.order() == 2)
		sigma = next(iter(s for s in automorphisms if not s.is_one()))
	else:
		sigma = PermutationGroupElement(sigma)
	sigma_af = sigma
	special_node = sigma_af(0)
	if special_node == 0:
		sigma_classical = sigma
		sigma_translation = coweight_lattice.zero()
		sigma_permutation = W.one()
		is_quasi_split = True
		is_split = all(i == sigma(i) for i in R.index_set())
	else:
		fund = We.fundamental_group()
		fund_action = fund.action(-special_node)
		sigma_classical = Permutation([fund_action(sigma(i)) for i in range(1,len(simple_roots)+1)])
		el = We_W0P(fund[special_node])
		sigma_translation = el.to_translation_right()
		sigma_permutation = el.to_classical_weyl()
		is_quasi_split = False
		is_split = False
	sigma_classical_orbits = set(tuple(sorted((sigma_classical.orbit(i)))) for i in R.index_set())
	ambient_space = R.ambient_space()
	vectors = [v for v in ambient_space.simple_roots().values()]
	'''
	if len(vectors) < ambient_space.dimension():
		vectors.append(sum(ambient_space.basis()).dense_coefficient_list())
	for b in ambient_space.basis().values():
		if ambient_space.submodule(vectors + [b]).rank() > len(vectors):
			vectors.append(b)
	'''
	vectors = [v.dense_coefficient_list() for v in vectors]
	#print(vectors)
	vectors += matrix(vectors).transpose().kernel().basis()
	#print(vectors)
	roots_to_ambient = matrix(vectors).transpose()
	roots_to_ambient_inverse = roots_to_ambient.inverse()
	sigma_classical_ambient = matrix([vectors[sigma_classical(i+1)-1] for i in range(len(vectors))]).transpose() * roots_to_ambient_inverse
	'''
	print(matrix([vectors[sigma_classical(i+1)-1] for i in range(len(vectors))]).transpose())
	print("*")
	print(roots_to_ambient_inverse)
	print("=")
	print(sigma_classical_ambient)
	'''
	sigma_classical_ambient_inv = sigma_classical_ambient.inverse()
def Frobenius_W(w):
	if is_split:
		return w
	'''
	result = W.from_reduced_word(sigma_classical(i) for i in w.reduced_word())
	print("w:")
	print(w.matrix())
	print("sigma(w):")
	print(result.matrix())
	print("sigma matrix:")
	print(sigma_classical_ambient)
	assert(result.matrix() == sigma_classical_ambient * w.matrix() * sigma_classical_ambient_inv)
	return result
	'''
	result = WeylGroupElement(W, sigma_classical_ambient * w.matrix() * sigma_classical_ambient_inv)
	'''
	print(w)
	print(" --- >")
	print(result)
	assert(result.matrix() == W.from_reduced_word(sigma_classical(i) for i in w.reduced_word()).matrix())
	'''
	return result
def Frobenius_Coweight(coweight):
	coeffs = coweight.dense_coefficient_list()
	fc = coweight.parent().basis()
	return sum(fc[sigma_classical(i)] * coeffs[i-1] for i in R.index_set())
def Frobenius_Average(coweight):
	result = coweight_space.zero()
	coeffs = coweight_space(coweight).dense_coefficient_list()
	fc = coweight_space.basis()
	for orbit in sigma_classical_orbits:
		coeff = sum(coeffs[i-1] for i in orbit) / len(orbit)
		result += sum(coeff*fc[i] for i in orbit)
	return result
def Frobenius_Parabolic_Closure(index_set):
	result = set()
	for i in index_set:
		result = result.union(sigma_classical.orbit(i))
	return tuple(sorted(result))
def Frobenius_Support(w):
	w = We_FW(w)
	fund = w.to_fundamental_group()
	w_af = w.to_affine_weyl_right()
	result = set(w_af.support())
	fundmap = fund.parent().action(fund.value())
	while True:
		new_elements = set(sigma_af(fundmap(i)) for i in result)
		if new_elements.issubset(result):
			break
		result = result.union(new_elements)
	return tuple(sorted(result))
def fastInverseWeylAction(groupElement,coweight):
	coeffs = coweight.dense_coefficient_list()
	while len(coeffs) < ambient_space.dimension():
		coeffs = list(coeffs)+[random.randint(-10,10)]
	coeffs = vector(coeffs)
	'''
	#print(coweight)
	#print("input coefficients: %s" % repr(coeffs))
	result = coweight.weyl_action(Wd(groupElement.inverse()))
	out_coeffs = result.dense_coefficient_list()
	if len(out_coeffs) < ambient_space.dimension():
		out_coeffs = list(out_coeffs)+[0]
	out_coeffs = vector(out_coeffs)
	print(result)
	print("output coefficients: %s" % repr(out_coeffs))
	root_lattice_el = sum(coeffs[i-1] * simple_roots[i] for i in R.index_set())
	print("w(%s) = %s" % (root_lattice_el, root_lattice_el.weyl_action(w.inverse())))
	print(roots_to_ambient_inverse.transpose() * groupElement.matrix() * roots_to_ambient.transpose() * vector(root_lattice_el.dense_coefficient_list()+[0]))
	print("ambient matrix:")
	print(groupElement.matrix())
	print("root matrix:")
	print(roots_to_ambient)
	print(roots_to_ambient * groupElement.matrix().transpose() * roots_to_ambient_inverse * coeffs)
	assert(out_coeffs == (roots_to_ambient * groupElement.matrix().transpose() * roots_to_ambient_inverse) * coeffs)
	return result
	'''
	#print(groupElement.matrix())
	#print(roots_to_ambient_inverse.transpose() * groupElement.matrix() * roots_to_ambient.transpose())
	coeffs = roots_to_ambient.transpose() * (groupElement.matrix().transpose() * (roots_to_ambient_inverse.transpose() * coeffs))
	#coeffs = (roots_to_ambient_inverse.transpose() * groupElement.matrix() * roots_to_ambient.transpose()) * coeffs
	if len(simple_roots) < ambient_space.dimension():
		coeffs = vector(list(coeffs)[:len(simple_roots)-ambient_space.dimension()])
	result = coweight.parent().from_vector(coeffs)
	#assert(result.dense_coefficient_list() == coweight.weyl_action(Wd(groupElement.inverse())).dense_coefficient_list())
	return result


def toCorootCoefficients(coweight):
	coweight = coweight_space(coweight)
	coeffs = coweight.dense_coefficient_list()
	return tuple(cartan_inverse.transpose() * sage.all.vector(coeffs))
def corootToCoweightLattice(coroot):
	return coweight_lattice.from_vector(cartan_matrix.transpose() * sage.all.vector(coroot.dense_coefficient_list()))
def corootToCoweightSpace(coroot):
	return coweight_space.from_vector(cartan_matrix.transpose() * sage.all.vector(coroot.dense_coefficient_list()))
def coweightLeq(coweight1,coweight2,checkKottwitzPoints=False):
	coeffs = toCorootCoefficients(coweight2-coweight1)
	if checkKottwitzPoints:
		return all(c >= 0 and c.is_integer() for c in coeffs)
	return all(c >= 0 for c in coeffs)
def parabolicAverage(index_set, coweight):
	if len(index_set) == 0:
		return coweight
	coweight = coweight_space(coweight)
	coeffs = coweight.dense_coefficient_list()
	cartan_submatrix = sage.all.matrix([[cartan_matrix[i-1][j-1] for j in index_set] for i in index_set])
	subtract = tuple(sage.all.vector(coeffs[i-1] for i in index_set) / cartan_submatrix)
	return coweight - sum(subtr * simple_coroots[i] for (i,subtr) in zip(index_set, subtract))
def convexHull(coweight):
	index_set = set()
	while True:
		coeffs = coweight.dense_coefficient_list()
		decr = [i+1 for (i,c) in enumerate(coeffs) if c<0]
		if len(decr) == 0:
			return coweight
		index_set = index_set.union(decr)
		coweight = parabolicAverage(index_set, coweight)
def randomExtendedWeylGroupElement(bound=3):
	permutation = W.random_element()
	translation = sum(random.randint(-bound,bound) * f for f in fundamental_coweights)
	return We_W0P.from_classical_weyl(permutation) * We_W0P.from_translation(translation)
def lengthFunctional(w):
	result = {}
	w = We_W0P(w)
	translation = w.to_translation_right()
	permutation = W(w.to_classical_weyl())
	for alpha in root_poset:
		alpha = alpha.element
		result[alpha] = translation.scalar(alpha) +int(not permutation.action(alpha.to_ambient()).is_positive_root())
		result[-alpha] = -result[alpha]
	return result
def lengthPositive(w,lfunct=None):
	w = We_W0P(w)
	translation = w.to_translation_right()
	permutation = W(w.to_classical_weyl())
	#if lfunct is None:
	#	lfunct = lengthFunctional(w)
	translation_dom,rw = translation.to_dominant_chamber(reduced_word=True)
	translation_dom = translation_dom.dense_coefficient_list()
	vmin = W.from_reduced_word(rw).coset_representative(tuple(i for i in R.index_set() if translation_dom[i-1] == 0),side='right')
	#assert(all(lfunct[alpha.element.weyl_action(vmin)] >= 0 for alpha in root_poset))
	worklist = [vmin]
	result = set(worklist)
	translation_replacement = ambient_space.zero()
	for (i,c) in enumerate(toCorootCoefficients(translation)):
		translation_replacement += c * ambient_space.alphacheck()[i+1]
	for alpha in simple_roots:
		assert(translation_replacement.scalar(alpha.to_ambient()) == translation.scalar(alpha))
	yield vmin
	while worklist:
		v = worklist.pop()
		for (i,alpha) in enumerate(simple_roots):
			valpha = v.action(alpha.to_ambient())
			lfun = translation_replacement.scalar(valpha) +int(valpha.is_positive_root()) - int(permutation.action(valpha).is_positive_root())
			assert(lfun >= 0)
			#assert(lfun <= w.length())
			if lfun > 0:
				continue
			v_new = v.apply_simple_reflection(i+1,side='right')
			if v_new in result:
				continue
			result.add(v_new)
			yield v_new
			worklist.append(v_new)
			if len(result) % 100 == 0:
				print("found %d length positive elements, worklist has length %d..." % (len(result),len(worklist)))
def bruhatSmallerElements(w,index_set = (),side='right'):
	found=set()
	candidates=[w.coset_representative(index_set, side=side)]
	while candidates:
		cand = candidates.pop()
		if any(i in index_set for i in cand.descents(side=side)):
			continue
		yield cand
		for next_cand in cand.bruhat_lower_covers():
			if next_cand in found:
				continue
			found.add(next_cand)
			candidates.append(next_cand)
_qbgAllDist_cache = None
_qbgAllWeight_cache = None
def qbgAllDistances_cached():
	global _qbgAllDist_cache
	if _qbgAllDist_cache is not None:
		return _qbgAllDist_cache
	qbg = W.quantum_bruhat_graph()
	_qbgAllDist_cache = {W(u1): {W(u2): v for (u2,v) in row.items()} for (u1,row) in qbg.distance_all_pairs().items()}
	return _qbgAllDist_cache
pqbg_coeff_data = None
all_max_projection_cache = {}
def prepareQbgDistances():
	global pqbg_coeff_data
	if pqbg_coeff_data is not None:
		return
	pqbg_coeff_data = {}
	for i in R.index_set():
		index_set = tuple(j for j in R.index_set() if i != j)
		pqbg = W.quantum_bruhat_graph(index_set)
		factor = 0
		coroot = simple_coroots[i]
		for root in root_poset:
			root = root.element
			if i in root.support():
				factor += coroot.scalar(root)
		weight = R.ambient_space().fundamental_weights()[i]
		lengthData = {weight: 0}
		workList = [weight]
		while workList:
			wt = workList.pop()
			for j in R.index_set():
				if not wt.has_descent(j):
					newWt = wt.simple_reflection(j)
					if not newWt in lengthData:
						lengthData[newWt] = lengthData[wt]+1
						workList.insert(0,newWt)
		#print(lengthData)
		pqbg_coeff_data[i] = {}
		for (u1,row) in pqbg.distance_all_pairs().items():
			output = {}
			u1wt = u1.action(weight)
			pqbg_coeff_data[i][u1wt] =output
			u1len = lengthData[u1wt]
			assert(u1.length() == u1len)
			for (u2,dist) in row.items():
				u2wt = u2.action(weight)
				u2len = lengthData[u2wt]
				coeff = (dist-u2len+u1len) / factor
				#print("PQBG%s dist %s i.e. %s to %s i.e. %s: %d coeff: %d" % (index_set, u1.reduced_word(), u1wt, u2.reduced_word(), u2wt,  dist, coeff))
				assert(coeff.is_integer())
				output[u2wt] = coeff
		#pqbg_coeff_data[i] = {W(u1): {W(u2): (d-u2.length() + u1.length()) / factor for (u2,d) in row.items()} for (u1,row) in pqbg.distance_all_pairs().items()}
def qbgWeight(start,end,returnDistance=False):
	global all_max_projection_cache,pqbg_coeff_data
	if pqbg_coeff_data is None:
		prepareQbgDistances()
	'''
	print("start:")
	print(start)
	print("end:")
	print(end)
	'''
	start=W(start)
	end=W(end)
	startProjections=None#all_max_projection_cache.get(start)
	fundamental_weights = R.ambient_space().fundamental_weights()
	if True:#startProjections is None:
		startProjections = {i: start.action(fundamental_weights[i]) for i in R.index_set()}
		#for i in R.index_set():
		#	assert(startProjections[i] == fundamental_coweights[i].weyl_action(Wd(start.coset_representative(tuple(j for j in R.index_set() if i!= j),side='right'))))
		#all_max_projection_cache[start] = startProjections
	endProjections=None#all_max_projection_cache.get(end)
	if True:#endProjections is None:
		endProjections = {i: end.action(fundamental_weights[i]) for i in R.index_set()}
		#endProjections = {i: end.coset_representative(tuple(j for j in R.index_set() if i != j),side='right') for i in R.index_set()}
		#all_max_projection_cache[end] = endProjections
	'''
	print("qbg weights %s -> %s" % (start.reduced_word(), end.reduced_word()))
	print([(startProjections[i], endProjections[i]) for i in R.index_set()])
	'''
	result = coroot_lattice.zero()
	for i in R.index_set():
		result += simple_coroots[i] * pqbg_coeff_data[i][startProjections[i]][endProjections[i]]
	if returnDistance:
		return result,sum(coroot_lattice(result).dense_coefficient_list())*2 - start.length() + end.length()
	return result
def qbgAllWeights_cached():
	global _qbgAllWeight_cache
	if _qbgAllWeight_cache is not None:
		return _qbgAllWeight_cache
	'''
	parabolic_index_sets = {}
	parabolic_distances = {}
	parabolic_factors = {}
	parabolic_projections = {}
	for i in R.index_set():
		index_set = tuple(j for j in R.index_set() if i != j)
		pqbg = W.quantum_bruhat_graph(index_set)
		distances = {W(u1): {W(u2): d for (u2,d) in row.items()} for (u1,row) in pqbg.distance_all_pairs().items()}
		factor = 0
		coroot = simple_coroots[i]
		for root in root_poset:
			root = root.element
			if i in root.support():
				factor += coroot.scalar(root)
				print("%s -> %d" % (root,coroot.scalar(root)))
		print("factor for %d: %d" % (i,factor))
		parabolic_index_sets[i] = index_set
		parabolic_distances[i] = distances
		parabolic_factors[i] = factor
		parabolic_projections[i] = {u: u.coset_representative(index_set,side='right') for u in W}
	result = {}
	for u1 in W:
		row = {}
		result[u1] = row
		for u2 in W:
			val = coroot_lattice.zero()
			for i in R.index_set():
				p1 = parabolic_projections[i][u1]
				p2 = parabolic_projections[i][u2]
				coeff = ( parabolic_distances[i][p1][p2] - p2.length() + p1.length() ) / factor
				assert(coeff.is_integer())
				val += simple_coroots[i] * coeff
			row[u2] = val
	#print(result)
	_qbgAllWeight_cache = result
	'''
	zero = coroot_lattice.zero()
	result = {u1: {u2: zero if u1 == u2 else None for u2 in W} for u1 in W}
	qbg = W.quantum_bruhat_graph()
	qbgArrowsByDomain = {u: list() for u in W}
	for (start,end,root) in qbg.edges():
		is_quantum = int(start.length() > end.length())
		qbgArrowsByDomain[start].append((end,root.associated_coroot() * is_quantum))
	worklist = [(u,u) for u in W]
	while worklist:
		start,end = worklist.pop()
		for (next_end,added_weight) in qbgArrowsByDomain[end]:
			if result[start][next_end] is None:
				result[start][next_end] = result[start][end] + added_weight
				assert(result[start][next_end] == qbgWeight(start,next_end))
				worklist.insert(0,(start,next_end))
	_qbgAllWeight_cache = result
	return result

def enumerateAdmissibleLocus(mu,index_set):
	mu = mu.to_dominant_chamber()
	mu_coeffs = mu.dense_coefficient_list()
	stabilizer = tuple([i for i in R.index_set() if mu_coeffs[i-1] == 0])
	conjugators = bruhatSmallerElements(W.long_element(), stabilizer, side='right')
	candidates = [We_FW(We_W0P.from_translation(mu.weyl_action(Wd(conj)))).coset_representative(index_set,side='left') for conj in conjugators]
	found = set(repr(c) for c in candidates)
	while candidates:
		cand = candidates.pop()
		if any(cand.has_descent(i,side='left') for i in index_set):
			continue
		yield cand
		fund = We_FW.from_fundamental(cand.to_fundamental_group())
		cand_af = cand.to_affine_weyl_right()
		for next_cand_af in cand_af.bruhat_lower_covers():
			next_cand = fund * We_FW.from_affine_weyl(next_cand_af)
			nc_repr = repr(next_cand)
			if nc_repr in found:
				continue
			found.add(nc_repr)
			candidates.append(next_cand)
class SigmaConjugacyClass:
	def __init__(self,newton,kottwitz_repr):
		self.newton_point = coweight_space(newton)
		self.kottwitz_repr = kottwitz_repr
		self.kottwitz_number = min(Frobenius_Parabolic_Closure((We_FW.from_translation(kottwitz_repr).to_fundamental_group().value(),)))
		coeffs = toCorootCoefficients(newton - Frobenius_Average(kottwitz_repr))
		self.J1_repr = set()
		self.J1 = set()
		for orbit in sigma_classical_orbits:
			c = sum(coeffs[i-1] for i in orbit)
			if not c.is_integer():
				self.J1_repr.add(next(iter(orbit)))
				self.J1 = self.J1.union(orbit)
		self.defect = len(self.J1_repr)
		newton_coeffs = coweight_space(newton).dense_coefficient_list()
		self.J2 = set(i for i in R.index_set() if newton_coeffs[i-1] == 0)
		self.is_basic = (len(self.J2) == len(simple_roots))
		self.newton2rho = 2*sum(toCorootCoefficients(newton))
	def __hash__(self):
		return hash((tuple(self.newton_point.dense_coefficient_list()),self.kottwitz_number))
	def __eq__(self,other):
		return isinstance(other,SigmaConjugacyClass) and self.newton_point == other.newton_point and self.kottwitz_number == other.kottwitz_number
	def __le__(self,other):
		return self.kottwitz_number == other.kottwitz_number and coweightLeq(self.newton_point, other.newton_point)
	def __lt__(self,other):
		return self.kottwitz_number == other.kottwitz_number and coweightLeq(self.newton_point, other.newton_point) and self.newton_point != other.newton_point
	def __repr__(self):
		return "(nu=%s, kappa=%d)" % (self.newton_point, self.kottwitz_number)
	def lambdaRepresentative(self):
		coeffs = toCorootCoefficients(self.newton_point - Frobenius_Average(self.kottwitz_repr))
		coeffs_floor = [floor(c) for c in coeffs]
		return self.kottwitz_repr + sum(c * alpha for (c,alpha) in zip(coeffs_floor, coweight_lattice.alpha()))
def W0P_hash(w):
	w = We_W0P(w)
	translation = coweight_space(w.to_translation_right())
	permutation = w.to_classical_weyl()
	return (permutation.matrix(), tuple(translation.dense_coefficient_list()))
classPolynomialQ = var('q')
classPolynomialRing = ZZ[classPolynomialQ]
def ADLVDimensions(w,verbose=None,outstream=None,computeClassPolynomials=False):
	if verbose is not None:
		outstream = outstream or sys.stdout
		outstream.write("%s%s" % (verbose, We_FW(w)))
		candidate_path_text = {W0P_hash(w): "\n"}
	worklist = [w]
	found = set([W0P_hash(w)])
	while worklist:
		cand = worklist.pop()
		for i in [0]+list(R.index_set()):
			if not cand.has_descent(i,side='left'):
				continue
			cand2 = cand.apply_simple_reflection(i,side='left')
			iright = sigma_af(i)
			if cand2.has_descent(iright,side='right'):
				verbose1=None
				verbose2=None
				if verbose is not None:
					text = candidate_path_text[W0P_hash(cand)]
					assert(text.startswith("\n") and text.endswith("\n"))
					outstream.write(candidate_path_text[W0P_hash(cand)])
					verbose1=verbose + " O "
					verbose2=verbose + " I "
				dim1 = ADLVDimensions(cand2,verbose1,outstream=outstream,computeClassPolynomials=computeClassPolynomials)
				dim2 = ADLVDimensions(cand2.apply_simple_reflection(iright,side='right'),verbose2,outstream=outstream,computeClassPolynomials=computeClassPolynomials)
				if computeClassPolynomials:
					zeroPoly = classPolynomialRing(0)
					for b in set(dim1).union(dim2):
						dim2[b] = classPolynomialRing(dim1.get(b,zeroPoly)*(classPolynomialQ-1) + dim2.get(b,zeroPoly)*classPolynomialQ)
				else:
					zerocounter = collections.Counter({})
					for b in set(dim1).union(dim2):
						dim2[b] = collections.Counter({d+1:v for (d,v) in (dim1.get(b,zerocounter)+dim2.get(b,zerocounter)).items()})
				return dim2
			cand2 = cand2.apply_simple_reflection(iright,side='right')
			cand2_hash = W0P_hash(cand2)
			if cand2_hash in found:
				continue
			found.add(cand2_hash)
			worklist.insert(0,cand2)
			if verbose is not None:
				candidate_path_text[cand2_hash] = candidate_path_text[W0P_hash(cand)]+"%s ~ %s via %d\n" % (verbose, repr(cand2),i)
	if verbose is not None:
		outstream.write("\n")
	w = We_W0P(w)
	permutation = W(w.to_classical_weyl())
	translation = coweight_lattice(w.to_translation_right())
	factor = 1
	prod = permutation
	perm_inv = permutation.inverse()
	summation = coweight_space(translation)
	assert(Frobenius_Coweight(fastInverseWeylAction(perm_inv,summation)) == fastInverseWeylAction(Frobenius_W(perm_inv),Frobenius_Coweight(summation)))
	while not prod.is_one() or not (sigma_af ** factor).is_one():
		factor += 1
		summation = translation + Frobenius_Coweight(fastInverseWeylAction(perm_inv,summation))
		prod = permutation*Frobenius_W(prod)
	nu = coweight_space(summation)/factor
	nu = nu.to_dominant_chamber()
	assert(nu == Frobenius_Coweight(nu))
	b = SigmaConjugacyClass(nu,translation)
	if verbose is not None:
		outstream.write("%sNP %s %s defect %d dimension %d\n" % (verbose, b.newton_point, toCorootCoefficients(b.newton_point), b.defect, w.length() - b.newton2rho))
	return {b: classPolynomialRing(classPolynomialQ ** (w.length()-b.newton2rho))} if computeClassPolynomials else {b: collections.Counter({w.length() - b.newton2rho:1})}
_basedDBGPaths_cache = {}
def listBasedDBGPaths(start, base):
	start = W(start)
	start_inv = start.inverse()
	base = W(base)
	key = (start.matrix(),base.matrix())
	if key in _basedDBGPaths_cache:
		return _basedDBGPaths_cache[key]
	result = {w.matrix(): [] for w in W}
	ascents = [i for i in R.index_set() if not base.has_descent(i,side='right')]
	if not ascents:
		result[start.matrix()].append([])
		return result
	index = next(iter(ascents))
	next_base = base.apply_simple_reflection(index, side='right')
	#assert(next_base.length() == base.length()+1)
	root = root_lattice.simple_roots()[index]
	root = root.weyl_action(base)
	next_start = start * W.from_reduced_word(root.associated_reflection())
	result1 = listBasedDBGPaths(start, next_base)
	result2 = listBasedDBGPaths(next_start, next_base)
	for (end, paths) in result1.items():
		result[end] += paths
	coroot = coweight_lattice(root.associated_coroot())
	for (end, paths) in result2.items():
		#edge = (start, next_start, root, not root.weyl_action(start).is_positive_root())
		is_positive = any(c>0 for c in toCorootCoefficients(fastInverseWeylAction(start_inv,coroot)))
		edge = (start, next_start, root, not is_positive)
		paths = [ [ edge ] + path for path in paths ]
		result[end] += paths
	_basedDBGPaths_cache[key] = result
	return result
def listDBGPaths(start,end,target=None):
	if target is None:
		target=end
	base = end.inverse()*target
	return [[(end,start,root,not quantum) for (start,end,root,quantum) in reversed(p)] for p in listBasedDBGPaths(end,base)[start.matrix()]]
_reduced_word_cache = {}
def reduced_word_cached(w):
	key = w.matrix()
	if key not in _reduced_word_cache:
		_reduced_word_cache[key] = w.reduced_word()
	return _reduced_word_cache[key]
def enumeratePathWeightsWithBound(path, bound,frobeniusAverage=False):
	result = {}
	minimalWeight = sum(is_quantum * coroot_lattice(root.associated_coroot()) for (start,end,root,is_quantum) in path)
	result[minimalWeight] = 1
	for edge in path:
		root = edge[2]
		coroot = coroot_lattice(root.associated_coroot())
		for (wt, mult) in list(result.items()):
			k=0
			while True:
				k += 1
				newWt = wt + k*coroot
				if not coweightLeq(Frobenius_Average(newWt) if frobeniusAverage else newWt, bound):
					break
				result[newWt] = result.get(newWt,0) + mult
	return result
def IHProduct(w1, w2):
	w1_word = w1.to_affine_weyl_right().reduced_word()
	w1_fun = We_W0P.from_fundamental(w1.to_fundamental_group())
	result = collections.Counter({(w2,0):1})
	for letter in reversed(w1_word):
		new_result = collections.Counter()
		for ((w,e),n) in result.items():
			new_result[(w.apply_simple_reflection(letter,side='left'),e)] += n
			if w.has_descent(letter,side='left'):
				new_result[(w,e+1)]+=n
		result=new_result
	return collections.Counter({(w1_fun*w,e):n for ((w,e),n) in result.items()})
def argparseAddCommonArguments(parser):
	parser.add_argument('CartanType', default='A2', help='Cartan type of the finite root system, such as A2, B3, F4 etc.')
	parser.add_argument('--frobenius', default='split',help="Frobenius action. Should be 'split', 'quasi-split' or a permutation in cycle notation. The word 'quasi-split' is only available if there is a unique non-trivial automorphism of the finite Dynkin diagram")
	parser.add_argument('--random-seed', default=1, type=int)
def initializeFromArgparse(args):
	global D
	D = DynkinDiagram(args.CartanType)
	InitializeRootSystem(D, args.frobenius)
	setRandomSeed(args.random_seed)
'''
def listDbgWeightsWithBound(start,end,bound,target=None,frobeniusAverage=False):
	result = {}
	for path in listDBGPaths(start,end,target):
		for (wt,num) in enumeratePathWeightsWithBound(path,bound,frobeniusAverage).items():
			result[wt] = result.get(wt,0)+num
	return result
'''
if __name__ == '__main__':
	import itertools
	D = DynkinDiagram(['D', 5])
	print(D)
	InitializeRootSystem(D,'split')
	'''
	for (start,end) in itertools.product(W,W):
		for target in W:
			print("Calculating wts(%s => %s --> %s)" % (start.reduced_word(), end.reduced_word(), target.reduced_word()))
			paths = listDBGPaths(start,end,target)
			for p in paths:
				print("\tpath of length %d" % len(p))
				for edge in p:
					print("\t\t%s -> %s via root %s (quantum: %s)" % (edge[0].reduced_word(), edge[1].reduced_word(), edge[2], edge[3]))
				bound = sum(40*c for c in coweight_lattice.basis())
				for (wt,num) in enumeratePathWeightsWithBound(p,bound).items():
					print("\t%s %s: %d" % (wt,toCorootCoefficients(wt),num))
			el1 = start.inverse()*target
			el2 = end.inverse()*target
			#print(repr(el1.reduced_word())+","+repr(el2.reduced_word()))
			assert( (len(paths)>0) == (el2.bruhat_le(el1) ))
			if paths:
				assert(max(len(p) for p in paths) == el1.length()-el2.length())
	'''
	coweight = coweight_lattice.an_element()
	print("%s -> %s. avg: %s" % (coweight, Frobenius_Coweight(coweight), Frobenius_Average(coweight)))
	w = W.an_element()**2
	print("%s -> %s" % (w.reduced_word(), Frobenius_W(w).reduced_word()))
	w_af = We_W0P.from_classical_weyl(w) * We_W0P.from_translation(coweight)
	print(w_af)
	print(w_af.length())
	w0 = We_W0P.from_classical_weyl(W.long_element())
	rho = sum(fundamental_coweights)
	w2 = We_W0P.from_translation(rho)
	w1 = w2*w0
	print("w1: %s length %d" % (w1, w1.length()))
	print("w2: %s length %d" % (w2, w2.length()))
	print("w1w2: %s length %d" % (w1*w2, (w1*w2).length()))
	print("w2w1: %s length %d" % (w2*w1, (w2*w1).length()))
	prod = IHProduct(w1,w2**10)
	print(prod)
	minValue = 10
	for w,e in prod:
		minValue = min(minValue, min(w.to_translation_right().dense_coefficient_list()))
	print("Minimal value: %d" % minValue)
	'''
	total=0
	for (b,d) in sorted(ADLVDimensions(w_af,"",computeClassPolynomials=True).items()):
		print("%s, %s: %s = %s deg %d" % (b.newton_point, toCorootCoefficients(b.newton_point), d,d.expand(), d.degree(classPolynomialQ)))
		total+=d * classPolynomialQ ** b.newton2rho
	print("total: %s" % repr(total.expand()))
	'''
