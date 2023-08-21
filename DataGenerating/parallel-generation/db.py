import functions
import importlib
import hashlib
import sage.all
import sys
import argparse

dbHandler = None
dbCursor = None
dbTablename = None
def openDatabase(modulename = 'MySQLdb', host='localhost', port=3306, user='root', password='', db='adlv',tablename='adlv'):
	global dbHandler, dbCursor, dbTablename
	db_module = importlib.import_module(modulename)
	dbHandler = db_module.connect(host=host, port=port, user=user, password=password, db=db)
	dbHandler.autocommit(True)
	dbTablename = tablename
	dbCursor = dbHandler.cursor()

def fetch(structure, data=()):
	#print("FETCH: %s %s" % (structure,data))
	dbCursor.execute(structure.replace('[TABLENAME]', dbTablename), data)
	rows = dbCursor.fetchall()
	return [ {desc[0]:value for (desc,value) in zip(dbCursor.description, row)} for row in rows ]
def exec(structure, data=(), expectAffected=None):
	#print("EXEC: %s %s" % (structure,data))
	affected=dbCursor.execute(structure.replace('[TABLENAME]', dbTablename), data)
	if expectAffected is not None and expectAffected != affected:
		raise Exception("Affected %d rows instead of %d" % (affected, expectAffected))
	return dbCursor.lastrowid
def createTable(extra = ''):
	exec('''
CREATE TABLE IF NOT EXISTS [TABLENAME]
(w_encoded VARCHAR(255) NOT NULL,
newton_pairings VARCHAR(255) NOT NULL,
is_complete BOOL NOT NULL,
is_empty BOOL NULL,
dimension INT NULL,
components INT NULL,
polynomial VARCHAR(255) NULL,
PRIMARY KEY (w_encoded, newton_pairings)%s)''' % extra)

def encodeW(w):
	key = str(functions.R.cartan_type())
	key += ";" + repr(functions.sigma_af)
	key += ";" + repr(w.to_classical_weyl().reduced_word())
	key += ";" + repr(functions.coweight_lattice(w.to_translation_right()).to_vector())
	#return hashlib.sha1(key).hexdigest()
	return key
def decodeW(key):
	parts = key.split(";")
	assert(len(parts) == 4)
	assert(parts[0] == str(functions.R.cartan_type()))
	assert(parts[1] == repr(functions.sigma_af))
	assert(parts[2].startswith("[") and parts[2].endswith("]"))
	assert(parts[3].startswith("(") and parts[3].endswith(")"))
	classical = functions.W.from_reduced_word(int(s) for s in parts[2][1:-1].split(","))
	mu_coeffs = [int(c) for c in parts[3][1:-1].split(",")]
	assert(len(mu_coeffs) == len(functions.fundamental_coweights))
	mu = sum(functions.coweight_lattice(c * cow) for (c,cow) in zip(mu_coeffs, functions.fundamental_coweights))
	return functions.We_W0P.from_classical_weyl(classical) * functions.We_W0P.from_translation(mu)
def generateNewtonPairings(b):
	newton = b.newton_point
	pairings = tuple(newton.scalar(alpha) for alpha in functions.simple_roots)
	return ','.join(str(p) for p in pairings)
def fetchInfo(w, b):
	wenc = encodeW(w)
	newton_pairings = generateNewtonPairings(b)
	statement = '''SELECT is_complete, is_empty, dimension, components, polynomial FROM [TABLENAME]'''
	statement += ''' WHERE w_encoded = %s AND newton_pairings = %s'''
	rows = fetch(statement,(wenc,newton_pairings))
	if len(rows) == 0:
		statement = '''SELECT COUNT(*) AS cnt FROM [TABLENAME] WHERE w_encoded = %s AND is_complete = TRUE LIMIT 1'''
		rows = fetch(statement,(wenc,))
		row = rows[0]
		return {'empty': None if row['cnt'] == 0 else True, 'dimension': None, 'components': None, 'polynomial': None}
	assert(len(rows) == 1)
	result = rows[0]
	polynomial = None
	if result['polynomial'] is not None:
		coefficients = [int(c) for c in result['polynomial'].split(',')]
		polynomial = functions.classPolynomialRing(sum(c * (functions.classPolynomialQ ** i) for (i,c) in enumerate(coefficients)))
	return {'empty': bool(result['is_empty']) if result['is_empty'] is not None else None,
		'dimension': int(result['dimension']) if result['dimension'] is not None else None,
		'components': int(result['components']) if result['components'] is not None else None,
		'polynomial': polynomial}
def fetchCompleteInfo(w):
	wenc = encodeW(w)
	statement = 'SELECT newton_pairings, dimension, components, polynomial'
	statement += ' FROM [TABLENAMe] WHERE w_encoded = %s AND is_complete = TRUE AND is_empty = FALSE'
	rows = fetch(statement,(wenc,))
	if len(rows) == 0:
		return None
	result = {}
	mu = w.to_translation_right()
	for row in rows:
		data = {desc[0]: value for (desc,value) in zip(dbCursor.description, row)}
		newton_vector_coefficients = [sage.all.QQ(s) for s in data['newton_pairings'].split(',')]
		newton_vector = sum(functions.coweight_space(cow) * coef for (cow,coef) in zip(functions.fundamental_coweights, newton_vector_coefficients))
		b = functions.SigmaConjugacyClass(newton_vector,mu)
		res = {}
		if data['polynomial'] is None:
			res['polynomial'] = None
		else:
			coefficients = [int(c) for c in data['polynomial'].split(',')]
			polynomial = functions.classPolynomialRing(sum(c * (functions.classPolynomialQ ** i) for (i,c) in enumerate(coefficients)))
			res['polynomial'] = polynomial
		res['dimension'] = data['dimension']
		res['components'] = data['components']
		result[b] = res
	return result
def insertSingleInfo(w, b, empty=None, dimension=None, components=None, polynomial=None, complete=None):
	polynomial_str = None if polynomial is None else ','.join(str(c) for c in polynomial.coefficients(sparse=False))
	wenc = encodeW(w)
	newton_pairings = generateNewtonPairings(b)
	statement = '''SELECT COUNT(*) AS cnt FROM [TABLENAME] WHERE w_encoded = %s AND newton_pairings = %s'''
	rows= fetch(statement,(wenc,newton_pairings))
	count = rows[0]['cnt']
	statements = []
	data = []
	if empty is not None:
		data.append(int(empty))
		statements.append('is_empty=%s')
	if dimension is not None:
		data.append(dimension)
		statements.append('dimension=%s')
	if components is not None:
		data.append(components)
		statements.append('components=%s')
	if polynomial is not None:
		data.append(polynomial_str)
		statements.append('polynomial=%s')
	if complete is not None:
		data.append(int(complete))
		statements.append('is_complete=%s')
	if count == 0:
		if complete is None:
			data.append(0)
			statements.append('is_complete=%s')
		statement1 = '''INSERT IGNORE INTO [TABLENAME] (w_encoded, newton_pairings'''
		statement2 = ') VALUES (%s, %s'
		for st in statements:
			st1,st2 = st.split('=')
			statement1 += ', '+st1
			statement2 += ', '+st2
		data = [wenc, newton_pairings] + data
		statement = statement1 + statement2 + ')'
		exec(statement, data)
	else:
		data += [wenc, newton_pairings]
		statement = '''UPDATE [TABLENAME] SET '''
		statement += ', '.join(statements)
		statement += ' WHERE w_encoded = %s AND newton_pairings = %s LIMIT 1'
		exec(statement, tuple(data))
def ADLVDimensions(w, b=None, classPolynomial=False,verbose=None,outstream=None):
	if verbose is not None:
		outstream = outstream or sys.stdout
		outstream.write("%s%s" % (verbose, functions.We_FW(w)))
		candidate_path_text = {functions.W0P_hash(w): "\n"}
	worklist = [w]
	found = set([functions.W0P_hash(w)])
	toUpdate = []
	result = None
	while worklist:
		cand = worklist.pop()
		if b is None:
			info = fetchCompleteInfo(cand)
			if info is not None:
				result = {}
				for (cl,data) in info.items():
					if classPolynomial:
						if data['polynomial'] is None:
							result = None
							break
						result[cl] = data['polynomial']
					else:
						if data['dimension'] is None or data['components'] is None:
							result = None
							break
						result[cl] = (data['dimension'],data['components'])
				if result is not None:
					break
		else:
			data = fetchInfo(cand,b)
			if classPolynomial and data['polynomial'] is not None:
				result = data['polynomial']
				break
			elif not classPolynomial and data['dimension'] is not None and data['components'] is not None:
				result = (data['dimension'],data['components'])
				break
			elif data['empty'] is True:
				result = functions.classPolynomialRing(0) if classPolynomial else (-1,0)
				break
		toUpdate.append(cand)
		for i in [0]+list(functions.R.index_set()):
			if not cand.has_descent(i,side='left'):
				continue
			cand2 = cand.apply_simple_reflection(i,side='left')
			iright = functions.sigma_af(i)
			if cand2.has_descent(iright,side='right'):
				verbose1=None
				verbose2=None
				if verbose is not None:
					text = candidate_path_text[functions.W0P_hash(cand)]
					assert(text.startswith("\n") and text.endswith("\n"))
					outstream.write(candidate_path_text[functions.W0P_hash(cand)])
					verbose1=verbose + " O "
					verbose2=verbose + " I "
				dim1 = ADLVDimensions(cand2,b,verbose=verbose1,outstream=outstream,classPolynomial=classPolynomial)
				dim2 = ADLVDimensions(cand2.apply_simple_reflection(iright,side='right'),b,verbose=verbose2,outstream=outstream,classPolynomial=classPolynomial)
				if classPolynomial and b is None:
					zeroPoly = functions.classPolynomialRing(0)
					result = {}
					q = functions.classPolynomialQ
					for cl in set(dim1).union(dim2):
						result[cl] = functions.classPolynomialRing(dim1.get(cl,zeroPoly)*(q-1) + dim2.get(cl,zeroPoly)*q)
				elif classPolynomial and b is not None:
					q = functions.classPolynomialQ
					result = dim1 * (q-1) + dim2 * q
				elif b is None:
					result = {}
					for cl in set(dim1).union(dim2):
						d1 = dim1.get(cl, (-2,0))
						d2 = dim2.get(cl, (-2,0))
						if d1[0] < d2[0]:
							result[cl] = (d2[0]+1,d2[1])
						elif d1[0] > d2[0]:
							result[cl] = (d1[0]+1,d1[1])
						else:
							result[cl] = (d1[0]+1, d1[1]+d2[1])
				else:
					if dim1[0] < dim2[0]:
						result = (dim2[0]+1,dim2[1])
					elif dim1[0] > dim2[0]:
						result = (dim1[0]+1,dim1[1])
					elif dim1 == (-1,0):
						result = (-1,0)
					else:
						result = (dim1[0]+1,dim1[1]+dim2[1])
				break
			if result is not None:
				break
			cand2 = cand2.apply_simple_reflection(iright,side='right')
			cand2_hash = functions.W0P_hash(cand2)
			if cand2_hash in found:
				continue
			found.add(cand2_hash)
			worklist.insert(0,cand2)
			if verbose is not None:
				candidate_path_text[cand2_hash] = candidate_path_text[functions.W0P_hash(cand)]+"%s ~ %s via %d\n" % (verbose, repr(cand2),i)
	if verbose and result is not None:
		outstream.write("%slookup successful!\n" % verbose)
	if result is None:
		if verbose is not None:
			outstream.write("\n")
		w = functions.We_W0P(w)
		permutation = functions.W(w.to_classical_weyl())
		translation = functions.coweight_lattice(w.to_translation_right())
		factor = 1
		prod = permutation
		perm_inv = permutation.inverse()
		summation = functions.coweight_space(translation)
		#assert(Frobenius_Coweight(fastInverseWeylAction(perm_inv,summation)) == fastInverseWeylAction(Frobenius_W(perm_inv),Frobenius_Coweight(summation)))
		while not prod.is_one() or not (functions.sigma_af ** factor).is_one():
			factor += 1
			summation = translation + functions.Frobenius_Coweight(functions.fastInverseWeylAction(perm_inv,summation))
			prod = permutation*functions.Frobenius_W(prod)
		nu = functions.coweight_space(summation)/factor
		nu = nu.to_dominant_chamber()
		assert(nu == functions.Frobenius_Coweight(nu))
		cl = functions.SigmaConjugacyClass(nu,translation)
		if verbose is not None:
			outstream.write("%sNP %s defect %d dimension %d\n" % (verbose, cl.newton_point, cl.defect, w.length() - cl.newton2rho))
		dim = w.length()-cl.newton2rho
		if b == cl:
			if classPolynomial:
				result = functions.classPolynomialRing(functions.classPolynomialQ ** dim)
			else:
				result = (dim,1)
		elif b is not None and b != cl:
			if classPolynomial:
				result = functions.classPolynomialRing(0)
			else:
				result = (-1,0)
		else:
			if classPolynomial:
				result = {cl: functions.classPolynomialRing(functions.classPolynomialQ ** dim)}
			else:
				result = {cl: (dim,1)}
	update = result if b is None else {b: result}
	complete = True if b is None else None
	for cand in toUpdate:
		for (cl,res) in update.items():
			empty=None
			dimension=None
			components=None
			polynomial=None
			if classPolynomial:
				if res == 0:
					empty=True
				else:
					empty=False
					dimension=res.degree()
					components=res.leading_coefficient()
				polynomial=res
			else:
				if res == (-1,0):
					empty=True
				else:
					empty=False
					dimension=res[0]
					components=res[1]
			insertSingleInfo(cand,cl,empty=empty,dimension=dimension,components=components,polynomial=polynomial,complete=complete)
	return result
def argparseAddCommonArguments(parser,tablename_default='adlv'):
	functions.argparseAddCommonArguments(parser)
	parser.add_argument('--db-host', default='localhost', help='Host parameter for the MySQL connection')
	parser.add_argument('--db-port', default=3306, type=int, help='Port parameter for the MySQL connection')
	parser.add_argument('--db-user', default='root', help='Username parameter for the MySQL connection')
	parser.add_argument('--db-pass', default='', help='Password parameter for the MySQL connection')
	parser.add_argument('--db-database', default='adlv', help='Name of the MySQL database to be used')
	parser.add_argument('--db-tablename', default=tablename_default, help='Name of the MySQL table to be used (will be created if not exists)')
	parser.add_argument('--db-pythonmodule', default='MySQLdb', help='Name of the python3 module to import for the database connection')
def initializeFromArgparse(args,createTableData=None):
	functions.initializeFromArgparse(args)
	openDatabase(host=args.db_host, port=args.db_port, user=args.db_user, password=args.db_pass, db=args.db_database, tablename=args.db_tablename)
	if createTableData is None or createTableData is True:
		createTable()
	elif createTableData is not False:
		createTable(createTableData)
if __name__ == '__main__':
	D = sage.all.DynkinDiagram('A2')
	functions.InitializeRootSystem(D)
	openDatabase()
	el = functions.W.from_reduced_word([1,2])
	w1 = functions.We_W0P.from_classical_weyl(el) * functions.We_W0P.from_translation(functions.coweight_lattice.alpha()[1])
	w2 = functions.We_W0P.from_classical_weyl(el) * functions.We_W0P.from_translation(functions.coweight_lattice.alpha()[2])
	dims1 = ADLVDimensions(w1, verbose='',classPolynomial=True)
	print(dims1)
	sys.exit(0)
	
	#dims2 = functions.ADLVDimensions(w2, '')
	for (b,cnt) in dims1.items():
		dim = max(cnt)
		irred = cnt[dim]
		#insertSingleInfo(w1, b, empty=False, complete=True, dimension=dim, components=irred)
		print(fetchInfo(w1,b))

