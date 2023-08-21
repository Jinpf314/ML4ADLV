import sage.all
import functions
import db
import csv
import sys
import argparse

parser = argparse.ArgumentParser(description='Exports data as CSV. Use the same Cartan Type and tablename as in the generation programs', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
db.argparseAddCommonArguments(parser, tablename_default='dataset2')
parser.add_argument('--limit', type=int, default=5000, help='How many elements should be exported')
args = parser.parse_args()
db.initializeFromArgparse(args, createTableData=False)

limit = args.limit
rank = len(functions.simple_roots)

statement = 'SELECT w_encoded, dimension, components FROM [TABLENAME] WHERE is_chosen=TRUE'
if limit is not None:
	statement += ' LIMIT %d' % limit

fieldnames=['index', 'l(w)', 'l(x)', 'l(y)', 'l(xy)', 'l(yx)', 'l(y*x)', 'l(y\lefttriangle x)']
for i in range(rank+1):
	fieldnames.append('\mu_%d' % (i+1))
for i in range(rank+1):
	for j in range(i+1,rank+1):
		fieldnames.append('\delta(x_%d%d)'% (i+1, j+1))
for i in range(rank+1):
	for j in range(i+1,rank+1):
		fieldnames.append('\delta(y^-1_%d%d)'% (i+1, j+1))
fieldnames += ['dimension', 'components']


rows = db.fetch(statement)
writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
writer.writeheader()
index = 0
for row in rows:
	w = db.decodeW(row['w_encoded'])
	dim = int(row['dimension'])
	comp = int(row['components'])
	w_classical = functions.W(w.to_classical_weyl())
	w_translation = functions.coweight_lattice(w.to_translation_right())
	mu, y_word = w_translation.to_dominant_chamber(reduced_word=True)
	y = functions.W.from_reduced_word(reversed(y_word))
	x = w_classical * y.inverse()
	assert(w == functions.We_W0P.from_classical_weyl(x) * functions.We_W0P.from_translation(mu) * functions.We_W0P.from_classical_weyl(y))
	#print(row['w_encoded'])
	#print(x.reduced_word(), mu, y.reduced_word())
	mu_ambient = sum(c*cow for (c,cow) in zip(functions.toCorootCoefficients(mu), functions.ambient_space.simple_coroots()))
	index += 1
	output = {'index': index, 'l(w)': w.length(), 'l(x)': x.length(), 'l(y)': y.length(), 'l(xy)': (x*y).length(), 'l(yx)': (y*x).length(), 'dimension': dim, 'components': comp}
	xinv = x.inverse()
	for i in range(rank+1):
		output['\mu_%d' % (i+1)] = mu_ambient[i]
		for j in range(i+1,rank+1):
			root = sum(functions.simple_roots[k+1] for k in range(i,j))
			root_as_coweight = functions.corootToCoweightLattice(root)
			ximg = functions.toCorootCoefficients(functions.fastInverseWeylAction(xinv, root_as_coweight))
			assert(root.is_positive_root())
			output['\delta(x_%d%d)' % (i+1,j+1)] = int(all(c <= 0 for c in ximg))
			yimg = functions.toCorootCoefficients(functions.fastInverseWeylAction(y, root_as_coweight))
			output['\delta(y^-1_%d%d)' % (i+1,j+1)] = int(all(c <= 0 for c in yimg))
	lengths = [(y*u).length() for u in functions.bruhatSmallerElements(x)]
	output['l(y*x)'] = max(lengths)
	output['l(y\lefttriangle x)'] = min(lengths)
	writer.writerow(output)
