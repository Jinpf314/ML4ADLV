import functions
import db
import sys
import random
import sage.all
import argparse

parser = argparse.ArgumentParser(description='Generation program for dataset 2. You can run this program several times in parallel to populate the database faster. Make sure to start the different processes with unique random seeds. The Cartan Type used in the original paper is A4.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
db.argparseAddCommonArguments(parser, tablename_default='dataset2')
parser.add_argument('--target', type=int, default=5000, help='How many elements should be generated')
args = parser.parse_args()
db.initializeFromArgparse(args, createTableData=', is_chosen BOOL NOT NULL DEFAULT FALSE')

target = args.target
rank = len(functions.simple_roots)

b_basic = functions.SigmaConjugacyClass(functions.coweight_space.zero(),functions.coweight_lattice.zero())
Wlist = list(functions.W)
while True:
	cnt = db.fetch('SELECT COUNT(*) AS cnt FROM [TABLENAME] WHERE newton_pairings = %s AND dimension IS NOT NULL AND is_chosen = TRUE', (db.generateNewtonPairings(b_basic),))[0]['cnt']
	if cnt >= target:
		break
	mu_indices = [0]*(rank+1)
	while sum(mu_indices) != 0 or len(set(mu_indices)) != rank+1:
		mu_indices = list(sorted(random.randint(-6,6) for _ in range(rank+1)))
	x = random.choice(Wlist)
	y = random.choice(Wlist)
	mu = sum(functions.fundamental_coweights[i+1] * (mu_indices[i+1] - mu_indices[i]) for i in range(rank))
	w = functions.We_W0P.from_classical_weyl(x) * functions.We_W0P.from_translation(mu) * functions.We_W0P.from_classical_weyl(y)
	print(w)
	dim,comp = db.ADLVDimensions(w, b_basic)
	print("dim: %d, components: %d" % (dim,comp))
	if comp > 0:
		db.exec('UPDATE [TABLENAME] SET is_chosen = TRUE WHERE w_encoded = %s AND newton_pairings = %s', (db.encodeW(w), db.generateNewtonPairings(b_basic)))

