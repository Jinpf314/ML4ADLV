run main.py£¬ load data types, variables, and functions. 
Users can use these functions to compute Lie theoretic data and dimension of affine Deligne-Lusztig variety

################################
First input type of group. A4, C5, D4, 2A5, etc.

Variables:

dictionary s, s[0], s[1], s[2], ..., s[n], affine simple reflections

dictionary W, W[i] = set of element in the extended Weyl group \tilde{W} with length i, i = 0,1, ..., 29

dictionary W0, W0[i] = set of element in the finite Weyl group W with length i, i = 0,1, ..., 29

set Phiplus = set of positive roots

set Phi = set of roots

vector two_rho = 2\rho

vector rho = \rho

element w_0: longest element in W

dictionary omegacheck, omega, alpha, alphacheck

dictionary tau, tau[i] = length zero elements in \tilde{W} corresponding to i.

#############################

Functions:

1. For an element w in \tilde{W}, we have


w.present(): print the translation part and finite part of w

w.finite(): return the finite part of w

w.rex(): return a reduced expression of w

w.len(): return the length of w

w.isstraight(): return if w is a straight element

w.Newton(): return the Newton vector of w

is_shrunken(w): return if w is in the shrunken Weyl chamber

partialprint(w): return the expression w = x t^{\mu} y

max_Newton(w): return the generic Newton point of w

is_cordial(w): return if w is a cordial element

BG(w): return the set of Newton points \nu with dim X_w(\mu) is non-empty 

dim_irr_print(w): print dim X_w(\nu) and number of J_b orbits of top-dim irreducible components of X_w(\nu) for all Newton points \nu.

2. Way to generate an element w:

exp(a list): for example, exp([0,1,2]) = s_0s_1s_2

trans(a list): for example, exp([1,2,3,4]) = t^{[1,2,3,4]}, a translation element

3. For two element w1, w2 in \tilde{W},

w1*w2: return the product of w1 and w2

w1==w2: return if w1 equals w2

4. For vectors or rational vectors v1, v2, 

v1.t: the value of v1

v1 + v2: return the sum of v1 and v2

v1 == v2: return if v1 equals v2

scalar(a,v1): return the scalar product of a and v1

