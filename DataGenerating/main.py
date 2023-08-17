# -*- coding: utf-8 -*-
"""
Created on Tue May 21 17:34:51 2019

@author: TOPOLOGY
"""


import numpy as np;  
import joblib
import time
import math
from fractions import Fraction
#########################33

####################   

def to_Bourbaki(nu):
    if TYPE == 'A':
        return(nu)
    if TYPE == 'D':
        return(rational_vect( 1/4* (nu.t - nu.t[::-1])[0:ORDER]))
    if TYPE == 'B' or TYPE == 'C':
        return(rational_vect( 1/2* (nu.t - nu.t[::-1])[0:ORDER]))
    
def from_Bourbaki(v):
    if TYPE == 'A':
        return(vect(v.t))
    if TYPE == 'D':
        nu = rational_vect(e0.t)
        for i in S0key:
            nu=nu+ scalar(  8*(v.t*to_Bourbaki(alpha[i]).t).sum(),omegacheck[i])
        return(vect(nu.t))
    if TYPE == 'B' or 'C':
        nu = rational_vect(e0.t)
        for i in S0key:
            nu=nu+ scalar(  2*(v.t*to_Bourbaki(alpha[i]).t).sum(),omegacheck[i])
        return(vect(nu.t))
def fraction_string(nu):
    nu = nu.t
    fraction_array = []
    for a in nu:
        fraction_array.append(Fraction(a).limit_denominator(1000))
    fraction_str = "[" + ", ".join(str(f.numerator) if f.denominator == 1 else str(f) for f in fraction_array) + "]"
    return(fraction_str)    
            
class vect:
    def __init__(self,t):
        self.t = np.array(t, dtype = int)
    def __add__(a,b):
        return(vect(a.t+b.t))
    def __sub__(a,b):
        return(vect(a.t-b.t))
    def minus(self):
        return(vect(-self.t))
    def __eq__(v1,v2):
        return( (v1.t==v2.t).all() )
    def __hash__(self):
        return(hash(self.t.tobytes()))
    def inner(v1,v2):
        return( (v1.t*v2.t).sum())
    def present(self):
        print(fraction_string(to_Bourbaki(self)))
        
        


class rational_vect:
    def __init__(self,t):
        self.t = np.array(t)
    def __add__(a,b):
        return(rational_vect(a.t+b.t))
    def __sub__(a,b):
        return(rational_vect(a.t-b.t))
    def minus(self):
        return(rational_vect(-self.t))
    def inner(v1,v2):
        return( (v1.t*v2.t).sum())
    def eq(v1,v2):
        nn = True
        for i in range(len(v1.t)):
            if abs(v1.t[i] - v2.t[i]) > 0.01:
                nn = False
                break
        return(nn)
    def present(self):
        print(fraction_string(to_Bourbaki(self)))

def vectisin(v,SET):
    for vv in SET:
        if rational_vect.eq(v,vv):
            return(True)
    return(False)    
def vectunion(S1,S2):
    
    for t in S2:
        if vectisin(t,S1)==(1==0):
            S1=S1|{t}
    return(S1)
    
def vectprint(A):
    for nu in A:
        print(nu.t)

def scalar(a,v):
    return(rational_vect(a*v.t))
def root_expand(v):
    expd = np.zeros(ORDER)
    for i in range(1,ORDER+1):
        expd[i-1] = (v.t * omega[i].t).sum()
    return(expd) 
def is_dominant(lam):
    nn=True
    for al in Phiplus:
        if ((lam.t * al.t).sum())< -0.001:
            nn = False;break
    return(nn)

def dominate(v):
    if TYPE != 'D':
        vt = np.array(v.t)
        vdom=np.array(sorted(vt,reverse = True))
        if type(v) == rational_vect:
            return(rational_vect(vdom))
        if type(v) == vect:
            return(vect(vdom))
    else:
        vt = np.array(v.t)
        vt1 = np.array(v.t)[0:ORDER]
        abs_set=set();minus_set=set()
        have_zero = False; odd_minus = False
        
        for i in range(ORDER):
            abs_set.add(abs(vt1[i]))
            if vt1[i] < -0.001:
                minus_set.add(i)
            odd_minus=(len(minus_set)%2==1)
            
        if odd_minus==False or have_zero == True:
            vdom = np.zeros(ORDER)
            for i in range(ORDER):
                vdom[i] = abs(vt1[i])
            vdom=np.array(sorted(vdom,reverse = True))
            vdom=np.append(vdom,-vdom[::-1])
            if type(v) == rational_vect:
                return(rational_vect(vdom))
            if type(v) == vect:
                return(vect(vdom))
        else:
            for i0 in range(ORDER):
                if abs(vt1[i0]) == min(abs_set):
                    break
            vdom = np.zeros(ORDER)
            for i in range(ORDER):
                if i!=i0:
                    vdom[i] = abs(vt1[i])
                else:
                    vdom[i] = -abs(vt1[i])
            vdom=np.array(sorted(vdom,reverse = True))
            vdom=np.append(vdom,-vdom[::-1])
            if type(v) == rational_vect:
                return(rational_vect(vdom))
            if type(v) == vect:
                return(vect(vdom))   
        
       
#def old_greater_than_or_equal_to(v1,v2):
#    return((root_expand(v1)>=root_expand(v2)).all())

def greater_than_or_equal_to(v1,v2):
    v1=v1.t;v2=v2.t;N=len(v1);a1=0;a2=0
    for i in range(N):
        a1 += v1[i];a2+=v2[i]
        if a1 < a2-0.001:
            return(False)
    return(True)
def height(v):
    h=0
    for i in S0key:
        h=h + (v.t,omega[i].t).sum()
    return(h) 
def breakpoint(nu):
    return(S0key-cent(nu))
class permutation:
    def __init__(self,f):
        self.f = np.array(f, dtype = int)
    def __add__(p1,p2):
        v1 = p1.f
        v2 = p2.f
        v = np.zeros(rank)
        for i in range(rank):
            v[i] = v1[math.floor(v2[i]+0.001)-1]
        return(permutation(v)) 
        
    def inv(self):
        v1 = self.f
        v = np.zeros(rank)
        
        for j in range(rank):
            v[math.floor(v1[j]+0.001)-1] = j + 1
        return(permutation(v))
    
    def __eq__(a,b):
        return( (a.f==b.f).all())
    
    def __hash__(self):
        return(hash(self.f.tobytes()))
        
    
    def order(self):
        if self == Id.f:
            return(1)
        i = 1
        wi = self
        while wi != Id.f:
            i = i + 1
            wi = wi + self 
        return(i) 
        
    def tocycle(p):
        pf = p.f
        c = []
        counter = [True] * rank
        for i in range(rank):
            if counter[i]:
                counter[i] = False
                ci = [i+1]
                j = int(pf[i])
                while j-1 != i:
                    counter[j-1] = False
                    ci.append(j)
                    j = int(pf[j-1])
                c.append(ci)
        return(c)
        
    def topartition(p):
        c = p.tocycle()
        par = []
        for ci in c:
            par.append(len(ci))
        return(tuple(sorted(par,reverse = True)))
        
def act(w,x):
    y = np.zeros(rank)
    for i in range(0,rank):
        ai = int(w.f[i])
        y[ai-1] = x.t[i]
    if type(x)==vect:
        return(vect(y))
    if type(x)==rational_vect:
        return(rational_vect(y))

        
def power(w,r):
    i=1
    ww = w
    while i < r:
        i = i+1
        ww = ww +w
    return(ww)    

def permutationprint(A):
    for a in A:
        print(a.f)

 
def permutation_to_matrix(p):
    
    MM = np.zeros(( rank , rank ))
    for i in range(1,rank + 1):
        MM[p.f[i-1]-1,i-1] = 1
    return( MM  )
    

def matrix_to_permutation(m):
    v = np.zeros(rank)
    for i in range(rank):
        for j in range(rank):
            if m[i,j] == 1:
                v[j] = i + 1
    return(permutation(v))

def cycle_to_permutation(c):
    p = np.zeros(rank)
    for ci in c:
        for i in range(len(ci)-1):
            p[ci[i]-1] = ci[i+1]
        p[ ci[-1]-1 ] = ci[0]
    return(permutation(p))
    
class affweyl:
    def __init__(self,t,f ):
        self.f = f
        self.t = t
    def inv(self):
        return(affweyl(vect(act(self.f.inv(), (self.t).minus()).t), self.f.inv()))
    
    def __eq__(self,other):
        return( self.t == other.t and self.f == other.f )     
    def __hash__(self):
        return(hash( (self.t,self.f) ) )
    def __add__(self,other):
        return(affweyl(self.t + act(self.f,other.t),self.f + other.f ))   
    def __mul__(self,other):
        return(self + other)
    
    def translation(self):
        return(affweyl(self.t,Id.f))
    
    def finite(self):
        return(affweyl(e0,self.f))
        
    def len(self):
        l=0
    
        for al in Phiplus:
            
            if vectisin(act(self.f.inv(),al),Phiplus) :
                
                l += abs((self.t.t*al.t).sum())   
                # omit dualroot
            else:
                l = l + abs((self.t.t*al.t).sum()-1)
            #print(l)    
        return(math.floor(l+0.001)) 
        
    def towp(w):
        c = w.f.tocycle()
        wp = []
        for ci in c:
            sumi = 0
            for j in ci:
                sumi += w.t.t[j-1]
            wp.append((len(ci),sumi))
        return(tuple(sorted(wp,reverse = True)))
    
    def Newton(self):
        if len(INPUT)==2 and TYPE != 'D':
            wp = self.towp()
            nuw = np.array([])
            for bc in wp:
                nuw = np.append(nuw, bc[1]/bc[0]*np.ones(bc[0]))
                
            return( dominate(rational_vect(nuw)))
        else:
            if len(INPUT)>2:
                u = affweyl(e0,self.f)
                if u == Id:
                    return(dominate((diamond(self.t))))
                else:
                    n_w = 2 * (u.f + sigma(u).f).order()
                    nu0 = self.t
                    nu1 = self.t
                    for i in range(1,n_w):
                        nu1 = act(u.f,sigma(nu1))
                        nu0 = nu0 + nu1
        
                return(dominate( scalar( 1/n_w ,nu0)  ))
            else:
                u = affweyl(e0,self.f)
                if u == Id:
                    return(dominate(self.t))
                else:
                    n_w =  (u.f).order()
                    nu0 = self.t
                    nu1 = self.t
                    for i in range(1,n_w):
                        nu1 = act(u.f, nu1 )
                        nu0 = nu0 + nu1
        
                return(dominate( scalar( 1/n_w ,nu0)  ))
    
    
    def isstraight(self):
        return(  abs(self.len() - (self.Newton().t*two_rho.t).sum())<0.01)
    def rex(self):
        rex = np.array([],dtype=int)  
        ww = self
        while ww != Id:
            for i in Skey:
                if isdecrease(ww,i):
                    rex = np.append(rex,i)
                    ww = ww + s[i]
                    break
        return(rex[::-1])
 
    def present(self):
        if TYPE == 'A':
            v = (fraction_string(self.t))
        if TYPE == 'D':
            v = rational_vect( 1/4* (self.t.t - self.t.t[::-1])[0:ORDER])
            v = (fraction_string(v))
        if TYPE == 'B' or TYPE == 'C':
            v = rational_vect( 1/2* (self.t.t - self.t.t[::-1])[0:ORDER])
            v = (fraction_string(v))
        print(f'translation part: {v}')
        print(f'finite part:{affweyl(e0,self.f).rex()}') 

def trans(v):
    return(affweyl(vect(v),Id.f))
        
def affine_Weyl(v,e):
    return( affweyl(from_Bourbaki(rational_vect(v)),Id.f)+exp(e))

def affweylprint(A):
    for w in A:
        print(w.rex())
def affweylpower(w,r):
    if r<0:
        print('error')
        return(None)
    if r==0:
        return(Id)
    i = 1
    ww = w
    while i < r:
        i = i + 1
        ww = ww + w
    return(ww)    

def isdecrease(w,j): # if w + s[j] < w
    mu = w.t
    u = w.f
    if j != 0:
        
        uj = act(u,alpha[j]).minus()
        if vectisin(uj,Phiplus):
            
            return(  (mu.t*uj.t).sum() < 1)
        else:
            
            return( (mu.t*uj.t).sum() < 0)
    else:
        utheta = act(u,theta)
        if vectisin(utheta,Phiplus):
            
            return( 1 +  (mu.t*utheta.t).sum() < 1)
        else:
            
            return( 1 +  (mu.t*utheta.t).sum() < 0)
    
    
def wp_to_st(wp):
    mu = np.array([])
    cycle = []
    j = 1
    for bc in wp:
        b = bc[0]
        c = bc[1]
        cycle.append(  list(range(j,b+j))  )
        j += b
        mub = np.zeros(b)
        for i in range(b):
            mub[i] = math.ceil( (i+1)*c/b) - math.ceil( (i)  *c/b)
        mu = np.append(mu,mub[::-1])    
    p = cycle_to_permutation(cycle)
    return(affweyl(vect(mu),p))
    
def a_seq(i,w):
    mu = w.t.t
    u = w.f
    N = 0
    for j in W0:
        N += len(W0[j])
    
    ai = np.zeros(N)
    uinvk = u
    for k in range(N):
        uinvk = uinvk + u.inv()
        ai[k] =  mu[  (uinvk.f[i-1]) -1]
    return(list(ai))

# reduce standard wst to in W_I(w)w
def minf(wp):
    u = Id
    w = wp_to_st(wp)
    aiw = []
    for i in range(rank):
        aiw.append(a_seq(i+1,w))
    for i in range(2,rank+1):
        i = rank + 1 - i
        for j in range(1,i +1):
            if aiw[j-1] < aiw[j]:
                temp = aiw[j]
                aiw[j] = aiw[j-1]
                aiw[j-1] = temp
                u = s[j] + u
    return(u + w + u.inv())



    
                     
def exp(exp):
    w = Id
    for i in range(0,len(exp)):
        u = exp[i]
        w = w + s[u]
        
    return(w)
    

def star(a,b):
    
    rexb = b.rex()
    for i in range(0,len(rexb)):
        if isdecrease(a,rexb[i])==False:
            a = a + s[rexb[i]]
    return(a)
def Demazure_operator(a,b):
    
    rexb=b.rex()
    for i in range(0,len(rexb)):
        if isdecrease(a,rexb[i]): 
            a=a+s[rexb[i]]
    return(a)           
class parabolic:
    def __init__(self,J,y):
        self.J=J
        self.y=y
    def __eq__(self,other):
        return(self.J==other.J and self.y == other.y)
    def __hash__(self):
        return(hash((self.J,self.y)))
        
##################################    partial conjugation
def ldecom(w,I):
    if I == set():
        return(     (Id,w)     )
    if len(I) == ORDER:
        return((w,Id))
    w1 = []
    w2 = w
    nn = True
    while nn == True:
        nn = False
        for j in I:
            
            if isdecrease(w2.inv(),j) :
                w1.append(s[j])
                w2 = s[j] + w2
                nn = True
                break
        
    a = Id
    for t in w1:
        a = a + t
    return(a,w2)

def cent(mu):
    I = set()
    for j in S0key:
        if abs((alpha[j].t*mu.t).sum())<0.01:
            I.add(j)
    return(I)                

# only work for type A, sorted         
def partial(w):
    wf = affweyl(e0,w.f)
    for i in W0:
        for x in W0[i]:
            xwt = act(x.f,w.t)
            if is_dominant(xwt ):        
                de = ldecom( x + wf, cent(xwt))
                return((x.inv()+de[0]), xwt, de[1]  )
            
def partialprint(w):
    new = partial(w)          
    print(new[0].rex(),new[1].t ,new[2].rex())

def rdecom(w,I):
    
    if I == set():
        return(     (w,Id)     )
    for j in I:
        if isdecrease(w,j):
            
            return(  ((rdecom(w + s[j],I)[0]),(rdecom(w + s[j],I)[1]) + s[j]) )
    return((w,Id) )    

       
def ldecomprint(w,I):
    print(ldecom(w,I)[0].rex(),ldecom(w,I)[1].rex())

def rdecomprint(w,I):
    print(rdecom(w,I)[0].rex(),rdecom(w,I)[1].rex())


def support(w):
    rex = w.rex()
    supp = set()
    for i in range(0,len(rex)):
        supp.add(rex[i])
    return(supp)

###################################################################        
    
#####################################################################   
#########################################################  DATA input
############################################  
available_type = {'A2','A3','A4','A5',
                  'B2','B3','B4','B5',
                  'C2','C3','C4','C5', 
                  'D4','D5',
                  '2A2','2A3','2A4','2A5',
                  '2D4','2D5'}

INPUT = input('please input type of the group (An, 2An, Bn, Cn, Dn, 2Dn): ') 
if INPUT[0] not in {'2','3'} and len(INPUT) >= 2:
    TYPE = INPUT[0]; ORDER = int(INPUT[1:])
else:
    TYPE = INPUT[1]; ORDER = int(INPUT[2:])
if TYPE == 'A':
    rank = ORDER + 1
    semisimplerank = ORDER
    Skey =set(range(0,ORDER+1))
    S0key =Skey-{0}
    e0 = vect((np.zeros(rank)))
    alpha=dict();alphacheck=dict();omega=dict();omegacheck=dict();
    s = dict()
    for j in range(1,rank):
        v = np.zeros(rank); v[j-1] = 1; v[j] = -1
        alpha[j] = rational_vect(v)
        alphacheck[j] = vect(alpha[j].t)
    for j in range(1,rank):
        v = np.zeros(rank)
        for k in range(1,j+1):
            v[k-1] = 1
        omega[j] = rational_vect(v)
        omegacheck[j] = vect(v)
    omega[rank] = rational_vect(np.ones(rank))
    for j in range(1,rank):
        f = np.array(range(1,rank+1))
        f[j-1] = j+1; f[j] = j;
        s[j] = affweyl(e0, permutation(f))
    S0 = set()
    for i in S0key:
        S0 = S0|{s[i]}
    v = np.zeros(rank); v[0] = 1; v[-1] = -1
    theta = rational_vect(v); thetacheck = vect(v)
    f=np.array(range(1,rank+1)); f[0]=rank; f[rank-1]=1
    s[0] = affweyl(thetacheck,permutation(f))
    S={s[0]}|S0
    Phiplus = set()
    for i in range(1,rank):
        for j in range(i+1, rank+1):
            v = np.zeros(rank);v[i-1]=1;v[j-1]=-1
            Phiplus=vectunion(Phiplus,{rational_vect(v)})
    
    Id = affweyl(e0, permutation((np.arange(1,rank+1))))
    two_rho = rational_vect(e0.t)            
    for al in Phiplus:
        two_rho = two_rho + al
    rho = scalar(1/2,two_rho)
    w_0 = affweyl(e0,permutation(np.array(range(1,rank+1)[::-1])))
    tau = dict(); tau[0]=Id;
    for j in range(1,rank):
        fj1=np.array(range(1,j+1))[::-1]
        fj2=np.array(range(j+1,rank+1))[::-1]
        fj=np.append(fj1,fj2)
        tau[j] = affweyl(omegacheck[j],Id.f)+affweyl(e0,permutation(fj))+w_0    
    Omega = set(tau.values())

if TYPE == 'B':
    rank = 2*ORDER + 1
    semisimplerank = ORDER
    Skey =set(range(0,ORDER+1))
    S0key =Skey-{0}
    e0 = vect((np.zeros(rank)))

    alpha=dict();alphacheck=dict();omega=dict();omegacheck=dict();
    s = dict()
    for j in range(1,ORDER):
        v=np.zeros(rank);v[j-1]=1;v[j]=-1;v[rank-j-1]=1;v[rank-j]=-1
        alpha[j] = rational_vect(1/2*v)
        alphacheck[j] = vect(v)
    v=np.zeros(rank);v[ORDER-1]=2;v[ORDER+1]=-2
    alphacheck[ORDER]=vect(v);alpha[ORDER]=rational_vect(1/4*v)
    
    for j in range(1,ORDER+1):
        v = np.zeros(rank)
        for k in range(1,j+1):
            v[k-1]=1;v[rank-k]=-1
        omega[j] = rational_vect(1/2*v)
        omegacheck[j] = vect(v)
    omega[ORDER]=rational_vect(1/2*np.append(np.ones(ORDER),np.zeros(ORDER+1)))
    
    
    
    for j in range(1,ORDER+1):
        f = np.array(range(1,rank+1))
        if j < ORDER:
            f[j-1]=j+1;f[j]=j;f[rank-j-1]=rank-j+1;f[rank-j]=rank-j
        if j == ORDER:
            f[j-1]=j+2;f[j+1]=j;
        s[j] = affweyl(e0, permutation(f))
    S0 = set()
    for i in S0key:
        S0 = S0|{s[i]}
        
    v = np.zeros(rank);v[0]=1;v[-1]=-1;v[1]=1;v[-2]=-1
    theta=rational_vect(1/2*v);thetacheck=vect(v)
    f=np.array(range(1,rank+1));f[0]=rank-1;f[1]=rank;f[-2]=1;f[-1]=2
    s[0] = affweyl(thetacheck,permutation(f))
    S={s[0]}|S0
    
    Phiplus = set()
    for i in range(1,ORDER):
        for j in range(i+1, ORDER+1):
            v=np.zeros(rank);
            v[i-1]=1/2;v[j-1]=-1/2
            v[rank-j]=1/2;v[rank-i]=-1/2
            Phiplus=vectunion(Phiplus,{rational_vect(v)})
            v=np.zeros(rank);
            v[i-1]=1/2;v[j-1]=1/2
            v[rank-j]=-1/2;v[rank-i]=-1/2
            Phiplus=vectunion(Phiplus,{rational_vect(v)})
    for k in range(1,ORDER+1):
        v=np.zeros(rank);v[k-1]=1/2;v[rank-k]=-1/2
        Phiplus=vectunion(Phiplus,{rational_vect(v)})
    
    Id = affweyl(e0, permutation((np.arange(1,rank+1))))
    two_rho = rational_vect(e0.t)            
    for al in Phiplus:
        two_rho = two_rho + al
    rho = scalar(1/2,two_rho)
    w_0 = affweyl(e0,permutation(np.array(range(1,rank+1)[::-1])))
    
    tau = dict(); tau[0]=Id
    f=np.array(range(1,rank+1));f[0]=rank;f[-1]=1
    tau[1] = affweyl(omegacheck[1],Id.f)+affweyl(e0,permutation(f))    
    Omega = set(tau.values())   

if TYPE == 'C':
    rank = 2*ORDER
    semisimplerank = ORDER
    
    Skey =set(range(0,ORDER+1))
    S0key =Skey-{0}

    e0 = vect((np.zeros(rank)))
    alpha=dict();alphacheck=dict();omega=dict();omegacheck=dict();
    s = dict()
    for j in range(1,ORDER):
        v=np.zeros(rank);v[j-1]=1;v[j]=-1;v[rank-j-1]=1;v[rank-j]=-1
        alpha[j] = rational_vect(1/2*v)
        alphacheck[j] = vect(v)
    v=np.zeros(rank);v[ORDER-1]=1;v[ORDER]=-1
    alphacheck[ORDER]=vect(v);alpha[ORDER]=rational_vect(v)
    
    for j in range(1,ORDER):
        v = np.zeros(rank)
        for k in range(1,j+1):
            v[k-1]=1;v[rank-k]=-1
        omega[j] = rational_vect(1/2*v)
        omegacheck[j] = vect(v)
    v0=np.zeros(ORDER);v1=np.ones(ORDER)
    omega[ORDER]=rational_vect(1/2*np.append(v1,-v1))
    omegacheck[ORDER]=vect(np.append(v1,v0))
    
    
    for j in range(1,ORDER+1):
        f = np.array(range(1,rank+1))
        if j < ORDER:
            f[j-1]=j+1;f[j]=j;f[rank-j-1]=rank-j+1;f[rank-j]=rank-j
        if j == ORDER:
            f[j-1]=j+1;f[j]=j;
        s[j] = affweyl(e0, permutation(f))
    S0 = set()
    for i in S0key:
        S0 = S0|{s[i]}
        
    v = np.zeros(rank);v[0]=1;v[-1]=-1 
    theta=rational_vect( v);thetacheck=vect(v)
    f=np.array(range(1,rank+1));f[0]=rank;f[-1]=1
    s[0] = affweyl(thetacheck,permutation(f))
    S={s[0]}|S0
    
    Phiplus = set()
    for i in range(1,ORDER):
        for j in range(i+1, ORDER+1):
            v=np.zeros(rank);
            v[i-1]=1/2;v[j-1]=-1/2
            v[rank-j]=1/2;v[rank-i]=-1/2
            Phiplus=vectunion(Phiplus,{rational_vect(v)})
            v=np.zeros(rank);
            v[i-1]=1/2;v[j-1]=1/2
            v[rank-j]=-1/2;v[rank-i]=-1/2
            Phiplus=vectunion(Phiplus,{rational_vect(v)})
    for k in range(1,ORDER+1):
        v=np.zeros(rank);v[k-1]=1;v[rank-k]=-1
        Phiplus=vectunion(Phiplus,{rational_vect(v)})
    
    Id = affweyl(e0, permutation((np.arange(1,rank+1))))

    two_rho = rational_vect(e0.t)            
    for al in Phiplus:
        two_rho = two_rho + al
    rho = scalar(1/2,two_rho)
    w_0 = affweyl(e0,permutation(np.array(range(1,rank+1)[::-1])))
    
    tau = dict(); tau[0]=Id
    f1=np.array(range(1,ORDER+1))
    f2=np.array(range(ORDER+1,rank+1))
    f = np.append(f2,f1)
    tau[ORDER] = affweyl(omegacheck[ORDER],Id.f)+affweyl(e0,permutation(f))    
    Omega = set(tau.values())   

if TYPE == 'D':
    rank = 2*ORDER
    semisimplerank = ORDER
    Skey =set(range(0,ORDER+1))
    S0key =Skey-{0}
    
    e0 = vect((np.zeros(rank)))
    alpha=dict();alphacheck=dict();omega=dict();omegacheck=dict();
    s = dict()
    for j in range(1,ORDER):
        v=np.zeros(rank);v[j-1]=1;v[j]=-1;v[rank-j-1]=1;v[rank-j]=-1
        alpha[j] = rational_vect(1/4*v)
        alphacheck[j] = vect(2*v)
    v=np.zeros(rank);v[ORDER-2]=1;v[ORDER-1]=1
    v[ORDER]=-1;v[ORDER+1]=-1
    alphacheck[ORDER]=vect(2*v);alpha[ORDER]=rational_vect(1/4*v)
    
    for j in range(1,ORDER-1):
        v = np.zeros(rank)
        for k in range(1,j+1):
            v[k-1]=1;v[rank-k]=-1
        omega[j] = rational_vect(1/4*v)
        omegacheck[j] = vect(2*v)
    v = np.zeros(rank)
    for k in range(1,ORDER+1):
        v[k-1]=1;v[rank-k]=-1
    omega[ORDER]=rational_vect(1/8*v)
    omegacheck[ORDER]=vect(v)
    v[ORDER-1]=-1;v[ORDER]=1
    omega[ORDER-1]=rational_vect(1/8*v)
    omegacheck[ORDER-1]=vect(v)
    
    
    
    for j in range(1,ORDER+1):
        f = np.array(range(1,rank+1))
        if j < ORDER:
            f[j-1]=j+1;f[j]=j;f[rank-j-1]=rank-j+1;f[rank-j]=rank-j
        if j == ORDER:
            f[j-2]=j+1;f[j-1]=j+2;f[j]=j-1;f[j+1]=j
        s[j] = affweyl(e0, permutation(f))
    S0 = set()
    for i in S0key:
        S0 = S0|{s[i]}
        
    v = np.zeros(rank);v[0]=1;v[-1]=-1;v[1]=1;v[-2]=-1
    theta=rational_vect(1/4*v);thetacheck=vect(2*v)
    f=np.array(range(1,rank+1));f[0]=rank-1;f[1]=rank;f[-1]=2;f[-2]=1
    s[0] = affweyl(thetacheck,permutation(f))
    S={s[0]}|S0
    
    Phiplus = set()
    for i in range(1,ORDER):
        for j in range(i+1, ORDER+1):
            v=np.zeros(rank);
            v[i-1]=1/4;v[j-1]=-1/4
            v[rank-j]=1/4;v[rank-i]=-1/4
            Phiplus=vectunion(Phiplus,{rational_vect(v)})
            v=np.zeros(rank);
            v[i-1]=1/4;v[j-1]=1/4
            v[rank-j]=-1/4;v[rank-i]=-1/4
            Phiplus=vectunion(Phiplus,{rational_vect(v)})
    
    Id = affweyl(e0, permutation((np.arange(1,rank+1))))
    two_rho = rational_vect(e0.t)            
    for al in Phiplus:
        two_rho = two_rho + al
    rho = scalar(1/2,two_rho)
    
    tau = dict(); tau[0]=Id
    f=np.array(range(1,rank+1))
    f[0]=rank;f[-1]=1
    f[ORDER-1]=ORDER+1
    f[ORDER]=ORDER
    tau[1]=affweyl(omegacheck[1],Id.f)+affweyl(e0,permutation(f))    
    f1=np.array(range(1,ORDER+1))
    f2=np.array(range(ORDER+1,rank+1))
    if ORDER%2==0:
        w_0 = affweyl(e0,permutation(np.array(range(1,rank+1)[::-1])))
        f = np.append(f2,f1)
        tau[ORDER]=trans(omegacheck[ORDER].t)+affweyl(e0,permutation(f)) 
        ff = np.array(f)
        ff[0]=ORDER;ff[-1]=ORDER+1;
        ff[ORDER-1]=1;ff[ORDER]=rank
        tau[ORDER-1]=trans(omegacheck[ORDER-1].t)+affweyl(e0,permutation(ff)) 
    if ORDER%2!=0:
        w0f=np.array(range(1,rank+1)[::-1])
        w0f[ORDER-1]=ORDER;w0f[ORDER]=ORDER+1
        w_0 = affweyl(e0,permutation(w0f))
        f = np.append(f2,f1)
        f[ORDER-1]=1;f[ORDER]=rank
        tau[ORDER]=trans(omegacheck[ORDER].t)+affweyl(e0,permutation(f)) 
        ff = np.append(f2,f1)
        ff[0]=ORDER;ff[-1]=ORDER+1;
        tau[ORDER-1]=trans(omegacheck[ORDER-1].t)+affweyl(e0,permutation(ff)) 
        
    Omega = set(tau.values())   

def diamond(mu):
    if len(INPUT) == 2:
        return(mu)        
    if len(INPUT)>2:
        return(scalar(1/2,rational_vect(mu.t + sigma(mu).t)))  

def sigma(a):
    if len(INPUT) == 2:
        return(a)        
    if len(INPUT)>2:
    
        if TYPE == 'A':
            if type(a) == vect:
                return(vect( -a.t[::-1]   ))
            if type(a) == rational_vect:
                return(rational_vect( -a.t[::-1]   ))
            if type(a) == affweyl:
                at = a.translation(); af = a.finite()
                sigat = trans(sigma(at.t).t)
                rexaf=af.rex()
                sigaf=Id
                for i in range(0,len(rexaf)):
                    sigaf=sigaf+s[rank - rexaf[i]]
                return(sigat + sigaf)
        if TYPE == 'D':
            if type(a) == vect:
                at = np.array(a.t); at[ORDER]=-at[ORDER]
                at[ORDER-1]=-at[ORDER-1]
                return(vect(at))
            if type(a) == rational_vect:
                at = np.array(a.t); at[ORDER]=-at[ORDER]
                at[ORDER-1]=-at[ORDER-1]
                return(rational_vect(at))
            if type(a) == affweyl:
                at = a.translation(); af = a.finite()
                sigat = trans(sigma(at.t).t)
                rexaf=af.rex()
                sigaf=Id
                for i in range(0,len(rexaf)):
                    if rexaf[i] < ORDER -1:
                        sigaf=sigaf+s[rexaf[i]]
                    if rexaf[i] == ORDER-1:
                        sigaf=sigaf+s[ORDER]
                    if rexaf[i] == ORDER:
                        sigaf=sigaf+s[ORDER-1]
                return(sigat + sigaf)
            

####################################################################
##########  DL datum
class DLdatum():
    def __init__(self,w,n):
        self.w=w
        self.n=n
    #def __eq__(self,other):
    #    return(self.w == other.w and self.n==other.n)
    #def __hash__(self):
    #    return(hash( (self.w,self.n)  ))

            
  
    


################### defect
def nu_defect(nu):
    defectnu = 0
    for i in range(1,ORDER+1):
        ai =  (omega[i].t * nu.t).sum() + 0.001
        defectnu += ai - math.floor(ai)
    return(math.floor(2*defectnu+0.001))
def length(nu,mu):
    l = 0; mu_nu = rational_vect(mu.t-nu.t)
    for i in range(1,ORDER+1):
        ai =  vect.inner(omega[i],mu_nu)-0.001
        l +=  math.ceil(ai)
    return(math.floor(l+0.001))
def sigma_length(nu,mu):
    l = 0; mu_nu = mu.t-nu.t
    for omega_c in omega_sigma_set:
        l+=math.ceil( (mu_nu*omega_c.t).sum()-0.001)
    return(l)      


def defect(b):
    nu = (b.Newton())
    return(nu_defect(nu))
         
#def defect(b):
#    return(np.linalg.matrix_rank(permutation_to_matrix(Id.f) - permutation_to_matrix(b.f) ))
 

def eta_delta(w):
    parw = partial(w)
    return( parw[2] + parw[0])




def virtualdim(w,nu): 
    return(1/2*(w.len()+eta_delta(w).len()-nu_defect(nu)- (two_rho.t*nu.t).sum()))


def is_Coxeter(w):
    return(support(w) == S0key and w.len() == ORDER)
       
def is_partial_Coxeter(w):
    return(w.len() == len(support(w)))
       
############ Bruhat          
def Bruhat(x,y):
    if x.len() > y.len():
        return('Not <=')
    if x == y:
        return('<=')
    for j in Skey:
        if isdecrease(y.inv(),j):
            break
    if isdecrease(x.inv(),j):    
        return( Bruhat(s[j]+x,s[j]+y) )
    
    if isdecrease(x.inv(),j) == False: 

        return( Bruhat(x,s[j]+y) )        


def floor(nu):
    nu=nu.t;lam = np.zeros(rank)
    lam[0] = math.floor(nu[0])
    for i in range(1,rank):
        lam[i] = math.floor((omega[i+1].t * nu).sum() - (omega[i].t * lam).sum() +0.0001)
    return(rational_vect(lam))
##################    
def tauprint(a):
    print(  (a+tau.inv()).rex())

#############################   reduction for split case, tau may not be Id



def one_step(w):
    if INPUT in available_type and INPUT[0]=='A' and isminimal(w):
        return(None)
    O=set()
    Onew=set()
    Onewnew={w}
    while Onewnew != set():
        O = O|Onew
        Onew = Onewnew
        Onewnew=set()
        for w in Onew:
            lw = w.len()
            for j in Skey:
                jw = s[j] + w
                jwj = s[j] + w + sigma(s[j]); ljwj = jwj.len()
                if ljwj < lw:
                    return((jw, jwj))
                if ljwj == lw:
                    if jwj not in O and jwj not in Onew:
                        Onewnew.add(jwj)
    return (None)     
    

def numeration(w):
    O=set()
    Onew=set()
    Onewnew={w}
    while Onewnew!=set():
        O=O|Onew
        Onew=Onewnew
        Onewnew=set()
        for w in Onew:
            lw = w.len()
            for j in Skey:
                jwj = s[j]+w+sigma(s[j])
                if jwj.len()==lw:
                    if jwj not in O and jwj not in Onew:
                        Onewnew.add(jwj)
    # Omega conjugation    
    O = O|Onew 
    OO = set()
    for w in O:
        OO.add(w)
    for j in tau:
        for w in O:
            OO.add(tau[j]+w+sigma(tau[j]).inv())
    return(OO)

def tominimal(w):
    while True:
        onew = one_step(w)
        if onew == None:
            return(w)
            
        else:
            w = onew[1]
              

def reduction(w):
    P = set()
    Pnew = {DLdatum(w,0)}
    while len(Pnew) != len(P):
        P = Pnew
        Pnew = set()        
        for p in P :
            repw = one_step(p.w)
            if repw is None:
                Pnew.add(p)
            else:
                Pnew = Pnew | {DLdatum(repw[0],p.n+1),DLdatum(repw[1],p.n+1)}
    return(  Pnew    )
    

def reduction_step(w,d):
    if d<=0:
        return({DLdatum(w,0)})
    else:
        P = set(); r=0; Pnew = {DLdatum(w,0)}
        print(f'r={r}')
        while r<d and len(Pnew) != len(P):
            print(r)
            P = Pnew;
            Pnew = set()        
            for p in P :
                repw = one_step(p.w)
                if repw is None:
                    Pnew.add(p)
                else:
                    Pnew = Pnew | {DLdatum(repw[0],p.n+1),DLdatum(repw[1],p.n+1)}
            r+=1    
        return(  Pnew    )
    
def dim_dictw_dict_to30(w,nu):
    lw = w.len(); nu=(rational_vect(nu))
    dimset=set(); nurho = (nu.t*two_rho.t).sum()
    for dat in reduction_step(w,lw-29):
        datw = dat.w
        if datw.len()<30:
            dictdatw = dictw_dict_to30[datw]
            for nuw in dictdatw:
                if rational_vect.eq( nuw,nu):
                    dimset.add(dictdatw[nuw][0] - nurho + dat.n)
                    break
        else:
            if rational_vect.eq( to_Bourbaki(datw.Newton()),nu):
                dimset.add(datw.len() -nurho + dat.n)
    if dimset==set():
        return('empty')
    return(math.floor(max(dimset)+0.001))
    
    
def partial_reduction(w):
    P = set()
    Pnew = {DLdatum(w,0)}
    while len(Pnew) != len(P):
        P = Pnew
        Pnew = set()        
        for p in P :
            w = p.w; n = p.n; lw = w.len(); nn = True
            for j in S0key:
                jw = s[j] + w
                if jw.len() < w.len():
                    jwj = jw + s[j]
                    if jwj.len() < lw:
                        Pnew.add(DLdatum(jw,n+1));Pnew.add(DLdatum(jwj,n+1))
                        nn = False
                        break
                    else:
                        Pnew.add(DLdatum(jwj,n))
                        nn = False
                        break
            if nn:
                Pnew.add(p)            
                
    return(  Pnew    ) 
def reductionprint(w):
    #tau = tau2
    for dat in reduction(w):
        print(f'dat.w.rex() = {dat.w.rex()}',f'dat.n = {dat.n}' )
        print(f'dat.w.Newton() = {dat.w.Newton().t}' )
        print(f'dat.n + dat.w.len() = {dat.n + dat.w.len()}')
        print('\n')  


def dim(w,nu):
    dimset=set();nu=(rational_vect(nu))
    for dat in reduction(w):
        nudat = dat.w.Newton()
        if rational_vect.eq( to_Bourbaki(nudat),nu):
            nurho = (nudat.t*two_rho.t).sum()
            dimset.add(dat.w.len() -nurho + dat.n)
            
    if dimset==set():
        return('empty')
    return(math.floor(max(dimset)+0.001))  
def Idim(w,nu):
    dimset=set();nu=(rational_vect(nu))
    for dat in reduction(w):
        if rational_vect.eq((dat.w.Newton()),nu):
            dimset.add(dat.w.len() + dat.n)        
    if dimset==set():
        return('empty')
    return(math.floor(max(dimset)+0.001))  
def irr(w, nu):
    redw = reduction(w);dimw = set();nu=(rational_vect(nu))
    nurho = (from_Bourbaki(nu).t*two_rho.t).sum()
    for dat in redw:
        if rational_vect.eq(to_Bourbaki(dat.w.Newton()),nu):
            dimw.add(dat.w.len()-nurho+dat.n)
    if dimw==set():
        return(0)
    else:
        dimw = math.floor(max(dimw)+0.001); lc = 0
        for dat in redw:
            if rational_vect.eq(to_Bourbaki(dat.w.Newton()),nu):
                if math.floor(dat.w.len()-nurho+dat.n+0.001)==dimw:
                    lc+=1
        return(lc)  
    

def dim_irr_dict(w):
    if type(w)==affweyl:
        redw = reduction(w);dimdict = dict();irrdict = dict()
        BGw = set()
        for dat in redw:
             BGw = vectunion(BGw,{ dat.w.Newton() })
        for nu in BGw:
            dimdict[nu] = set(); irrdict[nu] = 0
        for dat in redw:
            nudat = dat.w.Newton()
            nutworho = (nudat.t*two_rho.t).sum()
            for nu in dimdict:
                if rational_vect.eq(nu,nudat):
                    dimdict[nu].add(dat.w.len() + dat.n - nutworho)
                    break
   
        for nu in dimdict:
            dimdict[nu] = math.floor(max(dimdict[nu])+0.001)
        for dat in redw:
            nudat = dat.w.Newton()
            nutworho = (nudat.t*two_rho.t).sum()
            for nu in dimdict:
                if rational_vect.eq(nu,nudat):
                    if abs((dat.w.len() -nutworho +dat.n)-dimdict[nu])<0.01:
                        irrdict[nu] += 1
                        break
                        
        dictw = dict()
        for nu in BGw:
            dictw[nu] = (dimdict[nu],irrdict[nu])
        return(dictw)
    if type(w)==dict:
        nuset=set();dimdict=dict();irrdict=dict()
        alldict=dict()
        for i in w:
            for a in w[i]:
                alldict[a] = dim_irr_dict(a)
                for nu in alldict[a]:
                    nuset=vectunion(nuset,{nu})
        for nu in nuset:
            dimdict[nu]=set();irrdict[nu]=0
        for a in alldict:
            dicta = alldict[a]
            for nu in dicta:
                for nunu in dimdict:
                    if rational_vect.eq(nu,nunu):
                        dimdict[nunu].add(dicta[nu][0])
        for nu in dimdict:
            dimdict[nu]=math.floor(max(dimdict[nu])+0.001)
        for a in alldict:
            dicta = alldict[a]
            for nu in dicta:
                for nunu in irrdict:
                    if rational_vect.eq(nu,nunu):
                        if dimdict[nunu]==dicta[nu][0]:
                            irrdict[nunu]+=dicta[nu][1]
        result = dict()
        for nu in nuset:
            result[nu] = (dimdict[nu],irrdict[nu])
        return(result)
            




def dim_irr_print(w):
    dictw = dim_irr_dict(w)
    for nu in dictw:
        print(f'Newton point = {fraction_string(to_Bourbaki(nu))},', f'dim = {dictw[nu][0]},', f'irr = {dictw[nu][1]}' )






########################################################
# cordial element

    
def max_Newton(w):
    one = one_step(w)
    if one is None:
        return (w.Newton())
    else:
        return (max_Newton(one[0] ))


def BG(w):
    BGw = set()
    redw = reduction(w)
    for dat in redw:
        BGw = vectunion(BGw,{dat.w.Newton()})
    return (BGw)



def is_cordial(w):
    nu = (max_Newton(w))
    return(  abs(w.len()- eta_delta(w).len() - (nu.t*two_rho.t).sum() + nu_defect(nu))<0.01)

def is_semicordial(w):
    numax = (max_Newton(w))
    D = virtualdim(w,numax) - w.len()+(numax.t*two_rho.t).sum()            
    BGw = BG(w)                      
    for nu0 in BGw:
        nn = True
        for nu in BGw:
            if rational_vect.eq(nu0,nu)==False:
                if greater_than_or_equal_to(nu0, nu):
                    nn = False
        if nn:
            break
    D0 = virtualdim(w, nu0)-dim(w,to_Bourbaki(nu0).t)
    if abs(D-D0)<0.01:
        return(True)
    else:
        #print(nu0.t)
        return(False)
    
    

    
###########################################
# W = joblib.load('A2_W(to30)')
####################### compute X(mu,b)
def give_tau(mu):
    if TYPE == 'A':
        N = mu.t.sum(); M = N//(ORDER+1)
        i = N % (ORDER+1)
        return(tau[i] + trans(M*np.ones(ORDER+1)))
    
    if TYPE == 'B':
        v = rational_vect( 1/2* (mu.t - mu.t[::-1])[0:ORDER])
        N = math.floor(v.t.sum()+0.001)
        if N%2==0:
            return(Id)
        if N%2==1:
            return(tau[1])
                    
    if TYPE == 'C':
        
        N= mu.t.sum(); r=N//ORDER
        return(affweylpower(tau[ORDER],r))
        
        
        
    if TYPE == 'D':
        v = rational_vect( 1/4* (mu.t - mu.t[::-1])[0:ORDER])
        N = math.floor(v.t.sum()+0.001); nn = True
        for a in v.t:
            if (v.t[0]  - math.floor(v.t[0]+0.001))> 0.001:
                nn = False;break
        if nn:
            if N%2 == 0:
                return(Id)
            else:
                return(tau[1])
        else:
            vORDER=1/4*(omegacheck[ORDER].t-omegacheck[ORDER].t[::-1])[0:ORDER]
            NORDER=vORDER.sum()
            if ORDER%2==0:
                if N%2 == NORDER%2 :
                    return(tau[ORDER])
                else:
                    return(tau[ORDER-1])
            else:
                if (int(N-1/2))%2 == (int(NORDER-1/2))%2 :
                    return(tau[ORDER])
                else:
                    return(tau[ORDER-1])
    
def Adm(mu):
    mu = from_Bourbaki(rational_vect(mu))
    return(Adm_old(mu,give_tau(mu)))
    
def Adm_old(mu,tau): #  w + tau <= trans(mu), w in Wa
    N = math.floor((mu.t*two_rho.t).sum()+0.001)
    Admmu = dict()
    for i in range(N+1):
        Admmu[i] = set()
    
    for i in W0:
        for u in W0[i]:
            Admmu[N].add(u + trans(mu.t) + u.inv())
    
    for i in range(1,N+1 ):
        j = N - i
        #print(j)
        for w in Admmu[j+1]:
            wrex = (w+tau.inv()).rex()
            for k in range(j+1):
                wrexk = np.delete(wrex, k)
                wk = exp(wrexk)+tau;lwk=wk.len()
                Admmu[lwk].add(wk)
    return(Admmu)
def is_decom(mu,nu):
    for i in S0key:
        if ((mu.t-nu.t)*omega[i].t).sum()<0.001:
            return(True)
    return(False)
    
def Admisin(w):
    nn = False
    for i in Ad:
        if w in Ad[i]:
            nn = True
            break
    return(nn)
    



######################################################
# quantum Bruhat graph
reflection_dict=joblib.load('DATA/'+TYPE+str(ORDER)+'/'+TYPE+str(ORDER)+'_reflection_dict')


def isedge(a,b):
    t=a.inv()+b
    if t in reflection_dict:
        if b.len()==a.len()+1:
            return(True)
        beta=reflection_dict[t]
                
        
        if b.len()==a.len()- vect.inner(beta,two_rho)+1:
            return(True)
    return(False)


   

def d_and_wt_and_ht(a,b):
    
    dist = {a:0}
    wt = {a:e0}
    queue = [a]
    i = 0
    
    while b not in queue and queue != []:
        
        q0 = queue[i]
        
        for t in reflection_dict:
            w = q0 + t
            
            if w not in queue:
                if abs(w.len()-q0.len()-1)<0.01:
                    queue += [w,]
                    dist[w] = dist[q0] + 1
                    wt[w] = wt[q0]
                if abs(w.len() - q0.len()+vect.inner(reflection_dict[t],two_rho)-1)<0.01:
                    dist[w] = dist[q0] + 1
                    queue += [w,]
                    wt[w] = wt[q0] + reflection_dict[t]
                                
                        
        i += 1
        
       
    for bb in dist:
        if bb == b:
            return(dist[bb],wt[bb].t , 1/2 * vect.inner(wt[bb],two_rho))   
            
            
def all_d_and_wt(a):
    dist = {a:0}
    wt = {a:e0}
    queue = [a]
    i = 0
    
    while len(queue)>i:
        
        q0 = queue[i]
        
        for t in reflection_dict:
            w = q0 + t
            
            if w not in queue:
                if abs(w.len()-q0.len()-1)<0.01:
                    queue += [w,]
                    dist[w] = dist[q0] + 1
                    wt[w] = wt[q0]
                if abs(w.len() - q0.len() + vect.inner(reflection_dict[t],two_rho)-1)<0.01:
                    dist[w] = dist[q0] + 1
                    queue += [w,]
                    wt[w] = wt[q0] + reflection_dict[t]
                                
                        
        i += 1
        
       
    return(dist,wt )

############ Bruhat cover, length positive, QBG

         


def delta(a):
    if vectisin(a,Phiplus):
        return(0)
    else:
        return(1)
        
def l(w,al):
    u = w.f;lam = w.t
    return(delta(al) + (lam.t * al.t).sum() - delta(act(u.inv(),al)))
def is_shrunken(w):
    for al in Phiplus:
        if abs(l(w,al))<0.001:
            return(False)
    return(True)       
def is_positive_presentation(dat,w):
    (x,mu,y) = dat
    if x + trans(mu.t) + y == w:
        nn = True
        for j in S0key:
            if delta(act(x.f,alpha[j]))+(mu.t*alpha[j].t).sum()-delta(act(y.inv().f,alpha[j]))<0:
                nn = False
                break
        if nn:
            return(True)
    return(False)




def is_P_alcove( w,J,z ):
    zwfz = z.inv()+affweyl(e0,w.f)+z
    nn = True
    for al in Phiplus:
        nal = False
        for i in S0key-J:
            if abs((omega[i].t*al.t).sum())>0.001:
                nal=True;break
        if nal:
            if delta(act(zwfz.f,al)) == 1:
                nn=False;break
            if  l(w,act(z.f,al)) < -0.001:
                nn=False;break
    return(nn)


###############################################
# Kostant p fuction and weight multiplicity
class group_algebra:
    def __init__(self,d):
        self.d = d
    def __add__(a,b):
        d = a.d
        for keyb in b.d:
            nb = False
            for key in d:
                if rational_vect.eq(key,keyb):
                    d[key]+=b.d[keyb];nb=True;break
            if nb==False:
                d[keyb]=b.d[keyb]
                
        return(group_algebra(d))
    def __mul__(a,b):
        d = dict()
        for keya in a.d:
            for keyb in b.d:
                keyab = keya + keyb
                nab=False
                for key in d:
                    if rational_vect.eq(key,keyab):
                        d[key]+=a.d[keya]*b.d[keyb]
                        nab=True;break
                if nab==False:
                    d[keyab]=a.d[keya]*b.d[keyb]
        return(group_algebra(d)) 
    def present(self):
        for nu in self.d:
            print(nu.t,self.d[nu])
            
class group_algebra_int:
    def __init__(self,d):
        self.d = d
    def __add__(a,b):
        d = a.d
        for keyb in b.d:
            if keyb in d:
                d[keyb]+=b.d[keyb]
            else:
                d[keyb]=b.d[keyb]
                
        return(group_algebra_int(d))
    def __mul__(a,b):
        d = dict()
        for keya in a.d:
            for keyb in b.d:
                keyab = keya + keyb
                if keyab in d:
                    d[keyab]+=a.d[keya]*b.d[keyb]
                else:
                    d[keyab]=a.d[keya]*b.d[keyb]
        return(group_algebra_int(d))                  
    def present(self):
        for nu in self.d:
            print(nu.t,self.d[nu])
            

# p function
'''
p = group_algebra_int({e0:1})
for al in Phiplus:
    fal = dict()
    for i in range(10):
        fal[scalar(-i,al)]=1
    p = p * group_algebra_int(fal)
joblib.dump(p,'DATA/A'+str(ORDER)+'/A'+str(ORDER)+'_Kostant_p')
'''
#Kostant_p = joblib.load('DATA/A'+str(ORDER)+'/A'+str(ORDER)+'_Kostant_p')


def MV(mu,lam):
    result=0
    for i in W0:                       
        for u in W0[i]:
            mulam = rational_vect(lam.t-act(u.f,mu).t+rho.t-act(u.f,rho).t)
            if greater_than_or_equal_to(e0, mulam):
                #print(u.rex(),mulam.t)
                for key in Kostant_p.d:
                    if rational_vect.eq(key,mulam):
                        result+=(-1)**i*Kostant_p.d[key]
                        #print((-1)**i*Kostant_p.d[key])
                        break
    return(result)

                         

    
                     
########################################
if INPUT in available_type:
    W0 = joblib.load('DATA/'+TYPE+str(ORDER)+'/'+ TYPE+str(ORDER)+'_W0')
if INPUT[0]=='A' and INPUT in available_type:
    O_length_dict = joblib.load('DATA/'+TYPE+str(ORDER)+'/'+ TYPE+str(ORDER)+'_O_length_dict')

def isminimal(w):
    if INPUT[0]=='A' and INPUT in available_type:
        
        wp = w.towp()
        if wp in O_length_dict:
            return(O_length_dict[wp] == w.len())
        
    return(one_step(w)==None)
reflection_dict=joblib.load('DATA/'+TYPE+str(ORDER)+'/'+ TYPE+str(ORDER)+'_reflection_dict')
    
