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
    
def vectprint(A):
    for nu in A:
        print(nu.t)

def scalar(a,v):
    return(rational_vect(a*v.t))
def root_expand(v):
    expd = np.zeros(semisimplerank)
    for i in range(1,semisimplerank+1):
        expd[i-1] = (v.t * omega[i].t).sum()
    return(expd) 
def is_dominant(lam):
    nn=True
    for al in Phiplus:
        if ((lam.t * al.t).sum())<-0.001:
            nn = False;break
    return(nn)
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
    
    
class permutation:
    def __init__(self,f):
        self.f = np.array(f, dtype = int)
    def __add__(p1,p2):
        v1 = p1.f
        v2 = p2.f
        v = np.zeros(rank)
        for i in range(rank):
            v[i] = v1[int(v2[i])-1]
        return(permutation(v)) 
        
    def inv(self):
        v1 = self.f
        v = np.zeros(rank)
        
        for j in range(rank):
            v[int(v1[j])-1] = j + 1
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
    return(vect(y))
def rational_act(w,x):
    y = np.zeros(rank)
    for i in range(0,rank):
        ai = int(w.f[i])
        y[ai-1] = x.t[i]
    return(rational_vect(y))               
       
def sigma_act(w,v):#matrix and list
    return(  act(w,  vect(  -v.t[::-1]) )   )

        
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
        return(affweyl(act(self.f.inv(), (self.t).minus()), self.f.inv()))
    
    def __eq__(self,other):
        return( self.t == other.t and self.f == other.f )     
    def __hash__(self):
        return(hash( (self.t,self.f) ) )
    def __add__(self,other):

        return(affweyl(self.t + act(self.f,other.t),self.f + other.f ))
    
    
    def __mul__(self,other):
        return(self + other)
    
    def Ad(x,w):
        return( x + w + x.inv())  
    
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
        return(int(l)) 
        
    def towp(w):
        c = w.f.tocycle()
        #weight = []
        
        wp = []
        for ci in c:
            sumi = 0
            for j in ci:
                sumi += w.t.t[j-1]
            
            
            wp.append((len(ci),sumi))
        return(tuple(sorted(wp,reverse = True)))
    
    def isdistinguish(self):
        wp = self.towp()
        nn = True
        for bc in wp:
            if gcd(int(bc[0]),int(bc[1])) != 1:
                nn = False
                break
        return(nn)        
    
    # only for split type A   
    def Newton(w):
        wp = w.towp()
        nuw = np.array([])
        for bc in wp:
            nuw = np.append(nuw, bc[1]/bc[0]*np.ones(bc[0]))
            
        return( np.array(sorted(nuw,reverse = True)) )
    
    def sigma_Newton(self):
        sigma = sigma0
        u = affweyl(e0,self.f)
        if u == Id:
            return(sorted(diamond(self.t).t,reverse = True))
        else:
            n_w = 2 * (u.f + sigma(u).f).order()
            #if n_w == 0 :
            nu0 = self.t
            nu1 = self.t
            #print(nu.t)
            for i in range(1,n_w):
                nu1 = sigma_act(u.f,nu1)
                #print(weylact(ww,self.t).t) 
                nu0 = nu0 + nu1
        return(np.array(sorted(scalar((1/n_w),nu0).t,reverse = True)) )
    
    
    
    def isstraight(self):
        return(  abs(self.len() - (self.Newton()*two_rho.t).sum())<0.01)
    
    
        
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
        print(self.t.t,affweyl(e0,self.f).rex())
        
    def p(self):
        print(self.t.t)
        print((self.f).tocycle())

def fraction_print(nu):
    fraction_array = []
    for a in nu:
        fraction_array.append(Fraction(a).limit_denominator(1000))
    fraction_string = "[" + ", ".join(str(f.numerator) if f.denominator == 1 else str(f) for f in fraction_array) + "]"
    print(fraction_string)    
    
def affweylprint(A):
    for w in A:
        print(w.rex())
        
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



def affweylpower(w,r):
    i = 1
    ww = w
    while i < r:
        i = i + 1
        ww = ww + w
    return(ww)    

def tau_affweylpower(w,r):
    i=1
    ww = w
    while i<r:
        i = i + 1
        ww = w + tau + ww + tau.inv()  
    return(ww)     
    
                     
def exp(exp):
    w = Id
    for i in range(0,len(exp)):
        u = exp[i]
        w = w + s[u]
        
    return(w)
    

def trans(v):
    return(affweyl(vect(v),Id.f))
    
def star(a,b):
    
    rexb = b.rex()
    for i in range(0,len(rexb)):
        if isdecrease(a,rexb[i])==False:
            a = a + s[rexb[i]]
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
    if len(I) == semisimplerank:
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
            if (sorted(xwt.t,reverse = True) == xwt.t).all():        
                de = ldecom( x + wf, cent(xwt))
                return((x.inv()+de[0]), xwt, de[1]  )
            
def partialprint(w):
    new = partial(w)          
    print(new[0].rex())
    print(new[1].t )
    print(new[2].rex())    
    

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
#########################################################  Data input
############################################  
available_type = {'A2','A3','A4','A5',
                  'B2','B3','B4','B5',
                  'C2','C3','C4','C5', 
                  'D4','D5','D6',
                  '2A2','2A3','2A4','2A5',
                  '2D4','2D5','2D6'}

INPUT = input('please input type of the group (An, 2An, Bn, Cn, Dn, 2Dn): ') 
if len(INPUT) == 2:
    TYPE = INPUT[0]; ORDER = int(INPUT[1])
else:
    SIGMA_TYPE = INPUT[0]; TYPE = INPUT[1]; ORDER = int(INPUT[2])
if TYPE == 'A':
    rank = ORDER + 1
    semisimplerank = ORDER
    if INPUT in available_type:
        alpha = joblib.load('A'+str(ORDER)+'_alpha')
        alphacheck = joblib.load('A'+str(ORDER)+'_alphacheck')
        omega = joblib.load('A'+str(ORDER)+'_omega')
        omegacheck = joblib.load('A'+str(ORDER)+'_omegacheck')
        #omega12 = rational_vect([1,0,-1])
        #Omega0_sigma = {omega12}
        Skey =set(range(0,semisimplerank+1))
        S0key =Skey-{0}
        e0 = vect((np.zeros(rank)))
        Id = affweyl(e0, permutation((np.arange(1,rank+1))))
        s = joblib.load('A'+str(ORDER)+'_s')
        S0 = set()
        for i in S0key:
            S0 = S0|{s[i]}
        W0 = joblib.load('A'+str(ORDER)+'_W0')
        Phiplus = joblib.load('A'+str(ORDER)+'_Phiplus')
        S={s[0]}|S0
        two_rho = joblib.load('A'+str(ORDER)+'_two_rho')
        rho = joblib.load('A'+str(ORDER)+'_rho')
        w_0 = joblib.load('A'+str(ORDER)+'_w_0')
        theta = joblib.load('A'+str(ORDER)+'_theta')
        thetacheck = joblib.load('A'+str(ORDER)+'_thetacheck')
        tau = joblib.load('A'+str(ORDER)+'_tau')
        Omega = set(tau.values())

        #Omega0_sigma = {omega12}
        O_length_dict = dict()
        O_length_dict[0] = joblib.load('A'+str(ORDER)+'_O_length_dict')
        for j in range(1,ORDER + 1):
            O_length_dict[j] = joblib.load('A'+str(ORDER)+'_O_length_dict_tau'+str(j))
        
    else:
        e0 = vect((np.zeros(rank)))
        Id = affweyl(e0, permutation((np.arange(1,rank+1))))
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
        Skey =set(range(0,semisimplerank+1))
        S0key =Skey-{0}
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
            tau[j] = trans(omegacheck[j].t)+affweyl(e0,permutation(fj))+w_0    
        Omega = set(tau.values())

if TYPE == 'B':
    rank = 2*ORDER + 1
    semisimplerank = ORDER
    if INPUT in available_type:
        alpha = joblib.load('B'+str(ORDER)+'_alpha')
        alphacheck = joblib.load('B'+str(ORDER)+'_alphacheck')
        omega = joblib.load('B'+str(ORDER)+'_omega')
        omegacheck = joblib.load('B'+str(ORDER)+'_omegacheck')
        #omega12 = rational_vect([1,0,-1])
        #Omega0_sigma = {omega12}
        Skey =set(range(0,semisimplerank+1))
        S0key =Skey-{0}
        e0 = vect((np.zeros(rank)))
        Id = affweyl(e0, permutation((np.arange(1,rank+1))))
        s = joblib.load('B'+str(ORDER)+'_s')
        S0 = set()
        for i in S0key:
            S0 = S0|{s[i]}
        W0 = joblib.load('B'+str(ORDER)+'_W0')
        Phiplus = joblib.load('B'+str(ORDER)+'_Phiplus')
        S={s[0]}|S0
        two_rho = joblib.load('B'+str(ORDER)+'_two_rho')
        rho = joblib.load('B'+str(ORDER)+'_rho')
        w_0 = joblib.load('B'+str(ORDER)+'_w_0')
        theta = joblib.load('B'+str(ORDER)+'_theta')
        thetacheck = joblib.load('B'+str(ORDER)+'_thetacheck')
        tau = joblib.load('B'+str(ORDER)+'_tau')
        Omega = set(tau.values())
        
    else:
        e0 = vect((np.zeros(rank)))
        Id = affweyl(e0, permutation((np.arange(1,rank+1))))
        alpha=dict();alphacheck=dict();omega=dict();omegacheck=dict();
        s = dict()
        for j in range(1,semisimplerank):
            v=np.zeros(rank);v[j-1]=1;v[j]=-1;v[rank-j-1]=1;v[rank-j]=-1
            alpha[j] = rational_vect(1/2*v)
            alphacheck[j] = vect(v)
        v=np.zeros(rank);v[semisimplerank-1]=2;v[semisimplerank+1]=-2
        alphacheck[semisimplerank]=vect(v);alpha[semisimplerank]=rational_vect(1/4*v)
        
        for j in range(1,semisimplerank+1):
            v = np.zeros(rank)
            for k in range(1,j+1):
                v[k-1]=1;v[rank-k]=-1
            omega[j] = rational_vect(1/2*v)
            omegacheck[j] = vect(v)
        omega[semisimplerank]=scalar(1/2, omega[semisimplerank])
        
        
        Skey =set(range(0,semisimplerank+1))
        S0key =Skey-{0}
        for j in range(1,semisimplerank+1):
            f = np.array(range(1,rank+1))
            if j < semisimplerank:
                f[j-1]=j+1;f[j]=j;f[rank-j-1]=rank-j+1;f[rank-j]=rank-j
            if j == semisimplerank:
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
        for i in range(1,semisimplerank):
            for j in range(i+1, semisimplerank+1):
                v=np.zeros(rank);
                v[i-1]=1/2;v[j-1]=-1/2
                v[rank-j]=1/2;v[rank-i]=-1/2
                Phiplus=vectunion(Phiplus,{rational_vect(v)})
                v=np.zeros(rank);
                v[i-1]=1/2;v[j-1]=1/2
                v[rank-j]=-1/2;v[rank-i]=-1/2
                Phiplus=vectunion(Phiplus,{rational_vect(v)})
        for k in range(1,semisimplerank+1):
            v=np.zeros(rank);v[k-1]=1/2;v[rank-k]=-1/2
            Phiplus=vectunion(Phiplus,{rational_vect(v)})
        
        
        two_rho = rational_vect(e0.t)            
        for al in Phiplus:
            two_rho = two_rho + al
        rho = scalar(1/2,two_rho)
        w_0 = affweyl(e0,permutation(np.array(range(1,rank+1)[::-1])))
        
        tau = dict(); tau[0]=Id
        f=np.array(range(1,rank+1));f[0]=rank;f[-1]=1
        tau[1] = trans(omegacheck[1].t)+affweyl(e0,permutation(f))    
        Omega = set(tau.values())   

if TYPE == 'C':
    rank = 2*ORDER
    semisimplerank = ORDER
    if INPUT in available_type:
        alpha = joblib.load('C'+str(ORDER)+'_alpha')
        alphacheck = joblib.load('C'+str(ORDER)+'_alphacheck')
        omega = joblib.load('C'+str(ORDER)+'_omega')
        omegacheck = joblib.load('C'+str(ORDER)+'_omegacheck')
        #omega12 = rational_vect([1,0,-1])
        #Omega0_sigma = {omega12}
        Skey =set(range(0,semisimplerank+1))
        S0key =Skey-{0}
        e0 = vect((np.zeros(rank)))
        Id = affweyl(e0, permutation((np.arange(1,rank+1))))
        s = joblib.load('C'+str(ORDER)+'_s')
        S0 = set()
        for i in S0key:
            S0 = S0|{s[i]}
        W0 = joblib.load('C'+str(ORDER)+'_W0')
        Phiplus = joblib.load('C'+str(ORDER)+'_Phiplus')
        S={s[0]}|S0
        two_rho = joblib.load('C'+str(ORDER)+'_two_rho')
        rho = joblib.load('C'+str(ORDER)+'_rho')
        w_0 = joblib.load('C'+str(ORDER)+'_w_0')
        theta = joblib.load('C'+str(ORDER)+'_theta')
        thetacheck = joblib.load('C'+str(ORDER)+'_thetacheck')
        tau = joblib.load('C'+str(ORDER)+'_tau')
        Omega = set(tau.values())
        
    else:
        e0 = vect((np.zeros(rank)))
        Id = affweyl(e0, permutation((np.arange(1,rank+1))))
        alpha=dict();alphacheck=dict();omega=dict();omegacheck=dict();
        s = dict()
        for j in range(1,semisimplerank):
            v=np.zeros(rank);v[j-1]=1;v[j]=-1;v[rank-j-1]=1;v[rank-j]=-1
            alpha[j] = rational_vect(1/2*v)
            alphacheck[j] = vect(v)
        v=np.zeros(rank);v[semisimplerank-1]=1;v[semisimplerank]=-1
        alphacheck[semisimplerank]=vect(v);alpha[semisimplerank]=rational_vect(v)
        
        for j in range(1,semisimplerank+1):
            v = np.zeros(rank)
            for k in range(1,j+1):
                v[k-1]=1
            omega[j] = rational_vect(v)
            omegacheck[j] = vect(2*v)
        omegacheck[semisimplerank]=vect(omega[semisimplerank].t)
        
        
        Skey =set(range(0,semisimplerank+1))
        S0key =Skey-{0}
        for j in range(1,semisimplerank+1):
            f = np.array(range(1,rank+1))
            if j < semisimplerank:
                f[j-1]=j+1;f[j]=j;f[rank-j-1]=rank-j+1;f[rank-j]=rank-j
            if j == semisimplerank:
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
        for i in range(1,semisimplerank):
            for j in range(i+1, semisimplerank+1):
                v=np.zeros(rank);
                v[i-1]=1/2;v[j-1]=-1/2
                v[rank-j]=1/2;v[rank-i]=-1/2
                Phiplus=vectunion(Phiplus,{rational_vect(v)})
                v=np.zeros(rank);
                v[i-1]=1/2;v[j-1]=1/2
                v[rank-j]=-1/2;v[rank-i]=-1/2
                Phiplus=vectunion(Phiplus,{rational_vect(v)})
        for k in range(1,semisimplerank+1):
            v=np.zeros(rank);v[k-1]=1;v[rank-k]=-1
            Phiplus=vectunion(Phiplus,{rational_vect(v)})
        
        
        two_rho = rational_vect(e0.t)            
        for al in Phiplus:
            two_rho = two_rho + al
        rho = scalar(1/2,two_rho)
        w_0 = affweyl(e0,permutation(np.array(range(1,rank+1)[::-1])))
        
        tau = dict(); tau[0]=Id
        f1=np.array(range(1,semisimplerank+1))
        f2=np.array(range(semisimplerank+1,rank+1))
        f = np.append(f2,f1)
        tau[semisimplerank] = trans(omegacheck[semisimplerank].t)+affweyl(e0,permutation(f))    
        Omega = set(tau.values())   

if TYPE == 'D':
    rank = 2*ORDER
    semisimplerank = ORDER
    if INPUT in available_type:
        alpha = joblib.load('D'+str(ORDER)+'_alpha')
        alphacheck = joblib.load('D'+str(ORDER)+'_alphacheck')
        omega = joblib.load('D'+str(ORDER)+'_omega')
        omegacheck = joblib.load('D'+str(ORDER)+'_omegacheck')
        #omega12 = rational_vect([1,0,-1])
        #Omega0_sigma = {omega12}
        Skey =set(range(0,semisimplerank+1))
        S0key =Skey-{0}
        e0 = vect((np.zeros(rank)))
        Id = affweyl(e0, permutation((np.arange(1,rank+1))))
        s = joblib.load('D'+str(ORDER)+'_s')
        S0 = set()
        for i in S0key:
            S0 = S0|{s[i]}
        W0 = joblib.load('D'+str(ORDER)+'_W0')
        Phiplus = joblib.load('D'+str(ORDER)+'_Phiplus')
        S={s[0]}|S0
        two_rho = joblib.load('D'+str(ORDER)+'_two_rho')
        rho = joblib.load('D'+str(ORDER)+'_rho')
        w_0 = joblib.load('D'+str(ORDER)+'_w_0')
        theta = joblib.load('D'+str(ORDER)+'_theta')
        thetacheck = joblib.load('D'+str(ORDER)+'_thetacheck')
        tau = joblib.load('D'+str(ORDER)+'_tau')
        Omega = set(tau.values())
        
    else:
        e0 = vect((np.zeros(rank)))
        Id = affweyl(e0, permutation((np.arange(1,rank+1))))
        alpha=dict();alphacheck=dict();omega=dict();omegacheck=dict();
        s = dict()
        for j in range(1,semisimplerank):
            v=np.zeros(rank);v[j-1]=1;v[j]=-1;v[rank-j-1]=1;v[rank-j]=-1
            alpha[j] = rational_vect(1/4*v)
            alphacheck[j] = vect(2*v)
        v=np.zeros(rank);v[semisimplerank-2]=1;v[semisimplerank-1]=1
        v[semisimplerank]=-1;v[semisimplerank+1]=-1
        alphacheck[semisimplerank]=vect(2*v);alpha[semisimplerank]=rational_vect(1/4*v)
        
        for j in range(1,semisimplerank-1):
            v = np.zeros(rank)
            for k in range(1,j+1):
                v[k-1]=1;v[rank-k]=-1
            omega[j] = rational_vect(1/4*v)
            omegacheck[j] = vect(2*v)
        v = np.zeros(rank)
        for k in range(1,semisimplerank+1):
            v[k-1]=1;v[rank-k]=-1
        omega[semisimplerank]=rational_vect(1/8*v)
        omegacheck[semisimplerank]=vect(v)
        v[semisimplerank-1]=-1;v[semisimplerank]=1
        omega[semisimplerank-1]=rational_vect(1/8*v)
        omegacheck[semisimplerank-1]=vect(v)
        
        
        Skey =set(range(0,semisimplerank+1))
        S0key =Skey-{0}
        for j in range(1,semisimplerank+1):
            f = np.array(range(1,rank+1))
            if j < semisimplerank:
                f[j-1]=j+1;f[j]=j;f[rank-j-1]=rank-j+1;f[rank-j]=rank-j
            if j == semisimplerank:
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
        for i in range(1,semisimplerank):
            for j in range(i+1, semisimplerank+1):
                v=np.zeros(rank);
                v[i-1]=1/4;v[j-1]=-1/4
                v[rank-j]=1/4;v[rank-i]=-1/4
                Phiplus=vectunion(Phiplus,{rational_vect(v)})
                v=np.zeros(rank);
                v[i-1]=1/4;v[j-1]=1/4
                v[rank-j]=-1/4;v[rank-i]=-1/4
                Phiplus=vectunion(Phiplus,{rational_vect(v)})
        
        two_rho = rational_vect(e0.t)            
        for al in Phiplus:
            two_rho = two_rho + al
        rho = scalar(1/2,two_rho)
        w_0 = affweyl(e0,permutation(np.array(range(1,rank+1)[::-1])))
        
        tau = dict(); tau[0]=Id
        f=np.array(range(1,rank+1))
        f[0]=rank;f[-1]=1
        f[semisimplerank-1]=semisimplerank+1
        f[semisimplerank]=semisimplerank
        tau[1]=trans(omegacheck[1].t)+affweyl(e0,permutation(f))    
        f1=np.array(range(1,semisimplerank+1))
        f2=np.array(range(semisimplerank+1,rank+1))
        f = np.append(f2,f1);f[0]=semisimplerank;f[-1]=semisimplerank+1
        tau[semisimplerank-1] = trans(omegacheck[semisimplerank-1].t)+affweyl(e0,permutation(f)) 
        f = np.append(f2,f1);f[semisimplerank-1]=1;f[semisimplerank]=rank
        tau[semisimplerank] = trans(omegacheck[semisimplerank].t)+affweyl(e0,permutation(f)) 
        Omega = set(tau.values())   

            
'''
def idty(x):
    return(x)
def sigma0(x):
    rexx=x.rex()
    sigx=Id
    for i in range(0,len(rexx)):
        if rexx[i]==0:
            sigx=sigx+s[0]
        if rexx[i]==1:
            sigx=sigx+s[2]
        if rexx[i]==2:
            sigx=sigx+s[1]  
    return(sigx)
def sigma00(x):
    return(x)    
def sigma10(x):
    return(tau1+x+tau1.inv())    
def sigma20(x):
    return(tau2+x+tau2.inv())    
def sigma01(x):
    return(  sigma00(sigma0(x)) )
def sigma11(x):
    return(  sigma10(sigma0(x)) )
def sigma21(x):
    return(  sigma20(sigma0(x)))
def diamond(mu):
    return(scalar(1/2,rational_vect(mu.t-mu.t[::-1])))  
'''

        

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
    for i in range(1,semisimplerank+1):
        ai =  (omega[i].t * nu.t).sum() + 0.001
        defectnu += ai - math.floor(ai)
    return(int(2*defectnu+0.001))
def length(nu,mu):
    l = 0; mu_nu = rational_vect(mu.t-nu.t)
    for i in range(1,semisimplerank+1):
        ai =  vect.inner(omega[i],mu_nu)-0.001
        l +=  math.ceil(ai)
    return(int(l+0.001))
def sigma_length(nu,mu):
    l = 0; mu_nu = mu.t-nu.t
    for omega_c in omega_sigma_set:
        l+=math.ceil( (mu_nu*omega_c.t).sum()-0.001)
    return(l)      


def defect(b):
    nu = rational_vect(b.Newton())
    return(nu_defect(nu))
         
#def defect(b):
#    return(np.linalg.matrix_rank(permutation_to_matrix(Id.f) - permutation_to_matrix(b.f) ))
 

def eta_delta(w):
    parw = partial(w)
    return( parw[2] + parw[0])

def maxa(w):
    aw = set()
    for a in numeration(w):
        aw.add(1/2*(a.len()+eta_delta(a).len()))
    aw = max(aw)
    return(aw)
def mina(w):
    aw = set()
    for a in numeration(w):
        aw.add(1/2*(a.len()+eta_delta(a).len()))
    aw = min(aw)
    return(aw)


def virtualdim(w,nu): 
    return(1/2*(w.len()+eta_delta(w).len()-nu_defect(nu)- (two_rho.t*nu.t).sum()))


def is_Coxeter(w):
    return(support(w) == S0key and w.len() == semisimplerank)
       
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
    if isminimal(w):
        return(None)
    O=set()
    Onew=set()
    Onewnew={w}
    while Onewnew != set():
        O = O|Onew
        Onew = Onewnew
        Onewnew=set()
        for w in Onew:
            for j in Skey:
                wj = w + s[j]
                dewj = isdecrease(w,j)
                dejwj = isdecrease(wj.inv(),j)
                if dewj and dejwj:
                    #partialprint(w)
                    #print(j)
                    return((s[j]+w, s[j] + wj))
                if dewj != dejwj:
                    jwj = s[j] + wj
                    if jwj not in O and jwj not in Onew:
                        Onewnew.add(jwj)
    return (None)     

def sigma_one_step(w):
    sigma = sigma0
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
            for j in Skey:
                wj = w + s[j]
                if isdecrease(w,j) != isdecrease(wj.inv(),j):
                    jwj = s[j] + wj
                    if jwj not in O and jwj not in Onew:
                        Onewnew.add(jwj)
    # Omega conjugation    
    O = O|Onew 
    OO = set()
    for w in O:
        OO.add(w)
    for ta in Omega0:
        for w in O:
            OO.add(ta+w+ta.inv())
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
    

def sigma_reduction(w):
    P = set();sigma = sigma0
    Pnew = {DLdatum(w,0)}
    while len(Pnew) != len(P):
        P = Pnew
        Pnew = set()        
        for p in P :
            repw = sigma_one_step(p.w)
            if repw is None:
                Pnew.add(p)
            else:
                Pnew = Pnew | {DLdatum(repw[0],p.n+1),DLdatum(repw[1],p.n+1)}
    return(  Pnew    )
    
def reduction_step(w,d):
    P = set(); r=0; Pnew = {DLdatum(w,0)}
    while r<d and len(Pnew) != len(P):
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
        print(f'dat.w.Newton() = {dat.w.Newton()}' )
        print(f'dat.n + dat.w.len() = {dat.n + dat.w.len()}')
        print('\n')  


def dim(w,v = e0):
    dimset=set()
    for dat in reduction(w):
        if rational_vect.eq(rational_vect(dat.w.Newton()),v):
            dimset.add(dat.w.len() - (v.t*two_rho.t).sum() + dat.n)
            
    if dimset==set():
        return('empty')
    return(math.floor(max(dimset)+0.001))  
def Idim(w,v = e0):
    dimset=set()
    for dat in reduction(w):
        if rational_vect.eq(rational_vect(dat.w.Newton()),v):
            dimset.add(dat.w.len() + dat.n)        
    if dimset==set():
        return('empty')
    return(math.floor(max(dimset)+0.001))  
def irr(w, v = e0):
    redw = reduction(w);dimw = set()
    for dat in redw:
        if rational_vect.eq(rational_vect(dat.w.Newton()),v):
            dimw.add(dat.w.len()-(v.t*two_rho.t).sum()+dat.n)
    if dimw==set():
        return(0)
    else:
        dimw = math.floor(max(dimw)+0.001); lc = 0
        for dat in redw:
            if rational_vect.eq(rational_vect(dat.w.Newton()),v):
                if math.floor(dat.w.len()-(v.t*two_rho.t).sum()+dat.n+0.001)==dimw:
                    lc+=1
        return(lc)  
def dim_irr_dict(w):
    redw = reduction(w);dimdict = dict();irrdict = dict()
    BGw = set()
    for dat in redw:
        BGw = vectunion(BGw,{rational_vect(dat.w.Newton())})
    for nu in BGw:
        dimdict[nu] = set(); irrdict[nu] = 0
    for dat in redw:
        nudat = rational_vect(dat.w.Newton())
        for nu in dimdict:
            nutworho = (nu.t*two_rho.t).sum()
            if rational_vect.eq(nu,nudat):
                dimdict[nu].add(dat.w.len() + dat.n - nutworho)
                break
    for nu in dimdict:
        dimdict[nu] = math.floor(0.001+max(dimdict[nu]))
    for dat in redw:
        nudat = rational_vect(dat.w.Newton())
        for nu in dimdict:
            nutworho = (nu.t*two_rho.t).sum()
            if rational_vect.eq(nu,nudat):
                if abs((dat.w.len() -nutworho +dat.n)-dimdict[nu])<0.001:
                    irrdict[nu] += 1
                    break
    dictw = dict()
    for nu in BGw:
        dictw[nu] = (dimdict[nu],irrdict[nu])
    return(dictw)

def dim_irr_sigma_dict(w):
    redw = reduction(w);dimdict = dict();irrdict = dict()
    BGw = set()
    for dat in redw:
        BGw = vectunion(BGw,{rational_vect(dat.w.sigma_Newton())})
    for nu in BGw:
        dimdict[nu] = set(); irrdict[nu] = 0
    for dat in redw:
        nudat = rational_vect(dat.w.sigma_Newton())
        for nu in dimdict:
            nutworho = (nu.t*two_rho.t).sum()
            if rational_vect.eq(nu,nudat):
                dimdict[nu].add(dat.w.len() + dat.n - nutworho)
                break
    for nu in dimdict:
        dimdict[nu] = int(max(dimdict[nu]))
    for dat in redw:
        nudat = rational_vect(dat.w.sigma_Newton())
        for nu in dimdict:
            nutworho = (nu.t*two_rho.t).sum()
            if rational_vect.eq(nu,nudat):
                if abs((dat.w.len() -nutworho +dat.n)-dimdict[nu])<0.001:
                    irrdict[nu] += 1
                    break
    dictw = dict()
    for nu in BGw:
        dictw[nu] = (dimdict[nu],irrdict[nu])
    return(dictw)


def dim_irr_print(w):
    dictw = dim_irr_dict(w)
    for nu in dictw:
        print(f'Newton point = {nu.t}', f'dim = {dictw[nu][0]}', f'irr = {dictw[nu][1]}' )

def dim_irr_sigma_print(w):
    dictw = dim_irr_sigma_dict(w)
    for nu in dictw:
        print(f'Newton point = {nu.t}', f'dim = {dictw[nu][0]}', f'irr = {dictw[nu][1]}' )







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
        BGw = vectunion(BGw,{rational_vect(dat.w.Newton())})
    return (BGw)

def BG_sigma(w):
    BGw = set()
    redw = sigma_reduction(w)
    for dat in redw:
        BGw = vectunion(BGw,{rational_vect(dat.w.sigma_Newton())})
    return (BGw)


def iscordial(w):
    nu = rational_vect(max_Newton(w))
    return(  abs(w.len()- eta_delta(w).len() - (nu.t*two_rho.t).sum() + nu_defect(nu))<0.01)

def issemicordial(w):
    numax = rational_vect(max_Newton(w))
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
    D0 = virtualdim(w, nu0)-dim(w,nu0)
    if abs(D-D0)<0.01:
        return(True)
    else:
        #print(nu0.t)
        return(False)
    
    
def triangle_operator(a,b):
    
    rexb=b.rex()
    for i in range(0,len(rexb)):
        if isdecrease(a,rexb[i]): 
            a=a+s[rexb[i]]
    return(a)  
    
###########################################
# W = joblib.load('A2_W(to30)')
####################### compute X(mu,b)

def Adm(mu,tau): #  w + tau <= trans(mu), w in W
    N = int(trans(mu.t).len())
    Admmu = dict()
    for i in range(N+1):
        Admmu[i] = set()
    
    for i in W0:
        for u in W0[i]:
            Admmu[N].add(u + trans(mu.t) + u.inv())
    
    for i in range(1,N+1 ):
        j = N - i
        print(j)
        for w in Admmu[j+1]:
            wrex = (w+tau.inv()).rex()
            for k in range(j+1):
                wrexk = np.delete(wrex, k)
                wk = exp(wrexk)+tau;lwk=wk.len()
                Admmu[lwk].add(wk)
    return(Admmu)
    
    
def Admisin(w):
    nn = False
    for i in Ad:
        if w in Ad[i]:
            nn = True
            break
    return(nn)
    



######################################################
# quantum Bruhat graph
reflection_dict=joblib.load('A2_reflection_dict')


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
'''
w1 = trans([-1,-1,1,1])+exp([2,1,3,2])
parw1 = partial(w1); x1 = parw1[0]; mu1 = parw1[1]; y1 = parw1[2]
print(f'x1, mu1,y1 = {x1.rex(),mu1.t,y1.rex()}')
print(f'L1 = {coroot_expand(mu1)}')
print('\n')
for t in reflection_dict:
    beta = reflection_dict[t]
    for i in range(-10,10):
        it = trans( scalar(i,beta).t)+t
        if (it+w1).len() == w1.len()+1:
            w2  = it + w1
            partialprint(w2)
            u2 = affweyl(e0,w2.f);mu = act(u2.inv().f,w2.t)
            
            for i in W0:
                for y2 in W0[i]:
                    x2 = u2 + y2.inv() 
                    mu2 = act(y2.f,mu)
                    nn = True
                    for j in S0key:
                        if (mu2.t*alpha[j].t).sum()+sign(act(y2.inv().f,alpha[j]))-sign(act(x2.f,alpha[j])) < 0:
                            nn = False
                            break
                    if nn:
                        print(f'x2, mu2,y2 = {x2.rex(),mu2.t,y2.rex()}')
                        print(f'L = { coroot_expand(vect(mu2.t - d_and_wt_and_ht(y2.inv(),y1.inv())[1] - d_and_wt_and_ht(x1,x2)[1]))}')
            print('\n')
# 固定 w1, x1, mu1, y1, 找到所有w2, x2, mu2, y2 并计算L(w2)
            
'''            


def delta(a):
    if a in Phiplus:
        return(0)
    else:
        return(1)
def l(w,al):
    u = w.f;lam = w.t
    return(delta(al) + (lam.t * al.t).sum() - delta(act(u.inv(),al)))
def isshrunken(w):
    for al in Phiplus:
        if l(w,al) == 0:
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
    (x,mu,y)=(z,act(z.inv().f,w.t),z.inv()+affweyl(e0,w.f))
    yx = y+x
    nn = True
    for al in Phiplus:
        nal = False
        for i in S0key-J:
            if (omega[i].t*al.t).sum()!=0:
                nal=True;break
        if nal:
            if delta(act(yx.f,al)) == 1:
                nn=False;break
            if  l(w,act(z.f,al)) < 0:
                nn=False;break
    return(nn)
    
    
    
########################################

#W = joblib.load('A2_Wa_to30')
#Adm_mu1 = joblib.load('A2_Admmu1')
#BG_mu1 = joblib.load('A2_BGmu1')