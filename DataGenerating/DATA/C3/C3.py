# -*- coding: utf-8 -*-
"""
Created on Tue May 21 17:34:51 2019

@author: TOPOLOGY
"""

import numpy as np;  
import joblib
import time
import math
import concurrent.futures

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
    def dualroot(self):
        return(scalar(2/(self.t * self.t).sum() , self)    )

        
def vectisin(v,SET):
    for vv in SET:
        if rational_vect.eq(v,vv):
            return(True)
    return(False)    
def vectunion(S1,S2):
    
    for t in S2:
        if vectisin(t,S1) == False:
            S1=S1|{t}
    return(S1)    
    
def vectprint(A):
    for nu in A:
        print(nu.t)

  
def scalar(a,v):
    return(rational_vect(a*v.t))

def is_dominant(lam):
    nn=True
    for al in Phiplus:
        if ((lam.t * al.t).sum())<-0.001:
            nn = False;break
    return(nn) 
def coroot_expand(v):
    expd = np.zeros(semisimplerank)
    for i in range(1,semisimplerank+1):
        expd[i-1] = (v.t * omega[i].t).sum()
    return(expd)      
def greater_than_or_equal_to(v1,v2):
    return((coroot_expand(v1) >= coroot_expand(v2)).all())
def height(v):
    h=0
    for i in S0key:
        h=h + (v.t,omega[i].t).sum()
    return(h) 


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
       
        
def power(w,r):
    i=1
    ww = w
    while i < r:
        i = i+1
        ww = ww +w
    return(ww)    


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
        v = other.rex()
        aa = self
        for i in range(0,len(v)):
            if isdecrease(aa,v[i])==False:
                aa=aa+s[v[i]]
        return(aa)
    def Ad(x,w):
        return( x + w + x.inv())      
        
    def len(self):
        l=0
        for al in Phiplus:
            if vectisin(rational_act(self.f.inv(), al ), Phiplus) :
                l += abs( (self.t.t * al.t).sum() )
            else:
                l = l + abs( (self.t.t * al.t).sum() - 1 )
                
        return(l) 
        
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
    
    
    def isstraight(self):
        return(  abs(self.len() - (self.Newton()*two_rho.t).sum())<0.01)
    
    
        
    def rex(self):
        rex = np.array([])  
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

def affweylprint(A):
    for w in A:
        print(w.rex())
        
def isdecrease(w,j): # if w + s[j] < w
    mu = w.t
    u = w.f
    if j != 0:
        uj = rational_act(u,alpha[j]).minus()
        if vectisin(uj, Phiplus):
            return(  (mu.t * uj.t).sum() < 1)
        else:
            
            return( (mu.t * uj.t).sum() < 0)
    else:
        utheta = rational_act(u,theta)
        if vectisin(utheta, Phiplus):
            return( 1 + (mu.t * utheta.t).sum() < 1)
        else:
            
            return( 1 + (mu.t * utheta.t).sum() < 0)
    
    
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
    print(new[0].rex(),new[1].t,new[2].rex())
    

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

#####################################################################   
#########################################################  Data input
############################################    
    
rank = 6
semisimplerank = 3


alpha = dict()
alpha[1] = rational_vect([1/2, -1/2, 0, 0, 1/2, -1/2])
alpha[2] = rational_vect([0, 1/2, -1/2, 1/2, -1/2, 0])
alpha[3] = rational_vect([0, 0, 1,-1, 0, 0])



def Kotteq(w,nubar):
    return(abs( (w.t.t).sum() - (nubar.t).sum())<0.001)
    #return(True)


omega=dict()
omega[3] = rational_vect([1,1,1,0,0,0])
omega[2] = rational_vect([1,1,0,0,0,0])
omega[1] = rational_vect([1,0,0,0,0,0])

omegacheck=dict()
omegacheck[3] = vect([1,1,1,0,0,0])
omegacheck[2] = vect([1,1,0,0,-1,-1])
omegacheck[1] = vect([1,0,0,0,0,-1])

e0 = vect( np.zeros(rank) )
theta = rational_vect([1,0,0,0,0,-1])
thetacheck = vect([1,0,0,0,0,-1])

###################################################



#######################################################################
Skey =set(range(0,semisimplerank+1))
S0key =Skey-{0}
Id = affweyl(e0, permutation((np.arange(1,rank+1))))
alphacheck = dict()
for i in S0key:
    alphacheck[i] = vect(alpha[i].dualroot().t)


s = joblib.load('C3_s')
S0 = set()
for i in S0key:
    S0.add(s[i])
S = {s[0]} | S0
W0 = joblib.load('C3_W0')

Phiplus = joblib.load('C3_Phiplus')
two_rho = joblib.load('C3_two_rho')
rho = joblib.load('C3_rho')
w_0 = affweyl(e0, permutation(np.arange(1,rank+1)[::-1]) )



mu3 = omegacheck[3]
y03 = exp([1,2,1]) + w_0
tau3 = trans(mu3.t)+y03
Omega0 = {tau3}
Omega = {Id,tau3}



####################################################################
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

        
########################################### compute W(n)
#W=joblib.load('C3_W(n)')
'''
t0 = time.time()
W = {0:{Id},1:S}
for i in range(1,100):
    if i not in W:
        W[i] = set()
        for x in W[i-1]:
            for j in Skey:
                if isdecrease(x,j) == False:
                    W[i].add(x + s[j])
        print(i,len(W[i]))            
print(time.time()-t0) 
     
joblib.dump(W,'C3_W_to30')    '''
 



################### defect


def nu_defect(nu):
    defectnu = 0
    for i in range(1,semisimplerank+1):
        ai =  (omega[i].t * nu.t).sum() + 0.001
        defectnu += ai - math.floor(ai)
    return(int(2*defectnu+0.001))




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


def isCoxeter(w):
    return(support(w) == S0key and w.len() == semisimplerank)
       
    
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

##################    
def tauprint(a):
    print(  (a+tau.inv()).rex())

#############################   reduction for split case, tau may not be Id

def one_step(w):
    # tau = Id
    #if isminimal(w):
    #    return(None)
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
                    return((s[j]+w, s[j] + wj))
                if dewj != dejwj:
                    jwj = s[j] + wj
                    if jwj not in O and jwj not in Onew:
                        Onewnew.add(jwj)
    return (None)     

def partial_one_step(w):
    O=set()
    Onew=set()
    Onewnew={w}
    while Onewnew != set():
        O = O|Onew
        Onew = Onewnew
        Onewnew=set()
        for w in Onew:
            for j in S0key:
                wj = w + s[j]
                dewj = isdecrease(w,j)
                dejwj = isdecrease(wj.inv(),j)
                if dewj and dejwj:
                    return((s[j]+w, s[j] + wj))
                if dewj != dejwj:
                    jwj = s[j] + wj
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
            
def partial_tominimal(w):
    while True:
        onew = partial_one_step(w)
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
            if repw is not None:
                Pnew = Pnew | {DLdatum(repw[0],p.n+1),DLdatum(repw[1],p.n+1)}
 
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
        dimdict[nu] = max(dimdict[nu])
    for dat in redw:
        nudat = rational_vect(dat.w.Newton())
        for nu in dimdict:
            if rational_vect.eq(nu,nudat):
                if (dat.w.len() +dat.n)==dimdict[nu]:
                    irrdict[nu] += 1
                    break
    dictw = dict()
    for nu in BGw:
        dictw[nu] = (dimdict[nu],irrdict[nu])
    return(dictw)

def Idim_sigma_dict(w):
    redw = sigma_reduction(w);Idimdict = dict() 
    BGw = set()
    for dat in redw:
        BGw = vectunion(BGw,{rational_vect(dat.w.sigma_Newton())})
    for nu in BGw:
        Idimdict[nu] = set() 
    for dat in redw:
        nudat = rational_vect(dat.w.sigma_Newton())
        for nu in Idimdict:
            if rational_vect.eq(nu,nudat):
                Idimdict[nu].add(dat.w.len() + dat.n)
                break
    for nu in Idimdict:
        Idimdict[nu] = max(Idimdict[nu])
    return(Idimdict)


def dim_irr_print(w):
    dictw = dim_irr_dict(w)
    for nu in dictw:
        print(nu.t,dictw[nu])







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
# W = joblib.load('C3_W(to30)')
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
    


   
### compute Adm dimdict
'''
Admmu3 = Adm(mu3,tau3)
joblib.dump(Admmu3,'C3_Admmu3')
    
BGmu3 = set()
for i in Admmu3:
    for w in Admmu3[i]:
        BGmu3 = vectunion(BGmu3,{rational_vect(w.Newton())})
joblib.dump(BGmu3,'C3_BGmu3')
vectprint(BGmu3)
dimdict3 = dict()
dimset = dict()
for nu in BGmu3:
    dimset[nu] = set()

for i in Admmu3:
    print(i)
    for w in Admmu3[i]:
        redw = reduction(w)
        for nu in BGmu3:
            for dat in redw:
                if rational_vect.eq(rational_vect(dat.w.Newton()),nu):
                    dimset[nu].add(dat.w.len() + dat.n)
for nu in dimset:
    dimdict3[nu] = max(dimset[nu])
    print(nu.t,dimdict3[nu] )
joblib.dump(dimdict3,'C3_dimdict3')

'''

######################################################

# quantum Bruhat graph
reflection_dict=joblib.load('C3_reflection_dict')


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





    
   
    
def delta(a):
    if vectisin(a,Phiplus):
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


##############################
#W = joblib.load('C3_Wa_to30')

#Admmu3 = joblib.load('C3_Admmu3')
#BGmu3 = joblib.load('C3_BGmu3')
