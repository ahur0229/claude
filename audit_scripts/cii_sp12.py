# Independent check of the Sp12/(Sp6 x Sp6) example: E=(L5 (x) M) + M', grading s, triple (e,d,f) in p,
# p-distinguishedness via d-weights on p^e, and ranks of e^2 on U, W.
import sympy as sp
from itertools import product
# L5 basis v0..v4; ev_k = v_{k+1}; d v_k=(2k-4)v_k; f v_k = k(5-k) v_{k-1}; Q(v_k,v_l)=(-1)^k delta_{k+l,4}
n5=5
E5=sp.zeros(5); D5=sp.zeros(5); F5=sp.zeros(5); Q5=sp.zeros(5)
for k in range(5):
    if k+1<5: E5[k+1,k]=1
    D5[k,k]=2*k-4
    if k-1>=0: F5[k-1,k]=k*(5-k)
    for l in range(5):
        if k+l==4: Q5[k,l]=(-1)**k
# check sl2 relations on L5
assert D5*E5-E5*D5==2*E5 and D5*F5-F5*D5==-2*F5 and E5*F5-F5*E5==D5
# Q5 symmetric? (-1)^k delta: Q(v_l,v_k)=(-1)^l = (-1)^(4-k)=(-1)^k -> symmetric
assert Q5==Q5.T
# invariance: Q(ev,w)+Q(v,ew)=0
assert (E5.T*Q5+Q5*E5)==sp.zeros(5)
om=sp.Matrix([[0,1],[-1,0]])
# E = L5 (x) M (dim 10) + M' (dim 2); form Omega = Q5 (x) om  +  om
Omega=sp.diag(sp.kronecker_product(Q5,om), om)
assert Omega.T==-Omega and Omega.det()!=0
I2=sp.eye(2)
e=sp.diag(sp.kronecker_product(E5,I2), sp.zeros(2))
d=sp.diag(sp.kronecker_product(D5,I2), sp.zeros(2))
f=sp.diag(sp.kronecker_product(F5,I2), sp.zeros(2))
# grading involution
S5=sp.diag(*[(-1)**k for k in range(5)])
s=sp.diag(sp.kronecker_product(S5,I2), -I2)
# checks: s symplectic, s^2=1, e in p (s e s = -e), d in h
assert s*s==sp.eye(12)
assert s.T*Omega*s==Omega
assert s*e*s==-e and s*d*s==d and s*f*s==-f
# e in sp(E): e^T Omega + Omega e = 0
for X in (e,d,f):
    assert X.T*Omega+Omega*X==sp.zeros(12)
# dims of U (s=+1), W (s=-1)
print("dim U, dim W:", sum(1 for i in range(12) if s[i,i]==1), sum(1 for i in range(12) if s[i,i]==-1))
# sp(E) basis: X with X^T Omega + Omega X = 0 -> X = Omega^{-1} S with S symmetric
Oinv=Omega.inv()
basis=[]
for i in range(12):
    for j in range(i,12):
        S=sp.zeros(12); S[i,j]=1; S[j,i]=1
        basis.append(Oinv*S)
print("dim sp12 =",len(basis))
# p = {X in sp: sXs=-X}
def vec(X): return sp.Matrix([X[i,j] for i in range(12) for j in range(12)])
# project basis onto odd part
pbasis=[]
for X in basis:
    Y=(X - s*X*s)/2
    if Y!=sp.zeros(12): pbasis.append(Y)
Pmat=sp.Matrix.hstack(*[vec(Y) for Y in pbasis])
rk=Pmat.rank()
print("dim p =",rk)
# get a basis of p via column space
cols=Pmat.columnspace()
pb=[sp.Matrix(12,12,list(c)) for c in cols]
# p^e = ker ad e on p ; and d-weights there
# Solve: [e, X]=0 for X in span(pb)
coeffs=sp.symbols('c0:%d'%len(pb))
X=sum((c*B for c,B in zip(coeffs,pb)), sp.zeros(12))
eqs=list(e*X-X*e)
sol=sp.linsolve(eqs, coeffs)
# build kernel basis by nullspace
M=sp.Matrix.hstack(*[vec(e*B-B*e) for B in pb])
ns=M.nullspace()
print("dim p^e =",len(ns))
pe=[sum((c*B for c,B in zip(list(v),pb)), sp.zeros(12)) for v in ns]
# ad d acts on p^e; compute eigenvalues
# Express ad d in basis pe
PE=sp.Matrix.hstack(*[vec(Y) for Y in pe])
A=sp.zeros(len(pe))
for j,Y in enumerate(pe):
    w=vec(d*Y-Y*d)
    solv=PE.solve_least_squares(w) if False else sp.Matrix(sp.linsolve((PE,w)).args[0]) 
    A[:,j]=solv
print("ad d eigenvalues on p^e:", A.eigenvals())
# ranks of e^2 on U and W
U=[i for i in range(12) if s[i,i]==1]; W=[i for i in range(12) if s[i,i]==-1]
e2=e*e
print("rank e^2|U =", e2.extract(range(12),U).rank(), " rank e^2|W =", e2.extract(range(12),W).rank())
# T(e)=tr(ad d | h^e)
hb=[]
for Xb in basis:
    Y=(Xb + s*Xb*s)/2
    if Y!=sp.zeros(12): hb.append(Y)
Hm=sp.Matrix.hstack(*[vec(Y) for Y in hb]); hcols=Hm.columnspace(); hbb=[sp.Matrix(12,12,list(c)) for c in hcols]
print("dim h =",len(hbb))
M2=sp.Matrix.hstack(*[vec(e*B-B*e) for B in hbb]); ns2=M2.nullspace()
he=[sum((c*B for c,B in zip(list(v),hbb)), sp.zeros(12)) for v in ns2]
HE=sp.Matrix.hstack(*[vec(Y) for Y in he]); A2=sp.zeros(len(he))
for j,Y in enumerate(he):
    w=vec(d*Y-Y*d); A2[:,j]=sp.Matrix(sp.linsolve((HE,w)).args[0])
print("dim h^e =",len(he)," T(e)=tr(ad d|h^e) =",A2.trace(), " dim p =",rk)
