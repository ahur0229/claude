import sys; sys.path.insert(0,'.')
from e7lie import *
import sympy as sp
def parse(d): return {k:Fr(v) for k,v in d.items()}
e15=parse({56:-1,59:-1,62:1,64:-1,65:1,66:1,67:1,68:1,77:-1,83:-1,89:1,95:-1,96:1,101:1,106:-1,111:1})
d15=parse({0:4,1:4,2:6,3:8,4:6,5:4})
cvec=[d15.get(i,0) for i in range(7)]
print("characteristic:",[sum(C[i][j]*cvec[j] for j in range(7)) for i in range(7)])
print("[d,e]=2e:",add(br(d15,e15),e15,-2)=={})
# find f: solve [e,f]=d, [d,f]=-2f with f in p
pw={}
for i in pidx:
    r=roots[i-7]; w=sum(C[k][j]*cvec[j] for k in range(7) for j in range(7) if r[k]!=0 and False)
# weight of E_alpha under d: alpha(d) = sum_i a_i alpha_i(d)
chars=[sum(C[i][j]*cvec[j] for j in range(7)) for i in range(7)]
def wt(i): return sum(roots[i-7][k]*chars[k] for k in range(7))
pm2=[i for i in pidx if wt(i)==-2]
syms=sp.symbols('f0:%d'%len(pm2))
fdict={}
for s,i in zip(syms,pm2): fdict[i]=s
# [e,f] = d: compute symbolically
out={}
for i,ci in e15.items():
    for j,s in fdict.items():
        for k,c in BR[(i,j)].items():
            out[k]=out.get(k,0)+sp.Rational(ci.numerator,ci.denominator)*s*sp.Rational(c.numerator,c.denominator)
eqs=[]
for k in range(DIM):
    target=sp.Rational(d15.get(k,0).numerator,d15.get(k,0).denominator) if k in d15 else 0
    eqs.append(sp.expand(out.get(k,0)-target))
sol=sp.solve(eqs,syms,dict=True)
print("solutions for f:",len(sol))
f15={i:Fr(int(sp.nsimplify(sol[0][s]).p),int(sp.nsimplify(sol[0][s]).q)) for s,i in zip(syms,pm2) if sol[0].get(s,0)!=0}
print("[e,f]-d:",add(br(e15,f15),d15,-1)=={}, "[d,f]+2f:",add(br(d15,f15),f15,2)=={})
# p^e weight zero kernel
p0=[i for i in pidx if wt(i)==0]
M=sp.zeros(DIM,len(p0))
for b,j in enumerate(p0):
    for k,c in br(e15,vec(j)).items(): M[k,b]=sp.Rational(c.numerator,c.denominator)
ns=M.nullspace()
print("dim weight-zero kernel of ad e in p:",len(ns))
zs=[{p0[i]:Fr(int(v[i].p),int(v[i].q)) for i in range(len(p0)) if v[i]!=0} for v in ns]
G=sp.Matrix([[sp.Rational(B(zi,zj).numerator,B(zi,zj).denominator) for zj in zs] for zi in zs])
print("Gram matrix of B on weight-zero kernel:",G,"det",G.det())
z1=parse({40:2,53:1,57:1,60:-1,63:1,115:1}); z2=parse({52:1,103:Fr(1,2),108:Fr(-1,2),112:Fr(-1,2),117:Fr(1,2),126:1})
print("manuscript z1,z2 in kernel:",br(e15,z1)=={},br(e15,z2)=={}," B(z1,z1),B(z2,z2),B(z1,z2):",B(z1,z1),B(z2,z2),B(z1,z2))
# T(e) and Kazhdan weights for this triple
# check dims of p^e and h^e
Me=sp.zeros(DIM,DIM)
for j in range(DIM):
    for k,c in br(e15,vec(j)).items(): Me[k,j]=sp.Rational(c.numerator,c.denominator)
print("dim p^e:",len(Me[:,pidx].nullspace())," dim h^e:",len(Me[:,hidx].nullspace()))
