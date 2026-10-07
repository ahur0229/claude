import sys; sys.path.insert(0,'.')
from e7lie import *
import sympy as sp
beta=[(2,2,3,4,3,2,1),(0,1,1,2,2,2,1),(0,0,0,0,0,0,1)]
for b1 in beta:
    for b2 in beta:
        if b1!=b2: assert ip(b1,b2)==0, "not strongly orthogonal?"
        # strongly orthogonal: sum/difference not roots
    pass
for b1 in beta:
    for b2 in beta:
        if b1!=b2:
            s=tuple(b1[i]+b2[i] for i in range(7)); t=tuple(b1[i]-b2[i] for i in range(7))
            assert s not in rootset and t not in rootset
print("beta_i strongly orthogonal: ok")
def a_elem(c):
    out={}
    for ci,bi in zip(c,beta):
        if ci!=0:
            out=add(out,vec(idx[bi]),Fr(ci)); out=add(out,vec(idx[tuple(-t for t in bi)]),Fr(ci))
    return out
def dims(c):
    s=a_elem(c)
    # centralizer: kernel of ad s on g
    M=sp.zeros(DIM,DIM)
    for j in range(DIM):
        for k,cc in br(s,vec(j)).items(): M[k,j]=sp.Rational(cc.numerator,cc.denominator)
    ns=M.nullspace()
    gs=len(ns)
    # h^s, p^s: ad s maps h->p, p->h; kernel splits
    Mh=M[:,hidx]; Mp=M[:,pidx]
    hs=len(Mh.nullspace()); ps=len(Mp.nullspace())
    return gs,hs,ps
for pat in [(2,3,5),(2,3,0),(2,2,3),(2,2,0),(2,0,0),(2,2,2)]:
    print(pat,dims(pat))
# root subsystem of ker ad s for pattern: number of roots alpha with alpha(s)=0?  s is not in a Cartan of the root basis; instead compute the rank/dim of [g^s,g^s] via derived algebra dimension
def derived_dim(c):
    s=a_elem(c)
    M=sp.zeros(DIM,DIM)
    for j in range(DIM):
        for k,cc in br(s,vec(j)).items(): M[k,j]=sp.Rational(cc.numerator,cc.denominator)
    ns=M.nullspace()
    basis=[{i:Fr(int(v[i].p),int(v[i].q)) for i in range(DIM) if v[i]!=0} for v in ns]
    cols=[]
    for i in range(len(basis)):
        for j in range(i+1,len(basis)):
            b=br(basis[i],basis[j])
            if b: cols.append(sp.Matrix([sp.Rational(b.get(k,0).numerator,b.get(k,0).denominator) if k in b else 0 for k in range(DIM)]))
    return sp.Matrix.hstack(*cols).rank()
for pat in [(2,3,5),(2,3,0),(2,2,3),(2,2,0),(2,0,0),(2,2,2)]:
    print(pat,"dim [g^s,g^s] =",derived_dim(pat))
