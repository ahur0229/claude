# Independent E7 Chevalley basis following the manuscript's stated conventions (Section 5.x "Exact coordinates").
from fractions import Fraction as Fr
import itertools, sys
edges=[(1,3),(3,4),(4,5),(5,6),(6,7),(2,4)]
n=7
C=[[0]*n for _ in range(n)]
for i in range(n): C[i][i]=2
for a,b in edges: C[a-1][b-1]=-1; C[b-1][a-1]=-1
simple=[tuple(1 if j==i else 0 for j in range(n)) for i in range(n)]
# generate positive roots by adjoining simple roots when inner product is -1
def ip(a,b): return sum(a[i]*C[i][j]*b[j] for i in range(n) for j in range(n))
pos=set(simple); frontier=list(simple)
while frontier:
    new=[]
    for r in frontier:
        for i in range(n):
            if ip(r,simple[i])==-1:
                rr=tuple(r[j]+simple[i][j] for j in range(n))
                if rr not in pos: pos.add(rr); new.append(rr)
    frontier=new
pos=sorted(pos,key=lambda r:(sum(r),r))   # height, then ascending lexicographic
assert len(pos)==63
roots=pos+[tuple(-x for x in r) for r in pos]
idx={}   # basis index: 0..6 Cartan, 7..69 positive, 70..132 negative
for k,r in enumerate(roots): idx[r]=7+k
rootset=set(roots)
DIM=133
def eps(a,b):
    s=sum(a[i]*b[i] for i in range(n))
    s+=sum(a[i]*b[j] for i in range(n) for j in range(i+1,n) if C[i][j]==-1)
    return -1 if s%2 else 1
# bracket table: brk[(i,j)] = dict index->coeff
def br_basis(i,j):
    """bracket of basis vectors i,j as dict"""
    if i<7 and j<7: return {}
    if i<7:
        r=roots[j-7]; c=sum(C[i][k]*r[k] for k in range(n))
        return {j:Fr(c)} if c else {}
    if j<7:
        d=br_basis(j,i); return {k:-v for k,v in d.items()}
    a=roots[i-7]; b=roots[j-7]
    s=tuple(a[k]+b[k] for k in range(n))
    if all(x==0 for x in s):
        return {k:Fr(-a[k]) for k in range(n) if a[k]!=0}
    if s in rootset:
        return {idx[s]:Fr(eps(a,b))}
    return {}
BR={}
for i in range(DIM):
    for j in range(DIM):
        BR[(i,j)]=br_basis(i,j)
def br(x,y):
    out={}
    for i,xi in x.items():
        for j,yj in y.items():
            for k,c in BR[(i,j)].items():
                out[k]=out.get(k,0)+xi*yj*c
    return {k:v for k,v in out.items() if v!=0}
def add(x,y,c=1):
    out=dict(x)
    for k,v in y.items(): out[k]=out.get(k,0)+c*v
    return {k:v for k,v in out.items() if v!=0}
def scal(c,x): return {k:c*v for k,v in x.items() if c*v!=0}
def vec(i): return {i:Fr(1)}
# theta
def a7(i): return roots[i-7][6] if i>=7 else 0
hidx=[i for i in range(DIM) if i<7 or a7(i)%2==0]
pidx=[i for i in range(DIM) if i>=7 and a7(i)%2!=0]
# invariant form
def B(x,y):
    s=Fr(0)
    for i,xi in x.items():
        for j,yj in y.items():
            if i<7 and j<7: s+=xi*yj*C[i][j]
            elif i>=7 and j>=7:
                a=roots[i-7]; b=roots[j-7]
                if all(a[k]+b[k]==0 for k in range(n)): s+=-xi*yj
    return s
if __name__=="__main__":
    print("dims h,p:",len(hidx),len(pidx))
    # Jacobi on all triples (sample fully: 133^3/6 ~ 390k)  -- do full
    import random
    bad=0; cnt=0
    for i in range(DIM):
        for j in range(i+1,DIM):
            for k in range(j+1,DIM):
                t=add(add(br(vec(i),br(vec(j),vec(k))),br(vec(j),br(vec(k),vec(i)))),br(vec(k),br(vec(i),vec(j))))
                cnt+=1
                if t: bad+=1
    print("Jacobi triples checked:",cnt," failures:",bad)
    # antisymmetry
    bad=0
    for i in range(DIM):
        for j in range(DIM):
            if add(br(vec(i),vec(j)),br(vec(j),vec(i))): bad+=1
    print("antisymmetry failures:",bad)
    # invariance of B on random triples
    bad=0
    for _ in range(20000):
        i,j,k=[random.randrange(DIM) for _ in range(3)]
        if B(br(vec(i),vec(j)),vec(k))!=B(vec(i),br(vec(j),vec(k))): bad+=1
    print("B-invariance failures (20000 random triples):",bad)
    # theta is an automorphism: check [p,p] in h, [h,p] in p
    bad=0
    for i in pidx:
        for j in pidx:
            if any(k in pidx for k in br(vec(i),vec(j))): bad+=1
    for i in hidx:
        for j in pidx:
            if any(k in hidx for k in br(vec(i),vec(j))): bad+=1
    print("grading failures:",bad)
