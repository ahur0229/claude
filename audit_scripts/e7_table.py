# Independent reconstruction of the EVII orbit trace table T(e) and Kazhdan weights from characteristics.
# Bourbaki E7 diagram: edges 13,34,45,56,67,24 ; theta(E_alpha)=(-1)^{a_7} E_alpha.
import itertools
edges=[(1,3),(3,4),(4,5),(5,6),(6,7),(2,4)]
n=7
C=[[0]*n for _ in range(n)]
for i in range(n): C[i][i]=2
for a,b in edges: C[a-1][b-1]=-1; C[b-1][a-1]=-1
# generate positive roots in simple-root coords
roots=set()
simple=[tuple(1 if j==i else 0 for j in range(n)) for i in range(n)]
frontier=list(simple); roots=set(simple)
def ip(r,i): return sum(C[i][j]*r[j] for j in range(n))  # <r, alpha_i^vee>
while frontier:
    new=[]
    for r in frontier:
        for i in range(n):
            # reflect: r - <r,alpha_i^vee> alpha_i
            k=ip(r,i)
            rr=list(r); rr[i]-=k; rr=tuple(rr)
            if all(x>=0 for x in rr) and any(x>0 for x in rr) and rr not in roots:
                roots.add(rr); new.append(rr)
    frontier=new
pos=sorted(roots,key=lambda r:(sum(r),r))
print("number of positive roots:",len(pos))
allroots=pos+[tuple(-x for x in r) for r in pos]
def weight(r,a): return sum(r[i]*a[i] for i in range(n))
# characteristics (first six = alpha_1..alpha_6 values, last = alpha_7 value)
chars={1:(1,0,0,0,0,0,0),2:(0,0,0,0,0,1,-2),3:(0,0,0,0,0,1,0),4:(1,0,0,0,0,0,-2),5:(1,0,0,0,0,1,-2),
6:(0,0,0,0,0,0,2),7:(0,0,0,0,0,0,-2),8:(0,0,0,0,0,2,-2),9:(2,0,0,0,0,0,-2),10:(0,2,0,0,0,0,-2),
11:(0,1,0,0,1,0,-2),12:(0,1,1,0,0,0,-3),13:(3,0,0,0,0,1,-2),14:(1,0,0,0,0,3,-6),15:(2,0,0,0,0,2,-4),
16:(2,0,0,0,0,2,-2),17:(4,0,0,0,0,0,-2),18:(0,0,0,0,0,4,-6),19:(2,0,0,0,0,2,-6),20:(2,2,0,0,0,2,-6),
21:(4,0,0,0,0,4,-6),22:(4,0,0,0,0,4,-10)}
from collections import Counter
for k,a in chars.items():
    h=Counter(); p=Counter()
    h[0]+=7  # Cartan
    for r in allroots:
        w=weight(r,a)
        if r[6]%2==0: h[w]+=1
        else: p[w]+=1
    # check: d in h requires weights integer; dimension checks
    assert sum(h.values())==79 and sum(p.values())==54
    # multiplicities of highest weights: m_j^h = h_j - p_{j+2}, m_j^p = p_j - h_{j+2}
    mh={j:h[j]-p[j+2] for j in range(0,40) if h[j]-p[j+2]!=0}
    mp={j:p[j]-h[j+2] for j in range(0,40) if p[j]-h[j+2]!=0}
    assert all(v>=0 for v in mh.values()) and all(v>=0 for v in mp.values()), (k,mh,mp)
    T=sum(j*m for j,m in mh.items())
    kaz=sorted([(j+2,m) for j,m in mp.items()])
    dimpe=sum(mp.values()); dimhe=sum(mh.values())
    # check kernel-difference: dim p^e - dim h^e = 54-79 = -25
    assert dimpe-dimhe==-25,(k,dimpe,dimhe)
    orbdim=54-dimpe
    print(k,a,"T(e)=",T,"Kazhdan:",",".join(f"{w}^{m}" for w,m in kaz),"dim orbit=",orbdim, "T-54=",T-54)
