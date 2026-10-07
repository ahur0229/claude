# Independent check of Appendix A E7/A7 table: Cartan matrix with edges 12,23,34,45,56,37.
import itertools
from collections import Counter
edges=[(1,2),(2,3),(3,4),(4,5),(5,6),(3,7)]
n=7
C=[[0]*n for _ in range(n)]
for i in range(n): C[i][i]=2
for a,b in edges: C[a-1][b-1]=-1; C[b-1][a-1]=-1
simple=[tuple(1 if j==i else 0 for j in range(n)) for i in range(n)]
roots=set(simple); frontier=list(simple)
while frontier:
    new=[]
    for r in frontier:
        for i in range(n):
            k=sum(C[i][j]*r[j] for j in range(n))
            rr=list(r); rr[i]-=k; rr=tuple(rr)
            if rr not in roots and any(x!=0 for x in rr):
                roots.add(rr); new.append(rr)
    frontier=new
R=sorted(roots)
print("roots:",len(R))
res={}
cands=0; cand2=0
for a in itertools.product([0,1,2],repeat=7):
    c=Counter()
    for r in R: c[sum(r[i]*a[i] for i in range(7))]+=1
    c[0]+=7
    m={j:c[j]-c[j+2] for j in range(0,60)}
    if any(v<0 for v in m.values()): continue
    if any(j%2==1 and m[j]%2==1 for j in m): continue
    cands+=1
    p=sum(m.values()); q=(p+7)//2
    assert (p+7)%2==0
    pos=[]
    for j in sorted(m):
        if j>0: pos+= [j]*m[j]
    if len(pos)<q: continue
    cand2+=1
    b=sum(v+2 for v in sorted(pos)[:q])-70
    res.setdefault(q,[]).append(b)
print("candidates after sign tests:",cands," after distinguished count:",cand2)
for q in sorted(res): print(q,len(res[q]),min(res[q]))
print("global min b:",min(min(v) for v in res.values()))
