# Pleasantness test for inner involutions theta=Ad(exp(pi i w_k^vee)) of simply connected E7:
# (G_ad)^theta is disconnected iff exists a root beta with (beta - theta_high)/2 in P^vee \ Q^vee  (see derivation).
import sys; sys.path.insert(0,'.')
from e7lie import C, roots, pos
import itertools
n=7
# highest root
theta_high=max(pos,key=sum)
print("highest root:",theta_high)
# fundamental coweights in simple coroot coordinates: inverse Cartan matrix
import sympy as sp
Cm=sp.Matrix(C); Ci=Cm.inv()
def test(k):
    # t = exp(pi i w_k^vee); w(t) = exp(pi i w(w_k^vee)); need w(w_k) - w_k - 2 w_7 in 2 Q^vee  (z=-1 = exp(2 pi i w_7^vee))
    wk=Ci[:,k-1]; w7=Ci[:,6]
    # W-orbit of w_k^vee: compute by reflecting
    orbit={tuple(wk)}; frontier=[tuple(wk)]
    while frontier:
        new=[]
        for v in frontier:
            vv=sp.Matrix(v)
            for i in range(n):
                # reflection s_i(v) = v - <v,alpha_i> alpha_i^vee ; <v,alpha_i> = (C v)_i
                c=(Cm*vv)[i]
                r=vv.copy(); r[i]-=c
                t=tuple(r)
                if t not in orbit: orbit.add(t); new.append(t)
        frontier=new
    found=False
    for v in orbit:
        diff=sp.Matrix(v)-wk-2*w7
        if all(x.is_integer and x%2==0 for x in diff): found=True; break
    return len(orbit),found
for k in [1,2,7]:
    size,found=test(k)
    print("node",k,": orbit size",size,"; exists w with w(t)=-t :",found," => (G_ad)^theta",("disconnected (not pleasant)" if found else "connected (pleasant)"))
