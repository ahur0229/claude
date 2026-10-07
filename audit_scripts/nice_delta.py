# Numeric spot-checks of Appendix A: compute delta_p(e) for chain configurations in equal AIII, CI, BDI.
import sympy as sp, itertools
def chain(p):
    E=sp.zeros(p); D=sp.zeros(p); F=sp.zeros(p)
    for t in range(p):
        if t-1>=0: E[t-1,t]=1       # e v_t = v_{t-1}
        D[t,t]=p-1-2*t
        if t+1<p: F[t+1,t]=(t+1)*(p-1-t)
    return E,D,F
def delta(E,D,F,J,pbasis_filter):
    """E,D,F matrices on V, J grading; Lie algebra g given by list of basis matrices; p = odd part under conj by J (or given filter)."""
    pass
def vec(M): return sp.Matrix(list(M))
def sl2_check(E,D,F):
    assert D*E-E*D==2*E and D*F-F*D==-2*F and E*F-F*E==D
def delta_from(gbasis, E, D, F, J, sign_p):
    """gbasis: list of matrices spanning g; p = {X in g: J X J^{-1} = sign_p X}... with sign_p = -1 for odd."""
    N=E.shape[0]
    # p basis
    pb=[]
    for X in gbasis:
        Y=(X - J*X*J.inv())/2
        if Y!=sp.zeros(N): pb.append(Y)
    P=sp.Matrix.hstack(*[vec(Y) for Y in pb]); cols=P.columnspace(); pb=[sp.Matrix(N,N,list(c)) for c in cols]
    dimp=len(pb)
    M=sp.Matrix.hstack(*[vec(E*Y-Y*E) for Y in pb]); ns=M.nullspace()
    pe=[sum((c*Y for c,Y in zip(list(v),pb)),sp.zeros(N)) for v in ns]
    PE=sp.Matrix.hstack(*[vec(Y) for Y in pe]); A=sp.zeros(len(pe))
    for j,Y in enumerate(pe):
        w=vec(D*Y-Y*D); A[:,j]=sp.Matrix(sp.linsolve((PE,w)).args[0])
    tr=A.trace()
    # distinguished? p^{e,d,f}=0
    M2=sp.Matrix.vstack(sp.Matrix.hstack(*[vec(E*Y-Y*E) for Y in pb]),sp.Matrix.hstack(*[vec(D*Y-Y*D) for Y in pb]),sp.Matrix.hstack(*[vec(F*Y-Y*F) for Y in pb]))
    dist=(len(M2.nullspace())==0)
    return tr+2*len(pe)-dimp, dist, dimp, len(pe)
# ---- equal AIII: V=V+ + V-, g=sl(V), p = odd part. chains with signs.
def aiii(config):
    # config: list of (p, eps)
    blocks=[chain(p) for p,_ in config]
    E=sp.diag(*[b[0] for b in blocks]); D=sp.diag(*[b[1] for b in blocks]); F=sp.diag(*[b[2] for b in blocks])
    Jd=[]
    for (p,eps) in config: Jd+= [eps*(-1)**t for t in range(p)]
    J=sp.diag(*Jd); N=E.shape[0]
    sl2_check(E,D,F)
    gb=[]
    for i in range(N):
        for j in range(N):
            X=sp.zeros(N); X[i,j]=1
            if i==j: X[N-1,N-1]-=1
            if X!=sp.zeros(N): gb.append(X)
    nplus=sum(1 for x in Jd if x==1); nminus=N-nplus
    return delta_from(gb,E,D,F,J,-1), (nplus,nminus)
print("== equal AIII ==")
for config in [[(2,1),(2,1)],[(2,1),(2,-1)],[(3,1),(1,-1)],[(3,1),(3,-1)],[(4,1),(2,-1)],[(4,1),(2,1)],[(3,1),(2,1),(1,-1)],[(5,1),(1,-1)],[(3,1),(1,-1),(1,-1),(1,1)]]:
    r,(a,b)=aiii(config)
    if a==b: print(config,"delta=",r[0],"distinguished=",r[1],"dim p=",r[2],"dim p^e=",r[3])
# ---- CI: V symplectic dim 2n with form b_m on chains, J grading with B(Ju,Jv)=-B(u,v)
def chain_form(p):
    b=sp.zeros(p)
    for i in range(p): b[i,p-1-i]=(-1)**i
    return b
def iso_alg(Omega):
    N=Omega.shape[0]; Oi=Omega.inv(); gb=[]
    for i in range(N):
        for j in range(i,N):
            S=sp.zeros(N); S[i,j]=1; S[j,i]=1   # symmetric -> for symplectic form: X = Oinv S ; for symmetric form need skew S
            gb.append(Oi*S)
    return gb
def iso_alg_sym(Omega):
    N=Omega.shape[0]; Oi=Omega.inv(); gb=[]
    for i in range(N):
        for j in range(i+1,N):
            S=sp.zeros(N); S[i,j]=1; S[j,i]=-1
            gb.append(Oi*S)
    return gb
def ci(config):
    # config: list of (size 2a, eps) with even sizes (plus possible odd sizes in pairs - skip)
    blocks=[chain(p) for p,_ in config]
    E=sp.diag(*[b[0] for b in blocks]); D=sp.diag(*[b[1] for b in blocks]); F=sp.diag(*[b[2] for b in blocks])
    Om=sp.diag(*[chain_form(p) for p,_ in config])
    assert Om.T==-Om
    for X in (E,D,F): assert X.T*Om+Om*X==sp.zeros(Om.shape[0])
    Jd=[]
    for (p,eps) in config: Jd+=[eps*(-1)**t for t in range(p)]
    J=sp.diag(*Jd)
    assert J.T*Om*J==-Om, "J must be anti-symplectic"
    gb=iso_alg(Om)
    return delta_from(gb,E,D,F,J,-1)
print("== CI == (formula: sum a_i(a_i+1) + sum_{i<j, eps equal} min(2a_i,2a_j))")
for config in [[(2,1)],[(4,1)],[(2,1),(2,1)],[(2,1),(2,-1)],[(4,1),(2,1)],[(4,1),(2,-1)],[(6,1)],[(4,1),(4,-1)],[(2,1),(2,1),(2,1)]]:
    r=ci(config); print(config,"delta=",r[0],"distinguished=",r[1],"dim p=",r[2])
# ---- BDI: orthogonal, chains odd sizes, J isometric
def bdi(config):
    blocks=[chain(p) for p,_ in config]
    E=sp.diag(*[b[0] for b in blocks]); D=sp.diag(*[b[1] for b in blocks]); F=sp.diag(*[b[2] for b in blocks])
    Om=sp.diag(*[chain_form(p) for p,_ in config])
    assert Om.T==Om
    for X in (E,D,F): assert X.T*Om+Om*X==sp.zeros(Om.shape[0])
    Jd=[]
    for (p,eps) in config: Jd+=[eps*(-1)**t for t in range(p)]
    J=sp.diag(*Jd)
    assert J.T*Om*J==Om
    gb=iso_alg_sym(Om)
    r=delta_from(gb,E,D,F,J,-1)
    nplus=sum(1 for x in Jd if x==1); return r,(nplus,len(Jd)-nplus)
print("== BDI == (formula: sum A_i(A_i+2i-k)+sum B_j(B_j+2j+k)+rs)")
for config in [[(3,1)],[(3,1),(1,-1)],[(3,1),(1,1)],[(5,1),(1,-1)],[(3,1),(3,-1)],[(5,1),(3,-1)],[(5,1),(1,1)],[(3,1),(1,-1),(1,-1)],[(3,1),(3,1)],[(7,1),(1,-1)],[(5,1),(3,1)],[(3,1),(1,-1),(1,1)]]:
    r,(a,b)=bdi(config); print(config,"blocks",(a,b),"k=",a-b,"delta=",r[0],"distinguished=",r[1],"dim p=",r[2])
