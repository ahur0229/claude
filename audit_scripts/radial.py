import sys; sys.path.insert(0,'.')
from e7lie import *
import sympy as sp, numpy as np, itertools
PR=33554393  # prime < 2^25 to avoid int64 overflow in 54-term and 1849-term sums
def modinv(a): return pow(int(a)%PR,PR-2,PR)
def fr2mod(x):
    x=Fr(x); return (x.numerator%PR)*modinv(x.denominator)%PR
def parse(d): return {k:Fr(v) for k,v in d.items()}
R=lambda a,b: Fr(a,b)
e17=parse({40:-1,45:-1,49:1,53:-1,54:1,56:1,57:1,59:1,60:-1,62:-1,63:1,64:-1,65:1,66:1,67:-1,68:1,70:1,77:-1,83:-1,89:1,95:-1,96:-1,101:-1,106:-1,111:1,115:-1})
d17=parse({0:6,1:5,2:8,3:10,4:7,5:4,6:1})
f17={7:R(-4,3),14:R(1,3),20:R(2,3),26:R(2,3),32:R(1,3),33:R(4,3),38:R(1,3),43:R(-1,3),48:R(-4,3),52:R(1,3),103:R(1,2),108:R(-3,2),112:R(-3,2),116:R(3,2),117:R(-1,2),119:R(-1,2),120:R(-1,2),122:R(-3,2),123:R(-1,2),125:R(-3,2),126:R(-3,2),127:R(1,2),128:R(-1,2),129:R(-1,2),130:R(-1,2),131:R(-1,2)}
e18=parse({52:-1,56:-1,59:1,62:-1,64:1,65:1,66:1,67:1,68:-1,69:-1,77:1,83:-1,89:1,95:1,96:-1,101:1,103:1,106:-1,108:-1,111:1,112:-1,116:-1,117:-1,120:-1,123:1,126:-1})
d18=parse({0:2,1:3,2:4,3:6,4:5,5:4,6:-1})
f18={14:R(-9,10),20:R(3,10),26:R(-3,10),32:R(-9,10),33:R(3,10),38:R(3,10),40:R(-9,10),43:R(9,10),45:R(-3,10),48:R(-3,2),49:R(-3,10),53:R(3,10),54:R(3,10),57:R(-3,10),60:R(3,10),63:R(9,10),115:R(7,5),119:R(1,5),122:R(-4,5),125:R(4,5),127:R(-1,5),128:R(2,5),129:R(-1,5),130:R(-7,5),131:R(4,5),132:R(-2,5)}
orbit=int(sys.argv[1])
e,d,f={17:(e17,d17,f17),18:(e18,d18,f18)}[orbit]
print("orbit",orbit,": [e,f]=d:",add(br(e,f),d,-1)=={}," [d,e]=2e:",add(br(d,e),e,-2)=={}," [d,f]=-2f:",add(br(d,f),f,2)=={})
cvec=[d.get(i,0) for i in range(7)]
chars=[sum(C[i][j]*cvec[j] for j in range(7)) for i in range(7)]
print(" characteristic:",chars)
def wt(i): return sum(roots[i-7][k]*chars[k] for k in range(7)) if i>=7 else 0
# ---- normal basis via RREF nullspace of ad f per weight (increasing weight), p basis in increasing index order
def admat_sym(x,src,tgt):
    M=sp.zeros(len(tgt),len(src)); tpos={t:a for a,t in enumerate(tgt)}
    for b,s in enumerate(src):
        for k,c in br(x,vec(s)).items(): M[tpos[k],b]=sp.Rational(c.numerator,c.denominator)
    return M
weights=sorted(set(wt(i) for i in pidx))
normal=[]
for w in weights:
    src=[i for i in pidx if wt(i)==w]
    M=admat_sym(f,src,hidx)
    for v in M.nullspace():
        normal.append(({src[i]:Fr(int(v[i].p),int(v[i].q)) for i in range(len(src)) if v[i]!=0},w))
print(" normal basis weights:",[w for _,w in normal], " count:",len(normal))
nvecs=[x for x,_ in normal]
kaz=[w+2 for _,w in normal]; kaz=[-w+2 for _,w in normal]  # Kazhdan weight = -d-weight + 2
print(" Kazhdan weights:",kaz)
# ---- pivot columns of ad e : h -> p
Me=admat_sym(e,hidx,pidx)
rref,piv=Me.rref()
mvecs=[vec(hidx[c]) for c in piv]
print(" dim [h,e] =",len(piv))
assert len(piv)==54-len(nvecs)
# h^e basis
he=[{hidx[i]:Fr(int(v[i].p),int(v[i].q)) for i in range(len(hidx)) if v[i]!=0} for v in Me.nullspace()]
print(" dim h^e =",len(he))
# ---- mod-p linear algebra
ppos={t:a for a,t in enumerate(pidx)}
def pcol(x):
    col=np.zeros(54,dtype=np.int64)
    for k,c in x.items():
        assert k in ppos
        col[ppos[k]]=fr2mod(c)
    return col
nm=len(mvecs); nn=len(nvecs)
W0=np.zeros((54,54),dtype=np.int64)
for a,m in enumerate(mvecs): W0[:,a]=pcol(br(m,e))
for j,nv in enumerate(nvecs): W0[:,nm+j]=pcol(nv)
Wk=[]
for k,nv in enumerate(nvecs):
    Wkk=np.zeros((54,54),dtype=np.int64)
    for a,m in enumerate(mvecs): Wkk[:,a]=pcol(br(m,nv))
    Wk.append(Wkk)
def mm(A,B): return (A@B)%PR
def inv_mod(A):
    n=A.shape[0]; M=np.concatenate([A%PR,np.eye(n,dtype=np.int64)],axis=1)
    for c in range(n):
        piv=next(r for r in range(c,n) if M[r,c]%PR!=0)
        M[[c,piv]]=M[[piv,c]]
        inv=modinv(M[c,c]); M[c]=(M[c]*inv)%PR
        for r in range(n):
            if r!=c and M[r,c]!=0: M[r]=(M[r]-M[r,c]*M[c])%PR
    return M[:,n:]
U0=inv_mod(W0)
assert np.array_equal(mm(U0,W0),np.eye(54,dtype=np.int64))
# series U(s) = sum_m U0 (-W1 U0)^m, keyed by ordered index sequences
DEG=4
Useq={():U0}
for m in range(1,DEG+1):
    for seq in itertools.product(range(nn),repeat=m):
        prev=Useq[seq[:-1]]
        Useq[seq]=mm(prev,(-mm(Wk[seq[-1]],U0))%PR)
# collect into monomials (sorted tuples)
from collections import defaultdict
def key(seq): return tuple(sorted(seq))
Umon=defaultdict(lambda: np.zeros((54,54),dtype=np.int64))
for seq,M in Useq.items(): Umon[key(seq)]=(Umon[key(seq)]+M)%PR
Pmon={k:M[nm:,:] for k,M in Umon.items()}
# Q = inverse ambient form on p (in pidx basis): B(E_a,E_-a)=-1 ; matrix is its own inverse
Bp=np.zeros((54,54),dtype=np.int64)
for i in pidx:
    for j in pidx:
        if all(roots[i-7][k]+roots[j-7][k]==0 for k in range(7)): Bp[ppos[i],ppos[j]]=PR-1
assert np.array_equal(mm(Bp,Bp),np.eye(54,dtype=np.int64))
Q=Bp
def monmul(A,B):  # product of monomial-keyed matrix series, truncated at DEG
    out=defaultdict(lambda: np.zeros(A[next(iter(A))].shape[0] if False else 0))
    out={}
    for ka,Ma in A.items():
        for kb,Mb in B.items():
            if len(ka)+len(kb)>DEG: continue
            k=tuple(sorted(ka+kb))
            out[k]=(out.get(k,0)+mm(Ma,Mb))%PR
    return out
Pt={k:M.T.copy() for k,M in Pmon.items()}
g=monmul(monmul(Pmon,{():Q}),Pt)
# check degrees: g_{ij} should have Kazhdan degree w_i+w_j-4 and ordinary degree <=3
maxdeg=max(len(k) for k,M in g.items() if np.any(M%PR))
print(" max ordinary degree in g:",maxdeg)
# Kazhdan homogeneity check
bad=0
for k,M in g.items():
    kd=sum(kaz[t] for t in k)
    for i in range(nn):
        for j in range(nn):
            if M[i,j]%PR and kd!=kaz[i]+kaz[j]-4: bad+=1
print(" Kazhdan-inhomogeneous g entries:",bad)
# ---- h
K=monmul(monmul(Umon,{():Q}),{k:M.T.copy() for k,M in Umon.items()})
# second derivative terms: S(s) = sum_{a,b} K_ab [m_a,[m_b,e+s]] + 2 sum_{a,j} K_{a,nm+j} [m_a,n_j]   (vector in p, series in s)
# [m_a,[m_b,e+s]] = [m_a,[m_b,e]] + sum_k s_k [m_a,[m_b,n_k]]
BB0=np.zeros((54,nm,nm),dtype=np.int64); BBk=np.zeros((nn,54,nm,nm),dtype=np.int64)
for a,ma in enumerate(mvecs):
    for b,mb in enumerate(mvecs):
        BB0[:,a,b]=pcol(br(ma,br(mb,e)))
        for k,nv in enumerate(nvecs): BBk[k,:,a,b]=pcol(br(ma,br(mb,nv)))
MN=np.zeros((54,nm,nn),dtype=np.int64)
for a,ma in enumerate(mvecs):
    for j,nv in enumerate(nvecs): MN[:,a,j]=pcol(br(ma,nv))
S={}
HDEG=2
for kk,KM in K.items():
    if len(kk)>HDEG: continue
    Kab=KM[:nm,:nm]; Kan=KM[:nm,nm:]
    v=np.einsum('iab,ab->i',BB0,Kab)%PR + 2*np.einsum('iaj,aj->i',MN,Kan)%PR
    S[kk]=(S.get(kk,0)+v)%PR
    for k in range(nn):
        if len(kk)+1>HDEG: continue
        k2=tuple(sorted(kk+(k,)))
        S[k2]=(S.get(k2,0)+np.einsum('iab,ab->i',BBk[k],Kab))%PR
h={}
for kp,PM in Pmon.items():
    for ks,Sv in S.items():
        if len(kp)+len(ks)>HDEG: continue
        k=tuple(sorted(kp+ks))
        h[k]=(h.get(k,0)-(PM@Sv))%PR
maxdegh=max(len(k) for k,v in h.items() if np.any(v%PR))
print(" max ordinary degree in h:",maxdegh)
bad=0
for k,v in h.items():
    kd=sum(kaz[t] for t in k)
    for i in range(nn):
        if v[i]%PR and kd!=kaz[i]-4: bad+=1
print(" Kazhdan-inhomogeneous h entries:",bad)
# ---- jet action. polynomials in x_0..x_10: dict exponent-tuple -> coeff mod p
def padd(P,Qp,c=1):
    out=dict(P)
    for k,v in Qp.items(): out[k]=(out.get(k,0)+c*v)%PR
    return {k:v for k,v in out.items() if v}
def pmul_x(P,i):
    out={}
    for k,v in P.items():
        kk=list(k); kk[i]+=1; out[tuple(kk)]=v
    return out
def pdiff(P,i):
    out={}
    for k,v in P.items():
        if k[i]>0:
            kk=list(k); kk[i]-=1; out[tuple(kk)]=(out.get(tuple(kk),0)+v*k[i])%PR
    return {k:v for k,v in out.items() if v}
def apply_op(coeffmons,P):
    """coeffmons: dict monomial(sorted tuple of s-indices)->scalar; operator sum c s^gamma with s_k -> -d/dx_k"""
    out={}
    for gam,c in coeffmons.items():
        if not c: continue
        Qp=P
        for k in gam: Qp=pdiff(Qp,k)
        out=padd(out,Qp,(c*((-1)**len(gam)))%PR)
    return out
def Rop(F):
    out={}
    for i in range(nn):
        for j in range(nn):
            cm={k:int(M[i,j]) for k,M in g.items() if M[i,j]%PR}
            if cm: out=padd(out,apply_op(cm,pmul_x(pmul_x(F,i),j)))
        cm={k:int(v[i]) for k,v in h.items() if v[i]%PR}
        if cm: out=padd(out,apply_op(cm,pmul_x(F,i)))
    return out
zero=tuple([0]*nn)
def mono(i):
    k=[0]*nn; k[i]=1; return tuple(k)
# sanity: R q0 = 108 where q0(s)=B(e+s,e+s): as a polynomial in s -> encode as operator? Instead check g dq0 = (w_i s_i): skip.
# ---- stabilizer vector fields v_z(s)=P(s)[z,e+s]
VDEG=3
vz=[]
for z in he:
    ze=pcol(br(z,e)); zn=[pcol(br(z,nv)) for nv in nvecs]
    v={}
    for kp,PM in Pmon.items():
        if len(kp)>VDEG: continue
        v[kp]=(v.get(kp,0)+PM@ze)%PR
        for k in range(nn):
            if len(kp)+1>VDEG: continue
            k2=tuple(sorted(kp+(k,)))
            v[k2]=(v.get(k2,0)+PM@zn[k])%PR
    vz.append(v)
print(" max degree in stabilizer fields:",max(len(k) for v in vz for k,val in v.items() if np.any(val%PR)))
# linear system on F = sum_{j=1}^{10} c_j x_j
cols=[]
for j in range(1,nn):
    F={mono(j):1}
    img=[]
    for v in vz:
        out={}
        for i in range(nn):
            cm={k:int(val[i]) for k,val in v.items() if val[i]%PR}
            if cm: out=padd(out,apply_op(cm,pmul_x(F,i)))
        img.append(out)
    cols.append(img)
# assemble matrix: rows indexed by (z-index, monomial)
keys=sorted({(zi,k) for col in cols for zi,out in enumerate(col) for k in out})
kpos={k:a for a,k in enumerate(keys)}
A=np.zeros((len(keys),nn-1),dtype=np.int64)
for c,col in enumerate(cols):
    for zi,out in enumerate(col):
        for k,v in out.items(): A[kpos[(zi,k)],c]=v
def nullspace_mod(A):
    M=A.copy()%PR; nrow,ncol=M.shape; pivs=[]; r=0
    for c in range(ncol):
        p=next((i for i in range(r,nrow) if M[i,c]%PR),None)
        if p is None: continue
        M[[r,p]]=M[[p,r]]; M[r]=(M[r]*modinv(M[r,c]))%PR
        for i in range(nrow):
            if i!=r and M[i,c]: M[i]=(M[i]-M[i,c]*M[r])%PR
        pivs.append(c); r+=1
    free=[c for c in range(ncol) if c not in pivs]
    basis=[]
    for fc in free:
        v=np.zeros(ncol,dtype=np.int64); v[fc]=1
        for i,pc in enumerate(pivs): v[pc]=(-M[i,fc])%PR
        basis.append(v)
    return basis
ker=nullspace_mod(A)
print(" dim of connected-stabilizer kernel on J_4:",len(ker))
# compare with manuscript's A,B0
if orbit==17:
    Aman={10:1,2:1,3:3,4:3,5:3,6:-1,7:1,8:1,9:1}; Bman={1:1}
else:
    Aman={10:R(3,2),2:1,3:R(-1,2),4:R(1,2),5:-1,6:R(-3,2),7:-1,9:R(1,2)}; Bman={1:1,10:R(-1,2),3:R(-1,2),4:R(1,2),6:R(1,2),8:-1,9:R(1,2)}
def tovec(dct):
    v=np.zeros(nn-1,dtype=np.int64)
    for j,c in dct.items(): v[j-1]=fr2mod(c)
    return v
def inspan(v,basis):
    M=np.array(basis+[v]); return len(nullspace_mod(M.T))>len(nullspace_mod(np.array(basis).T))
print(" manuscript A in kernel:",inspan(tovec(Aman),ker)," B0 in kernel:",inspan(tovec(Bman),ker))
def R3(dct):
    F={mono(j):fr2mod(c) for j,c in dct.items()}
    for _ in range(3): F=Rop(F)
    return F
RA=R3(Aman); RB=R3(Bman)
# rank of {R^3 A, R^3 B0}
allk=sorted(set(RA)|set(RB)); M=np.array([[RA.get(k,0) for k in allk],[RB.get(k,0) for k in allk]],dtype=np.int64)
print(" rank of R^3 on span{A,B0} (mod p):", 2-len(nullspace_mod(M.T)) if len(allk) else 0, " (nonzero monomials:",len(allk),")")
def coeff(F,expo): 
    k=[0]*nn
    for i,e_ in expo.items(): k[i]=e_
    return F.get(tuple(k),0)
if orbit==17:
    ents=[[coeff(RA,{2:4}),coeff(RA,{1:3,2:1})],[coeff(RB,{2:4}),coeff(RB,{1:3,2:1})]]
    claimed=[[R(3,500),0],[R(286,9),10]]
else:
    ents=[[coeff(RA,{1:4}),coeff(RA,{3:4})],[coeff(RB,{1:4}),coeff(RB,{3:4})]]
    claimed=[[R(-354111,125000),R(-763911,125000)],[R(-192,15625),R(-192,15625)]]
print(" matrix entries mod p:",ents," claimed mod p:",[[fr2mod(x) for x in row] for row in claimed])
