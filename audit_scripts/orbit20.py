import sys; sys.path.insert(0,'.')
from e7lie import *
import sympy as sp
def V(d):
    v=[0]*DIM
    for k,c in d.items(): v[k]=c
    return v
def parse(expr):
    # expr: dict index->coef
    return {k:Fr(v) for k,v in expr.items()}
e=parse({56:-1,111:1,62:-4,101:4,63:-4,103:4,106:-2,59:2})
f=parse({33:Fr(3,8),38:Fr(-3,4),40:-1,43:Fr(3,2),122:Fr(-3,2),125:Fr(3,4),126:1,127:Fr(-3,8)})
d=parse({0:6,1:8,2:10,3:14,4:10,5:6})
print("[e,f]-d:",add(br(e,f),d,-1))
print("[d,e]-2e:",add(br(d,e),e,-2))
print("[d,f]+2f:",add(br(d,f),f,2))
print("e in p:",all(k in pidx for k in e)," f in p:",all(k in pidx for k in f)," d in h:",all(k in hidx for k in d))
# characteristic alpha_i(d): d = sum c_i H_i ; alpha_i(d)=sum_j C[i][j] c_j
cvec=[d.get(i,0) for i in range(7)]
print("characteristic:",[sum(C[i][j]*cvec[j] for j in range(7)) for i in range(7)])
# ad matrices on p and h
def admat(x, src, tgt):
    M=sp.zeros(len(tgt),len(src))
    tpos={t:a for a,t in enumerate(tgt)}
    for b,s in enumerate(src):
        for k,c in br(x,vec(s)).items():
            M[tpos[k],b]=sp.Rational(c.numerator,c.denominator)
    return M
adf_p=admat(f,pidx,hidx)
ns=adf_p.nullspace()
print("dim p^f =",len(ns))
pf=[{pidx[i]:Fr(int(v[i].p),int(v[i].q)) for i in range(len(pidx)) if v[i]!=0} for v in ns]
# weights of ad d on p^f
for x in pf:
    dx=br(d,x)
    # find lambda with dx = lambda x
    k0=next(iter(x)); lam=dx.get(k0,0)/x[k0]
    assert add(dx,x,-lam)=={}
    print(" p^f vector weight",lam,":",{k:str(v) for k,v in x.items()})
# manuscript normal basis
nb=parse({7:1,132:1}); nu=parse({14:Fr(1,4),20:Fr(-1,2),26:1}); nv=parse({129:-1,130:Fr(-1,2),131:Fr(1,4)})
na=parse({33:-1,38:2,40:Fr(8,3),43:-4,122:4,125:-2,126:Fr(-8,3),127:1})
for name,x in [("n_b",nb),("n_u",nu),("n_v",nv),("n_a",na)]:
    print(name,"[f,n]=0:",br(f,x)=={}, " in p:",all(k in pidx for k in x), " weight:", (lambda dx: dx.get(next(iter(x)),0)/x[next(iter(x))])(br(d,x)))
# triple centralizer in h
ade_h=admat(e,hidx,pidx); adf_h=admat(f,hidx,pidx); add_h=admat(d,hidx,hidx)
M=sp.Matrix.vstack(ade_h,adf_h,add_h)
cent=M.nullspace()
print("dim triple centralizer in h:",len(cent))
centv=[{hidx[i]:Fr(int(v[i].p),int(v[i].q)) for i in range(len(hidx)) if v[i]!=0} for v in cent]
# center of the centralizer: z with [z,c]=0 for all c
# express brackets in terms of centralizer basis: solve linear system
Cmat=sp.Matrix.hstack(*cent)
rows=[]
syms=sp.symbols('c0:%d'%len(centv))
eqs=[]
for cvv in centv:
    # [z, c] where z=sum s_i centv_i
    for i,zi in enumerate(centv):
        pass
# build matrix of ad(centv_i) restricted: map z -> ([z,c_j])_j in coordinates of h
big=[]
for j,cj in enumerate(centv):
    cols=[]
    for i,zi in enumerate(centv):
        b=br(zi,cj); cols.append(sp.Matrix(V(b))[hidx,:] if False else sp.Matrix([b.get(k,0) for k in hidx]))
    big.append(sp.Matrix.hstack(*cols))
BIG=sp.Matrix.vstack(*big)
zs=BIG.nullspace()
print("dim center of triple centralizer:",len(zs))
z=zs[0]; zvec={}
for i,zi in enumerate(centv):
    zvec=add(zvec,zi,Fr(int(z[i].p),int(z[i].q)))
print("center generator:",{k:str(v) for k,v in zvec.items()})
for name,x in [("n_b",nb),("n_u",nu),("n_v",nv),("n_a",na)]:
    zx=br(zvec,x); k0=next(iter(x)); lam=zx.get(k0,0)/x[k0]; assert add(zx,x,-lam)=={}
    print(" center weight on",name,":",lam)
# ---- component sigma_e ----
# Chevalley map C: v_i -> -v_i (Cartan), E_alpha -> E_{-alpha}
def chev(x):
    out={}
    for k,c in x.items():
        if k<7: out[k]=out.get(k,0)-c
        else:
            r=roots[k-7]; out[idx[tuple(-t for t in r)]]=out.get(idx[tuple(-t for t in r)],0)+c
    return {k:v for k,v in out.items() if v!=0}
def expad(x,y,order=20):
    # exp(ad x) y  (ad x nilpotent for root vectors)
    out=dict(y); term=dict(y)
    for m in range(1,order):
        term=scal(Fr(1,m),br(x,term))
        if not term: break
        out=add(out,term)
    return out
def n_i(i,y):
    Ea=vec(idx[simple[i-1]]); Em=vec(idx[tuple(-t for t in simple[i-1])])
    return expad(Ea,expad(Em,expad(Ea,y)))
def w0rep(y):
    word=[1,2,3,4,5,6]*4+[1,2,3,4,5]+[1,2,3,4]+[1,3]+[1]
    for i in word: y=n_i(i,y)   # applied from left to right
    return y
def theta(x): return {k:(c if (k<7 or a7(k)%2==0) else -c) for k,c in x.items()}
def torsign(x):
    out={}
    for k,c in x.items():
        if k<7: out[k]=c
        else: out[k]=c*((-1)**roots[k-7][3])  # alpha_4 coefficient
    return out
def sigma(x): return torsign(theta(w0rep(chev(x))))
print("sigma fixes e:",add(sigma(e),e,-1)=={}," d:",add(sigma(d),d,-1)=={}," f:",add(sigma(f),f,-1)=={})
print("sigma^2=1 on e,nb,nu,nv,na:",[add(sigma(sigma(x)),x,-1)=={} for x in [e,nb,nu,nv,na]])
for name,x in [("n_b",nb),("n_u",nu),("n_v",nv),("n_a",na)]:
    sx=sigma(x)
    print(" sigma(",name,") =", {("n_b" if k==7 else k):str(v) for k,v in sx.items()} if False else None, "equals n_b?",add(sx,nb,-1)=={}, "n_u?",add(sx,nu,-1)=={}, "n_v?",add(sx,nv,-1)=={}, "n_a?",add(sx,na,-1)=={})
# swaps the two 27-spaces?
plus=[i for i in pidx if a7(i)==1]; minus=[i for i in pidx if a7(i)==-1]
ok=all(all(k in minus for k in sigma(vec(i))) for i in plus) and all(all(k in plus for k in sigma(vec(i))) for i in minus)
print("sigma swaps J27 and J27^*:",ok)
# sigma inverts the central torus generator
sz=sigma(zvec); print("sigma(z) = -z:",add(sz,zvec,1)=={})
# sigma preserves B?
import random
bad=0
for _ in range(2000):
    i,j=random.randrange(DIM),random.randrange(DIM)
    if B(sigma(vec(i)),sigma(vec(j)))!=B(vec(i),vec(j)): bad+=1
print("sigma preserves B failures:",bad)
# sigma is a Lie algebra automorphism (sample)
bad=0
for _ in range(3000):
    i,j=random.randrange(DIM),random.randrange(DIM)
    if add(sigma(br(vec(i),vec(j))),br(sigma(vec(i)),sigma(vec(j))),-1): bad+=1
print("sigma automorphism failures:",bad)
# ---- slice invariants ----
b_,u_,v_,a_=sp.symbols('b u v a')
# x = e + b nb + u nu + v nv + a na ; D_x = ad(x+) ad(x-) | J27 (plus)
def tosym(x): return {k:sp.Rational(c.numerator,c.denominator) for k,c in x.items()}
def symadd(*terms):
    out={}
    for coef,x in terms:
        for k,c in x.items(): out[k]=out.get(k,0)+coef*sp.Rational(c.numerator,c.denominator)
    return out
x=symadd((1,e),(b_,nb),(u_,nu),(v_,nv),(a_,na))
xp={k:c for k,c in x.items() if a7(k)==1}; xm={k:c for k,c in x.items() if a7(k)==-1}
def adsym(x, src, tgt):
    M=sp.zeros(len(tgt),len(src)); tpos={t:a for a,t in enumerate(tgt)}
    for b,s in enumerate(src):
        for k,c in x.items():
            for kk,cc in BR[(k,s)].items():
                M[tpos[kk],b]+=c*sp.Rational(cc.numerator,cc.denominator)
    return M
Dx=adsym(xp,plus,minus)*adsym(xm,minus,plus) if False else None
# ad(x-): J27 -> J27^*?? [x_-, J27] lands in h; careful: D_x = ad(x+) ad(x-) restricted to J27 means ad(x-): J27->h, ad(x+): h->J27
A1=adsym(xm,plus,hidx)   # ad(x-): J27 -> h
A2=adsym(xp,hidx,plus)   # ad(x+): h -> J27
Dx=sp.expand(A2*A1)
print("tr D_x =",sp.expand(Dx.trace()))
print("tr D_x^2 =",sp.expand((Dx*Dx).trace()))
