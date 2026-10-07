import sys; sys.path.insert(0,'.')
from e7lie import *
import sympy as sp, itertools
def parse(d): return {k:Fr(v) for k,v in d.items()}
e=parse({56:-1,111:1,62:-4,101:4,63:-4,103:4,106:-2,59:2})
nb=parse({7:1,132:1}); nu=parse({14:Fr(1,4),20:Fr(-1,2),26:1}); nv=parse({129:-1,130:Fr(-1,2),131:Fr(1,4)})
na=parse({33:-1,38:2,40:Fr(8,3),43:-4,122:4,125:-2,126:Fr(-8,3),127:1})
plus=[i for i in pidx if a7(i)==1]; minus=[i for i in pidx if a7(i)==-1]
print("dim J27, J27*:",len(plus),len(minus))
def wt6(i):  # E6 weight vector (alpha_1..alpha_6)(E_alpha) for root index i
    r=roots[i-7]; return tuple(sum(C[k][j]*r[j] for j in range(7)) for k in range(6))
def invariant_cubic(space):
    pos={i:a for a,i in enumerate(space)}
    # weight-zero cubic monomials (multisets of 3 indices)
    monos=[]
    for m in itertools.combinations_with_replacement(space,3):
        w=tuple(sum(wt6(i)[k] for i in m) for k in range(6))
        if all(x==0 for x in w): monos.append(m)
    print(" weight-zero cubic monomials:",len(monos))
    cs=sp.symbols('c0:%d'%len(monos))
    xs=sp.symbols('x0:%d'%len(space))
    Pol=sum(c*sp.prod([xs[pos[i]] for i in m]) for c,m in zip(cs,monos))
    eqs=[]
    for i in range(1,7):
        for sgn in (1,-1):
            r=tuple(sgn*t for t in simple[i-1]); X=vec(idx[r])
            # derivation: (X.P)(x) = sum_j (X x)_j dP/dx_j where X acts on the vector x = sum x_j E_j: X x = sum_j x_j [X,E_j]
            der=0
            for j in space:
                b=br(X,vec(j))
                for k,c in b.items():
                    if k in pos: der+= xs[pos[j]]*sp.Rational(c.numerator,c.denominator)*sp.diff(Pol,xs[pos[k]])
                    else:
                        assert False,"left the space"
            der=sp.expand(der)
            eqs+= sp.Poly(der,*xs).coeffs()
    sol=sp.linsolve(eqs,cs)
    sol=list(sol)[0]
    free=[s for s in cs if any(s in sp.sympify(v).free_symbols for v in sol)]
    print(" dimension of invariant cubic space:",len(free))
    Pinv=Pol.subs(dict(zip(cs,sol)))
    return Pinv, xs, pos, free
Np,xs,pos,free=invariant_cubic(plus)
Nm,ys,posm,freem=invariant_cubic(minus)
b_,u_,v_,a_=sp.symbols('b u v a')
def slice_sub(xsyms,pos_,parts):
    x={}
    for coef,vecd in parts:
        for k,c in vecd.items():
            if k in pos_: x[k]=x.get(k,0)+coef*sp.Rational(c.numerator,c.denominator)
    return {xsyms[pos_[k]]:v for k,v in x.items()}
subp=slice_sub(xs,pos,[(1,e),(b_,nb),(u_,nu),(v_,nv),(a_,na)])
subm=slice_sub(ys,posm,[(1,e),(b_,nb),(u_,nu),(v_,nv),(a_,na)])
Np_slice=sp.expand(Np.subs({x:0 for x in xs}|subp)); Nm_slice=sp.expand(Nm.subs({y:0 for y in ys}|subm))
print("N_+ on slice (up to scalar):",sp.factor(Np_slice))
print("N_- on slice (up to scalar):",sp.factor(Nm_slice))
# stabilizer dimensions on the branches e+n_u, e+n_v, and at e
def stabdim(x):
    M=sp.zeros(DIM,len(hidx))
    for b,j in enumerate(hidx):
        for k,c in br(x,vec(j)).items(): M[k,b]=sp.Rational(c.numerator,c.denominator)
    return len(M.nullspace())
print("dim h^{e} =",stabdim(e)," dim stabilizer at e+n_u:",stabdim(add(e,nu))," at e+n_v:",stabdim(add(e,nv)))
