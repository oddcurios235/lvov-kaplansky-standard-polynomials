"""Exact finite checks of the odd standard-polynomial identities.

The sample Jordan routine accepts rational eigenvalues; the article proves its
all-degree results symbolically over algebraically closed characteristic-zero fields.
"""
from __future__ import annotations
import argparse
import json
from collections import defaultdict
from pathlib import Path
import sympy as s
from graded import slice_grade, subset_standard


def em(n: int, i: int, j: int) -> s.Matrix:
    A=s.zeros(n); A[i,j]=1; return A


def standard(args: list[s.Matrix]) -> s.Matrix:
    """Independent full-matrix subset recurrence, including all signs."""
    n=args[0].rows
    dp={0:s.eye(n)}
    for mask in range(1,1<<len(args)):
        ids=[j for j in range(len(args)) if mask>>j&1]
        out=s.zeros(n)
        for pos,j in enumerate(ids):
            out+=(-1)**(len(ids)-1-pos)*dp[mask^(1<<j)]*args[j]
        dp[mask]=out
    return dp[(1<<len(args))-1]


def unit_standard(edges: list[tuple[int,int]], n: int) -> s.Matrix:
    """Independent enumeration of nonzero matrix-unit Euler trails."""
    by=defaultdict(list); out=s.zeros(n); full=(1<<len(edges))-1
    for j,(a,b) in enumerate(edges): by[a].append((j,b))
    def rec(start,cur,mask,sign):
        if mask==full:
            out[start,cur]+=sign; return
        for j,b in by[cur]:
            if not mask>>j&1:
                rec(start,b,mask|(1<<j),sign*(-1)**((mask>>(j+1)).bit_count()))
    for a in range(n):rec(a,a,0,1)
    return out


def hom_matrix(w: list[int], grade: int) -> s.Matrix:
    n=len(w);return s.Matrix(n,n,lambda i,j:w[i] if j==(i+grade)%n else 0)


def periodic(k: int,t: int):
    n=k+t+1
    F=[[int(i<n-1 and i%k==j) for i in range(n)] for j in range(k)]
    R=[[((-1)**((i-1)//k) if i and (i-1)%k==j else 0) for i in range(n)] for j in range(k)]
    return F+R[::-1],[1]*k+[-1]*k


def vandermonde_core(r: int,n: int):
    k=r-1
    F=[[(i+1)**j if i<n-1 else 0 for i in range(n)] for j in range(k)]
    R=[[(-1)**((i-1)//k)*i**j if i else 0 for i in range(n)] for j in range(k)]
    return F+R[::-1],[1]*k+[-1]*k


def positive_formula(r: int,n: int, factor=1):
    c=s.Matrix([1]+[2]*(r-2)+[1]); H=s.diag(*c)+c*c.T
    G=s.zeros(n);b=s.zeros(n,1)
    for a in range(n-r+1):
        for i in range(r):
            b[a+i]+=factor*c[i]
            for j in range(r):G[a+i,a+j]+=factor*H[i,j]
    E=s.diag(*[(-1)**i for i in range(n)])
    return G*E,E*b


def interval_vectors(n,p,q):
    I=[];j=q
    while True:
        I.append(j)
        if j==p:break
        j=(j+1)%n
    a=s.Matrix([1 if j in [p,q] else 2 if j in I else 0 for j in range(n)])
    b=s.Matrix([int(j not in I) for j in range(n)])
    return a,b


def cycle_bivector(D: s.Matrix) -> s.Matrix:
    """Actual signed cofactor bivector for cycle inputs then diagonal rows D."""
    n=D.cols;B=s.zeros(n)
    for i in range(n):
        for j in range(i+1,n):
            cols=[l for l in range(n) if l not in [i,j]]
            minor=D[:,cols].det()
            gamma=(-1)**(n*(n+1)//2+j-i-2)
            B[j,i]=gamma*minor;B[i,j]=-gamma*minor
    return B


def cycle_preimage(blocks: list[tuple[int,int]]):
    """Construct sampled s_(2n-1) inputs for a lower Jordan target.

    ``blocks`` supplies integer/rational eigenvalues and Jordan block sizes.
    Eigenvalues are converted by SymPy Rational, so this checker does not
    implement arbitrary algebraic or complex field inputs. Returns the target,
    ordered input matrices, permutation and auxiliary construction data.
    """
    n=sum(z for _,z in blocks); assert n>=3
    d=[];edges=[];zero_black=[];zero_white=[];nz=[];base=0
    for lam,z in blocks:
        d += [s.Rational(lam)]*z
        edges += [(base+j+1,base+j) for j in range(z-1)]
        if lam==0:
            zero_black += [base+j for j in range(z) if j%2==0]
            zero_white += [base+j for j in range(z) if j%2==1]
        else:nz+=list(range(base,base+z))
        base+=z
    T=sum(d);order=list(range(n))
    if T==0 and nz and zero_white:
        order=zero_black+[nz[0]]+zero_white+nz[1:]
    pos={old:new for new,old in enumerate(order)}
    d=s.Matrix([d[i] for i in order]);edges=[(pos[p],pos[q]) for p,q in edges]
    target=s.diag(*d)
    for p,q in edges:target[p,q]=1
    A=2*s.ones(n)-s.eye(n)
    v=A.inv()*d if d!=s.zeros(n,1) else s.ones(n,1)
    for tau in range((len(edges)+2)*(n-1)+3):
        w=s.Matrix([s.Integer(tau)**i for i in range(n)])
        if s.Matrix.hstack(v,w).rank()<2:continue
        B=v*w.T-w*v.T
        if all((interval_vectors(n,p,q)[0].T*B*interval_vectors(n,p,q)[1])[0]!=0 for p,q in edges):break
    else:raise AssertionError('Polynomial avoidance bound unexpectedly failed')
    rows=s.Matrix.vstack(v.T,w.T).nullspace()
    D=s.Matrix.vstack(*(z.T for z in rows))
    B0=cycle_bivector(D)
    assert B0.rank()==2 and s.Matrix.hstack(B,B0).rank()==2
    sol,params=(A*B0).gauss_jordan_solve(d)
    sol=sol.subs({z:0 for z in params})
    X=s.diag(*sol)
    for p,q in edges:
        a,b=interval_vectors(n,p,q);co=-2*(a.T*B0*b)[0]
        assert co!=0
        X[p,q]=1/co
    args=[em(n,i,(i+1)%n) for i in range(n)]+[s.diag(*D.row(i)) for i in range(n-2)]+[X]
    return target,args,order,tau,B0


def check_local_path(max_n):
    count=0
    for r in range(2,max_n+1):
        n=r+2;c=[1]+[2]*(r-2)+[1]
        edges=[(i,i+1) for i in range(r-1)]+[(i+1,i) for i in reversed(range(r-1))]
        endpoint=unit_standard(edges,n)
        assert endpoint==s.diag(*([(-1)**i*c[i] for i in range(r)]+[0,0]))
        for p in range(n):
            for q in range(n):
                M=unit_standard(edges+[(p,q)],n)
                if p==q:
                    exp=s.zeros(n)
                    if p<r:
                        for i in range(r):exp[i,i]=(-1)**p*c[i]*(c[p]+int(i==p))
                else:
                    co=0
                    if p<r and q<r and (p-q)%2==0:co=(-1)**p*c[p]*c[q]
                    elif p<r<=q:co=(-1)**p*c[p]
                    elif q<r<=p:co=(-1)**q*c[q]
                    exp=co*em(n,p,q)
                assert M==exp,(r,p,q,M,exp)
                count+=1
    return count


def check_cycle_table(max_n):
    count=0
    for n in range(3,max_n+1):
        cyc=[(i,(i+1)%n) for i in range(n)];A=2*s.ones(n)-s.eye(n)
        for i in range(n):
            for j in range(i+1,n):
                loops=[(l,l) for l in range(n) if l not in [i,j]]
                D=s.Matrix([[int(a==b) for b in range(n)] for a,_ in loops]);B=cycle_bivector(D)
                for p in [i,j]:
                    assert unit_standard(cyc+loops+[(p,p)],n)==s.diag(*list(A*B[:,p]))
                for p in range(n):
                    for q in range(n):
                        if p==q:continue
                        a,b=interval_vectors(n,p,q)
                        assert unit_standard(cyc+loops+[(p,q)],n)==-2*(a.T*B*b)[0]*em(n,p,q),(n,i,j,p,q)
                        count+=1
    return count


def check_periodic(max_k):
    records=[]
    for k in range(3,max_k+1):
        for t in range(1,k+1):
            n=k+t+1;w,g=periodic(k,t);L=slice_grade(w,g,2)
            L0=slice_grade(w,g,0);G,_=positive_formula(k+1,n)
            assert L0==G
            pos=L[:n-2,:n-2];neg=L[n-2:,n-2:]
            assert L[:n-2,n-2:]==s.zeros(n-2,2) and L[n-2:,:n-2]==s.zeros(2,n-2)
            for i in range(n-2):
                assert (-1)**i*pos[i,i]>=2
                for j in range(n-2):
                    if i!=j:assert (pos[i,j]==0 if (i-j)%k else abs(pos[i,j])==1)
            if t==k:assert neg==s.diag(-(2*k-1),-(2*k-1))
            elif t==1:
                assert list(neg.diagonal())==([3,-3] if k%2==0 else [-2,-2])
                assert abs(neg[0,1])<=1 and abs(neg[1,0])<=1
            else:assert neg==s.diag(*([4,-4] if n%2==0 else [-2,-2]))
            assert pos.det()!=0 and neg.det()!=0
            cross=None
            if k>=4:
                P=slice_grade(w,g,3)[:n-3,:n-3];inds=[0]+([k] if t>=3 else [])
                M=P.extract(inds,inds)
                if t<3:assert M==s.Matrix([[-2 if t==1 else -4]])
                else:
                    assert M[0,0]==-5 and M[1,1]==5*(-1)**k
                    assert abs(M[0,1])==abs(M[1,0])==1
                for i in inds:
                    for j in range(n-3):
                        if j not in inds:assert P[i,j]==P[j,i]==0
                assert M.det()!=0;cross=M.tolist()
            records.append({'k':k,'t':t,'grade_2_det':int(L.det()),'cross_parity_block':cross})
    return records


def main(output: Path,extended: bool):
    result={}
    # Cross-check two independent implementations on every input unit.
    w,g=periodic(2,2);args=[hom_matrix(a,b) for a,b in zip(w,g)]
    for h in range(5):
        L=slice_grade(w,g,h)
        for i in range(5):
            out=standard(args+[em(5,i,(i+h)%5)])
            expected=sum((L[j,i]*em(5,j,(j+h)%5) for j in range(5)),s.zeros(5))
            assert out==expected
    result['independent_full_matrix_comparisons']=25
    result['local_path_checks']=check_local_path(10 if extended else 8)
    result['cycle_coefficient_checks']=check_cycle_table(8 if extended else 6)
    examples=[]
    for n in range(3,9):
        cases=[[(0,n)],[(1,n)],[(1,n-1),(2,1)],[(1,n-1),(-(n-1),1)],
               [(1,1)]*n, [(1,2)]+[(1,1)]*(n-2),
               [(1,1)]*(n-1)+[(2,1)], [(0,2)]+[(0,1)]*(n-2)]
        if n>=4:cases.append([(0,n-2),(1,1),(-1,1)])
        if n>=6:cases.append([(0,n-4),(1,2),(-1,2)])
        for blocks in cases:
            T,A,order,tau,B=cycle_preimage(blocks)
            if n<=5:assert standard(A)==T
            examples.append({'blocks':blocks,'ordering':order,'avoidance_parameter':tau,'direct_full_evaluation':n<=5})
    result['minimal_size_Jordan_examples']=examples
    result['periodic_blocks']=check_periodic(7 if extended else 6)
    large=[]
    for r,n in [(2,2),(2,5),(3,3),(3,7),(4,4),(4,8),(5,5),(5,9)]:
        w,g=vandermonde_core(r,n);factor=s.prod(s.factorial(j) for j in range(1,r-1))**2
        L,b=positive_formula(r,n,factor)
        assert slice_grade(w,g,0)==L
        dp,_=subset_standard(w,g);assert s.Matrix(dp[-1])==b
        assert L.det()!=0 and all(z!=0 for z in b)
        large.append({'r':r,'n':n,'diagonal_det':int(L.det()),'endpoint_diagonal':list(map(int,b))})
    result['large_path_cores']=large
    tiny=[]
    for label,w,g,hs in [('s5_M5',*periodic(2,2),[0,2]),('s7_M6',*vandermonde_core(4,6),[0,2,3])]:
        mats=[slice_grade(w,g,h) for h in hs]
        assert all(M.det()!=0 for M in mats)
        aa=[hom_matrix(a,b) for a,b in zip(w,g)];nn=len(w[0])
        for h,M in zip(hs,mats):
            for i in range(nn):
                actual=standard(aa+[em(nn,i,(i+h)%nn)])
                expected=sum((M[j,i]*em(nn,j,(j+h)%nn) for j in range(nn)),s.zeros(nn))
                assert actual==expected
        tiny.append({'label':label,'weights':w,'grades':g,'restrictions':{str(h):M.tolist() for h,M in zip(hs,mats)},'determinants':[int(M.det()) for M in mats]})
    result['exceptional_certificates']=tiny
    # Check the 3-by-3 core used in the 4-by-4 transfer.
    core=[em(3,0,1),em(3,1,2),em(3,2,0),s.diag(-1,0,1)]
    alpha,beta,u,v,z=s.symbols('alpha beta u v z')
    X=s.diag((3*beta-4*alpha)/5,(2*beta-alpha)/5,0)+u*em(3,1,0)/2+v*em(3,2,1)/2-z*em(3,0,2)/4
    assert standard(core+[X])==s.Matrix([[alpha,0,z],[u,beta,0],[0,v,alpha]])
    assert standard(core)==s.diag(-1,2,-1)
    result['s5_M3_core_and_endpoint']='verified symbolically'
    output.write_text(json.dumps(result,indent=2,default=lambda x:int(x)),encoding='utf-8')
    print('All exact standard-polynomial checks passed.')
    print('Local path checks:',result['local_path_checks'])
    print('Cycle coefficient checks:',result['cycle_coefficient_checks'])
    print('Periodic core parameter pairs:',len(result['periodic_blocks']))
    for row in tiny:print(row['label'],row['determinants'])
    print('Certificate file:',output)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=Path('standard_certificates.json'))
    ap.add_argument('--extended',action='store_true')
    a=ap.parse_args();main(a.output,a.extended)
