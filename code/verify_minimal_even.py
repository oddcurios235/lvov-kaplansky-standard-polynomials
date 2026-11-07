#!/usr/bin/env python3
"""Exact finite tests for the uniform minimal even-standard seed.
Matrix units are zero-indexed pairs. A signed state recurrence and explicit
Euler-trail enumeration check the values; full differentials are tested in
specified finite matrix orders.
"""
from __future__ import annotations
import argparse, json, math, time
from functools import lru_cache
from pathlib import Path
import sympy as sp


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def std_units(edges: list[tuple[int,int]], n: int) -> list[list[int]]:
    m=len(edges); full=(1<<m)-1
    outgoing=[[] for _ in range(n)]
    for i,(a,b) in enumerate(edges): outgoing[a].append((i,b))
    @lru_cache(None)
    def walk(mask: int, vertex: int) -> tuple[int,...]:
        if mask==full: return tuple(int(j==vertex) for j in range(n))
        result=[0]*n
        for i,b in outgoing[vertex]:
            if (mask>>i)&1: continue
            sign=-1 if (mask>>(i+1)).bit_count()%2 else 1
            tail=walk(mask|(1<<i),b)
            for j,v in enumerate(tail): result[j]+=sign*v
        return tuple(result)
    return [list(walk(0,a)) for a in range(n)]


def std_trails(edges: list[tuple[int,int]], n:int) -> list[list[int]]:
    """Separate enumeration, permutation sign computed only on completed trails."""
    result=[[0]*n for _ in range(n)]; m=len(edges)
    def rec(start,current,word,remaining):
        if not remaining:
            inv=sum(word[i]>word[j] for i in range(m) for j in range(i+1,m))
            result[start][current]+=(-1)**inv
            return
        for k in remaining:
            if edges[k][0]==current:
                rec(start,edges[k][1],word+[k],remaining-{k})
    for start in range(n): rec(start,start,[],set(range(m)))
    return result


def seed(n:int):
    return [(0,0),(0,1)]+[e for j in range(2,n) for e in [(0,j),(j,0)]]


def unit(n,i,j,c=1):
    a=[[0]*n for _ in range(n)];a[i][j]=c;return a


def check(n:int, full_rank:bool):
    e=seed(n); p=n-2; c=math.factorial(p+1); b=math.factorial(p)
    trials=[('value',None,None,unit(n,0,1,c)),
            ('exit_reverse',1,(1,0),unit(n,1,0,-c))]
    if n>=3:
        diag=unit(n,1,1,b);diag[2][2]=-b
        trials += [('exit_to_leaf',1,(1,2),unit(n,1,2,-b)),
                   ('outgoing_reverse',2,(1,0),unit(n,2,0,-c)),
                   ('outgoing_triangle',2,(1,2),diag)]
    values={}
    for name,slot,replacement,expected in trials:
        ed=e.copy()
        if slot is not None: ed[slot]=replacement
        out=std_units(ed,n)
        require(out==expected,f'coefficient mismatch {n} {name}: {out} != {expected}')
        if n<=6:
            require(std_trails(ed,n)==out,f'independent trails {n} {name}')
        values[name]=out
    res={'n':n,'degree':2*n-2,'seed':e,'selected_values':values,'independent_trails':n<=6}
    if full_rank:
        cols=[]
        for slot in range(len(e)):
            for i in range(n):
                for j in range(n):
                    ed=e.copy();ed[slot]=(i,j)
                    out=std_units(ed,n)
                    require(sum(out[k][k] for k in range(n))==0,'trace derivative')
                    cols.append([v for row in out for v in row])
        J=sp.Matrix.hstack(*map(sp.Matrix,cols))
        piv=list(J.rref()[1]);rank=len(piv)
        require(rank==n*n-1,f'full differential rank {n}: {rank}')
        sub=J[:,piv];rows=list(sub.T.rref()[1]);minor=int(sub.extract(rows,range(rank)).det())
        require(minor!=0,'zero pivot minor')
        res.update(jacobian_rank=rank,minor_rows=rows,minor_columns=piv,minor=str(minor))
    return res


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--max-n',type=int,default=9)
    parser.add_argument('--rank-max-n',type=int,default=6)
    parser.add_argument('--output',default='minimal_even.json');args=parser.parse_args()
    t=time.monotonic();record={'status':'running','cases':[]}
    try:
        for n in range(2,args.max_n+1):
            a=check(n,n<=args.rank_max_n);record['cases'].append(a)
            print(n,'coefficient identities checked','rank='+str(a.get('jacobian_rank','not computed')),flush=True)
        record['status']='passed'
    except Exception as exc:
        record.update(status='failed',error=repr(exc));raise
    finally:
        record['seconds']=time.monotonic()-t
        Path(args.output).parent.mkdir(parents=True,exist_ok=True)
        Path(args.output).write_text(json.dumps(record,indent=2))
    print('PASSED',record['seconds'],flush=True)
if __name__=='__main__':main()
