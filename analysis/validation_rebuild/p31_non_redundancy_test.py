import sys; sys.path.insert(0,'.')
import numpy as np
from epc.models.cell_view_sorting import CellViewSorting
from epc.metrics.aggregation import MoransI, SegregationIndex, ClusterStats
def per_agent_dg(arrays):
    vals=arrays[0]; ideal={int(v):r for r,v in enumerate(np.sort(vals))}
    T=len(arrays); pos={int(v):np.empty(T,dtype=int) for v in vals}
    for t,a in enumerate(arrays):
        for i,v in enumerate(a): pos[int(v)][t]=i
    dg=[]
    for v in vals:
        p=pos[int(v)]; d=np.abs(p-ideal[int(v)]); mv=np.where(np.diff(p)!=0)[0]
        dg.append(float(np.sum(d[mv+1]>d[mv]))/len(mv) if len(mv) else 0.0)
    return np.array(dg)
def dgf(v): return [float(np.mean(v)),float(np.std(v)),float(np.quantile(v,.25)),float(np.quantile(v,.5)),float(np.quantile(v,.75))]
def rw(X,y,a=1.0):
    Xb=np.hstack([np.ones((len(X),1)),X]);A=Xb.T@Xb+a*np.eye(Xb.shape[1]);A[0,0]-=a;return np.linalg.solve(A,Xb.T@y)
def rp(X,w): return np.hstack([np.ones((len(X),1)),X])@w
def r2(y,yp): st=np.sum((y-y.mean())**2); return 1-np.sum((y-yp)**2)/st if st>0 else 0.0
mi,seg,cl=MoransI(),SegregationIndex(),ClusterStats(); algos=["bubble","insertion","selection"]; aidx={a:i for i,a in enumerate(algos)}
def build(n,nf,nseed):
    base,dgff,y=[],[],[]
    for algo in algos:
        for seed in range(nseed):
            m=CellViewSorting(n=n,algorithm=algo,n_frozen=nf,seed=seed)
            h=m.run_to_completion(max_rounds=3000); arr=[np.asarray(s["array"]) for s in h]
            oh=[0.0]*3; oh[aidx[algo]]=1.0
            mf=mi.compute(h,timestep=-1);sf=seg.compute(h,timestep=-1);cf=cl.compute(h,timestep=-1)
            base.append(oh+[mf["morans_i"],sf["segregation_index"],cf["cluster_count"],cf["mean_cluster_size"]])
            dgff.append(dgf(per_agent_dg(arr))); y.append(np.log(h[-1]["swap_count"]+1.0))
    return np.array(base),np.array(dgff),np.array(y)
def evaluate(Xb,Xd,y,reps=10):
    Xe=np.hstack([Xb,Xd]); deltas=[]; abl_gaps=[]
    for rep in range(reps):
        rng=np.random.default_rng(rep)
        Xds=Xd.copy()
        for c in range(Xds.shape[1]): rng.shuffle(Xds[:,c])
        Xa=np.hstack([Xb,Xds]); idx=rng.permutation(len(y)); fo=np.array_split(idx,5)
        db=[];de=[];da=[]
        for k in range(5):
            te=fo[k]; tr=np.concatenate([fo[j] for j in range(5) if j!=k])
            db.append(r2(y[te],rp(Xb[te],rw(Xb[tr],y[tr]))))
            de.append(r2(y[te],rp(Xe[te],rw(Xe[tr],y[tr]))))
            da.append(r2(y[te],rp(Xa[te],rw(Xa[tr],y[tr]))))
        deltas.append(np.mean(de)-np.mean(db)); abl_gaps.append(np.mean(de)-np.mean(da))
    return np.array(deltas),np.array(abl_gaps)
for (n,nf,ns) in [(40,3,120),(60,6,100)]:
    Xb,Xd,y=build(n,nf,ns)
    d,ag=evaluate(Xb,Xd,y,reps=10)
    print("config n=%d frozen=%d runs=%d: deltaR2=%.4f +/- %.4f | extended-ablation=%.4f +/- %.4f | survives=%s"%(
        n,nf,len(y),d.mean(),d.std(),ag.mean(),ag.std(), (d.mean()>0.01 and ag.mean()>0.01)))
