import numpy as np, json, time, sys
rng=np.random.default_rng(11)
def irls(X,y,it=8):
    B,n,p=X.shape; b=np.zeros((B,p))
    for _ in range(it):
        eta=np.einsum('bnp,bp->bn',X,b); mu=1/(1+np.exp(-eta)); w=mu*(1-mu)+1e-9
        H=np.einsum('bnp,bn,bnq->bpq',X,w,X)+1e-4*np.eye(p)
        g=np.einsum('bnp,bn->bp',X,y-mu)
        b=b+np.linalg.solve(H,g[...,None])[...,0]
    return b
def design(sm,si,g,inter):
    cols=[np.ones_like(sm),sm,si,g]+([sm*si] if inter else [])
    return np.stack(cols,-1)
def dev(X,y,b):
    mu=np.clip(1/(1+np.exp(-np.einsum('bnp,bp->bn',X,b))),1e-6,1-1e-6)
    return -2*(y*np.log(mu)+(1-y)*np.log(1-mu)).sum(-1)
def cig(sm,si,g,y,folds):
    tot=np.zeros(y.shape[0])
    for f in range(5):
        tr=folds!=f; te=~tr
        for inter,sgn in ((False,1),(True,-1)):
            b=irls(design(sm[:,tr],si[:,tr],g[:,tr],inter),y[:,tr])
            tot+=sgn*dev(design(sm[:,te],si[:,te],g[:,te],inter),y[:,te],b)
    return tot
def one(n,gamma,pa,pb,P=199):
    sm=(rng.random(n)<0.4).astype(float)
    si=(rng.random(n)<np.where(sm==1,pa,pb)).astype(float)
    g=(rng.random(n)<0.5).astype(float)
    eta=-1.6+0.7*sm+0.7*si+gamma*sm*si+0.3*g
    y=(rng.random(n)<1/(1+np.exp(-eta))).astype(float)
    folds=rng.permutation(np.arange(n)%5)
    T=cig(sm[None],si[None],g[None],y[None],folds)[0]
    # parametric bootstrap under fitted additive model, predictors fixed
    Xa=design(sm,si,g,False)[None]; ba=irls(Xa,y[None])[0]
    pa_hat=1/(1+np.exp(-(Xa[0]@ba)))
    ys=(rng.random((P,n))<pa_hat).astype(float)
    Tb=cig(np.tile(sm,(P,1)),np.tile(si,(P,1)),np.tile(g,(P,1)),ys,folds)
    p_boot=(1+(Tb>=T).sum())/(P+1)
    # permutation within (y,g)
    perm=np.tile(np.arange(n),(P,1))
    for yy in (0,1):
        for gg in (0,1):
            idx=np.where((y==yy)&(g==gg))[0]
            for b in range(P): perm[b,idx]=rng.permutation(idx)
    Tp=cig(np.tile(sm,(P,1)),si[perm],np.tile(g,(P,1)),np.tile(y,(P,1)),folds)
    p_perm=(1+(Tp>=T).sum())/(P+1)
    bi=irls(design(sm,si,g,True)[None],y[None])
    lrt=dev(Xa,y[None],ba[None])[0]-dev(design(sm,si,g,True)[None],y[None],bi)[0]
    from math import erfc,sqrt
    return p_boot,p_perm,erfc(sqrt(max(lrt,0)/2)),y.mean()
R=int(sys.argv[1]) if len(sys.argv)>1 else 150
res={}; t=time.time()
scen=[('moderate',0.6,0.3),('strong',0.85,0.1)]
for name,pa,pb in scen:
    gl=(0.0,0.5,1.0,1.5,2.0) if name=='moderate' else (0.0,)
    for n in (200,400):
        for gm in gl:
            a=np.array([one(n,gm,pa,pb) for _ in range(R)])
            k=f"{name}_{n}_{gm}"
            res[k]=dict(boot=float((a[:,0]<0.05).mean()),perm=float((a[:,1]<0.05).mean()),lrt=float((a[:,2]<0.05).mean()),pcr=float(a[:,3].mean()),R=R)
            print(k,res[k],round(time.time()-t),flush=True)
json.dump(res,open('/tmp/claude-0/-home-claude/b3f81ece-5cac-5452-b51c-d5e3b06ec265/scratchpad/sim/res2.json','w'),indent=1)
