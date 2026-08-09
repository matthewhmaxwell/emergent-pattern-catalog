"""ADAPTHETERO -- compositional substrate: non-stationary opponent (social) x heterogeneous body (morphology).
Two INFORMATION-ORTHOGONAL axes, each with an optimal action that VARIES BY STATE (the necessity condition
metabolism failed): opponent-block-switch (adapt #32) AND own-body-type (hetero #34).
A single red-or-blue agent on a grid with a decaying energy budget. Two food sources exist for EACH body type
(red-food-A/red-food-B, blue-food-A/blue-food-B), i.e. 2 spatial sites, each site serves BOTH colors but only the
color-matching agent can eat there. An opponent BLOCKS one of the two sites and SWITCHES which every PERIOD steps;
blocked-site food is unavailable to everyone. To survive the agent must (a) know its OWN BODY (only its color's food
nourishes it) AND (b) track WHICH site is currently open. Optimal target = the OPEN site, eat only if it matches body.
Neither leaks: body says nothing about which site is blocked; block-state says nothing about body.
Channels: blind, nobody(zero own body one-hot), freezeblock(freeze block-state to initial), memwipe."""
import numpy as np
G=7; T=80; NACT=5
E0=4.0; METAB=0.34; INTAKE=2.2; EMAX=5.0
PERIOD=20; CPERIOD=13
SITES=np.array([[1,3],[5,3]])          # two spatial food sites
MOVES=np.array([[0,0],[0,1],[0,-1],[1,0],[-1,0]])

class AdaptHetero:
    def __init__(self,B,seed,n_ag=1,blind=False,nobody=False,freezeblock=False,memwipe=False):
        self.B=B; self.rng=np.random.default_rng(seed)
        self.blind=blind; self.nobody=nobody; self.freezeblock=freezeblock; self.memwipe=memwipe
        self.reset()
    def reset(self):
        B=self.B
        self.ap=self.rng.integers(0,G,size=(B,2))
        self.e=np.full(B,E0,np.float32)
        self.body=self.rng.integers(0,2,size=B)                        # 0 red, 1 blue
        self.block0=self.rng.integers(0,2,size=B)                      # initially blocked SITE
        self.col0=self.rng.integers(0,2,size=B)                        # color served by site0 in phase 0
        self.alive=np.ones(B,bool); self.t=0; self._surv=np.zeros(B,np.float32)
        o=self.obs(); self.odim=o.shape[-1]; return o
    def _cur_block(self):
        phase=(self.t//PERIOD)%2
        b=np.where(phase==0,self.block0,1-self.block0)
        if self.freezeblock: b=self.block0
        return b
    def _site_color(self,si):
        # color flips on its OWN clock (CPERIOD), independent of the block clock (PERIOD)
        phase=(self.t//CPERIOD)%2
        base=self.col0 if si==0 else 1-self.col0
        return np.where(phase==0,base,1-base)
    def obs(self):
        B=self.B
        pos=self.ap/G
        rel=np.concatenate([(SITES[0][None]-self.ap)/G,(SITES[1][None]-self.ap)/G],axis=1).astype(np.float32)
        # per-site served color (dynamic, observable): 2 sites x nothing -> encode each site color as scalar
        c0=self._site_color(0)[:,None].astype(np.float32); c1=self._site_color(1)[:,None].astype(np.float32)
        sitecol=np.concatenate([c0,c1],axis=1)
        if self.blind: pos=np.zeros_like(pos); rel=np.zeros_like(rel); sitecol=np.zeros_like(sitecol)
        blk=self._cur_block(); blk1h=np.zeros((B,2),np.float32); blk1h[np.arange(B),blk]=1.0
        if self.blind: blk1h=np.zeros_like(blk1h)
        body1h=np.zeros((B,2),np.float32); body1h[np.arange(B),self.body]=1.0
        if self.nobody: body1h=np.zeros_like(body1h)
        en=(self.e/EMAX)[:,None].astype(np.float32)
        clock=np.full((B,1),self.t/T,np.float32)
        return np.concatenate([pos,rel,sitecol,blk1h,body1h,en,clock],axis=1)[:,None,:].astype(np.float32)
    def step(self,act):
        B=self.B; a=act[:,0]; self.ap=np.clip(self.ap+MOVES[a],0,G-1); self.e-=METAB
        blk=self._cur_block(); openidx=1-blk
        for si in range(2):
            on=(self.ap[:,0]==SITES[si][0])&(self.ap[:,1]==SITES[si][1])
            canfeed=on&(openidx==si)&(self._site_color(si)==self.body)&self.alive
            self.e=np.where(canfeed,np.minimum(self.e+INTAKE,EMAX),self.e)
        self.alive&=(self.e>0); self._surv+=self.alive.astype(np.float32); self.t+=1
        return np.zeros((B,1),np.float32)
    def score(self): return float(self._surv.mean()/T)

def rollout(net,env):
    import torch
    B=env.B; obs=env.obs()
    O=np.zeros((B,1,T,env.odim),np.float32);A=np.zeros((B,1,T),np.int64)
    LP=np.zeros((B,1,T),np.float32);V=np.zeros((B,1,T),np.float32);R=np.zeros((B,1,T),np.float32)
    h=[None]; prev=np.zeros(B,np.float32)
    for t in range(T):
        if env.memwipe: h[0]=None
        ot=torch.from_numpy(obs[:,0])[:,None,:]
        with torch.no_grad(): lg,v,h[0]=net(ot,h[0])
        d=torch.distributions.Categorical(logits=lg[:,0]); aa=d.sample()
        O[:,0,t]=obs[:,0];A[:,0,t]=aa.numpy();LP[:,0,t]=d.log_prob(aa).numpy();V[:,0,t]=v[:,0].numpy()
        env.step(A[:,:,t]); cur=env._surv.copy(); R[:,0,t]=cur-prev; prev=cur
        obs=env.obs()
    return O,A,LP,V,R
def gae(R,V,gamma=0.99,lam=0.95):
    BN,T=R.shape; adv=np.zeros_like(R); last=np.zeros(BN,np.float32)
    for t in reversed(range(T)):
        nv=V[:,t+1] if t+1<T else np.zeros(BN,np.float32)
        delta=R[:,t]+gamma*nv-V[:,t]; last=delta+gamma*lam*last; adv[:,t]=last
    return adv, adv+V
def _score(net,seed,B=400,**kw):
    import torch
    env=AdaptHetero(B,seed,n_ag=1,**kw); obs=env.obs(); h=[None]
    for t in range(T):
        if env.memwipe: h[0]=None
        with torch.no_grad(): lg,_,h[0]=net(torch.from_numpy(obs[:,0])[:,None,:],h[0])
        mv=lg[:,0].argmax(1).numpy()[:,None]; env.step(mv); obs=env.obs()
    return env.score()
