"""HETERO WORLD -- heterogeneous-bodies substrate (structural axis: OWN-MORPHOLOGY conditioning).

Control grid-easy. The NEW axis: agents have different BODIES (affordances), and must learn to condition behavior
on their OWN body type -- self-sorting by morphology, not by an identity label.

Distinct from division-of-labor (#22, identity-anchored, IDENTICAL agents self-assigning roles) and from noid
(zeros an identity ONE-HOT label): here type-0 harvests RED food efficiently (2x), type-1 harvests BLUE. Both see
both foods. The efficient strategy is each type specializing to ITS food. The discriminating channel is NOBODY:
zero the agent's own body-type observation. If load-bearing, the policy conditions on its own morphology to decide
what to harvest. Reward = total food value collected (shared); specialization is DISCOVERED not rewarded
(both foods give the same base reward; only the harvest RATE differs by body).

Ablation channels:
  blind    : zero other agents' positions
  nobody   : zero own body-type one-hot (agent doesn't know which body it is)   <- NOVEL channel
  noid     : zero own identity index (control: identity label, distinct from body affordance)
  memwipe  : reset recurrent hidden each step
"""
import numpy as np, torch, torch.nn as nn

N=9; T=60; H=96; NACT=5   # 0=stay,1-4=move
DIRS=np.array([[0,0],[1,0],[-1,0],[0,1],[0,-1]])
NRED=4; NBLUE=4; FAST=2.0; SLOW=0.5   # harvest-rate multipliers by matching body

class HeteroWorld:
    def __init__(self, B, seed, n_ag=4, blind=False, nobody=False, noid=False, memwipe=False):
        self.B=B; self.n=n_ag; self.blind=blind; self.nobody=nobody; self.noid=noid; self.memwipe=memwipe
        self.rng=np.random.default_rng(seed); self.reset()
    @property
    def odim(self):
        # own pos(2)+body one-hot(2)+id index(1) | red rel(NRED*2 nearest 3)| blue rel | others rel
        return 2+2+1 + 3*2 + 3*2 + (self.n-1)*2
    def reset(self):
        B,n=self.B,self.n
        self.ap=self.rng.integers(0,N,size=(B,n,2))
        self.body=np.tile(np.arange(n)%2,(B,1))   # alternate red/blue bodies
        self.red=self.rng.integers(0,N,size=(B,NRED,2)); self.rq=np.ones((B,NRED))
        self.blue=self.rng.integers(0,N,size=(B,NBLUE,2)); self.bq=np.ones((B,NBLUE))
        self._val=np.zeros(B); return self.obs()
    def _nearest3(self, items, qty, ap):
        B=self.B; d=items-ap[:,None,:]; dist=np.linalg.norm(d,axis=2)+(qty<=0)*999
        idx=np.argsort(dist,axis=1)[:,:3]
        rel=np.take_along_axis(items,idx[:,:,None],axis=1)-ap[:,None,:]
        return (rel/(N-1)).reshape(B,-1)
    def obs(self):
        B,n=self.B,self.n; out=np.zeros((B,n,self.odim),np.float32)
        for i in range(n):
            o=0
            out[:,i,o:o+2]=self.ap[:,i]/(N-1); o+=2
            if not self.nobody:
                out[np.arange(B),i,o+self.body[:,i]]=1.0
            o+=2
            out[:,i,o]=0.0 if self.noid else i/max(n-1,1); o+=1
            out[:,i,o:o+6]=self._nearest3(self.red,self.rq,self.ap[:,i]); o+=6
            out[:,i,o:o+6]=self._nearest3(self.blue,self.bq,self.ap[:,i]); o+=6
            for j in range(n):
                if j==i: continue
                rel=np.zeros((B,2),np.float32) if self.blind else (self.ap[:,j]-self.ap[:,i])/(N-1)
                out[:,i,o:o+2]=rel; o+=2
        return out
    def step(self, act):
        B,n=self.B,self.n
        for i in range(n):
            d=DIRS[act[:,i]]
            self.ap[:,i,0]=np.clip(self.ap[:,i,0]+d[:,0],0,N-1); self.ap[:,i,1]=np.clip(self.ap[:,i,1]+d[:,1],0,N-1)
        rew=np.zeros(B,np.float32)
        for i in range(n):
            body=self.body[:,i]
            # harvest red
            for f in range(NRED):
                on=np.all(self.ap[:,i]==self.red[:,f],axis=1)&(self.rq[:,f]>0)
                rate=np.where(body==0,FAST,SLOW)   # body-0 fast on red
                got=on*rate*self.rq[:,f]; self.rq[:,f]=np.clip(self.rq[:,f]-on*0.5,0,1); rew+=got*0.0
                # value credited = base 1.0 per harvest event, but harvest SPEED differs -> track collected value
                rew+=on.astype(np.float32)*np.where(body==0,1.0,0.0)
            for f in range(NBLUE):
                on=np.all(self.ap[:,i]==self.blue[:,f],axis=1)&(self.bq[:,f]>0)
                self.bq[:,f]=np.clip(self.bq[:,f]-on*0.5,0,1)
                rew+=on.astype(np.float32)*np.where(body==1,1.0,0.0)
        # regrow slowly
        self.rq=np.clip(self.rq+0.01,0,1); self.bq=np.clip(self.bq+0.01,0,1)
        self._val+=rew
        return np.repeat(rew[:,None],n,axis=1)
    def score(self): return float(self._val.mean())

class Pol(nn.Module):
    def __init__(self, odim):
        super().__init__(); self.recurrent=True
        self.gru=nn.GRU(odim,H,batch_first=True); self.pi=nn.Linear(H,NACT); self.v=nn.Linear(H,1)
    def forward(self,x,h=None):
        z,h=self.gru(x,h); return self.pi(z),self.v(z).squeeze(-1),h

def rollout(net, env, greedy=False):
    B,n=env.B,env.n; obs=env.obs(); h=[None]*n; memwipe=getattr(env,"memwipe",False)
    O=np.zeros((B,n,T,env.odim),np.float32); A=np.zeros((B,n,T),int)
    LP=np.zeros((B,n,T),np.float32); V=np.zeros((B,n,T),np.float32); R=np.zeros((B,n,T),np.float32)
    for t in range(T):
        act=np.zeros((B,n),int)
        for i in range(n):
            hin=None if memwipe else h[i]
            with torch.no_grad(): lg,v,h[i]=net(torch.from_numpy(obs[:,i])[:,None,:],hin)
            d=torch.distributions.Categorical(logits=lg[:,0]); a=lg[:,0].argmax(1) if greedy else d.sample()
            O[:,i,t]=obs[:,i]; A[:,i,t]=a.numpy(); LP[:,i,t]=d.log_prob(a).numpy(); V[:,i,t]=v[:,0].numpy(); act[:,i]=a.numpy()
        r=env.step(act); R[:,:,t]=r; obs=env.obs()
    return O,A,LP,V,R

def gae(R,V,gamma=0.98,lam=0.95):
    B,Tt=R.shape; adv=np.zeros((B,Tt),np.float32); last=np.zeros(B,np.float32)
    for t in reversed(range(Tt)):
        nv=V[:,t+1] if t+1<Tt else np.zeros(B,np.float32)
        delta=R[:,t]+gamma*nv-V[:,t]; last=delta+gamma*lam*last; adv[:,t]=last
    return adv, adv+V

def train(iters, seed=0, n_ag=4, B=128):
    torch.manual_seed(seed); od=HeteroWorld(1,0,n_ag=n_ag).odim
    net=Pol(od); opt=torch.optim.Adam(net.parameters(),lr=3e-3)
    for it in range(iters):
        env=HeteroWorld(B,1000+it,n_ag=n_ag); O,A,LP,V,R=rollout(net,env)
        O2=O.reshape(B*n_ag,T,od); A2=A.reshape(B*n_ag,T); LP2=LP.reshape(B*n_ag,T)
        V2=V.reshape(B*n_ag,T); R2=R.reshape(B*n_ag,T)
        adv,ret=gae(R2,V2); adv=(adv-adv.mean())/(adv.std()+1e-8)
        Ot,At,LPt=torch.from_numpy(O2),torch.from_numpy(A2),torch.from_numpy(LP2)
        advt,rett=torch.from_numpy(adv),torch.from_numpy(ret)
        for _ in range(4):
            lg,v,_=net(Ot); d=torch.distributions.Categorical(logits=lg)
            ratio=torch.exp(d.log_prob(At)-LPt); s1=ratio*advt; s2=torch.clamp(ratio,0.8,1.2)*advt
            loss=-torch.min(s1,s2).mean()+0.5*((v-rett)**2).mean()-0.02*d.entropy().mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(),0.5); opt.step()
        if it%100==0 or it==iters-1:
            ev=HeteroWorld(400,7000+it,n_ag=n_ag); rollout(net,ev,greedy=True)
            print(f"  iter {it:>3}: value/episode {ev.score():.2f}",flush=True)
    return net

def _score(net, seed, B=400, **kw):
    env=HeteroWorld(B,seed,**kw); rollout(net,env,greedy=True); return env.score()
