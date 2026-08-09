"""TERRAFORM WORLD -- niche-construction substrate (structural axis: DURABLE agent-shaped landscape).

Control stays grid-easy (move + dig). The NEW axis: agents durably MODIFY the world's PHYSICS (wall passability),
and the modification COMPOUNDS -- a dug shortcut speeds up every future goal-reaching, for the digger itself.

Distinct from stigmergy (competency #24): stigmergy marks are a COORDINATION SIGNAL other agents read; niche
construction changes the TRAVERSAL COST of the world and benefits even a SINGLE agent (no reader needed). The
discriminating channel is REVERT: dug walls restore each step (no persistence). If revert is load-bearing, the
behavior conditions on DURABLE SELF-MODIFICATION of the environment -- a content source no catalog behavior loads.

Grid with a central wall band pierceable only by digging. Goals respawn on the far side; agent must cross the wall
repeatedly. Digging a hole once lets all future crossings use it. Reward = ONLY goals reached; digging is never
rewarded, so infrastructure-building must be DISCOVERED (degenerate-novelty guard).

Ablation channels:
  blind    : zero other agents' positions
  revert   : agent's wall modifications do NOT persist (terrain resets to initial each step)  <- NOVEL channel
  noterrain: agent can't see terrain (walls invisible in obs)
  memwipe  : reset recurrent hidden each step
"""
import numpy as np, torch, torch.nn as nn

N=9; T=120; H=96; NACT=6   # 0=stay,1-4=move,5=dig(adjacent-forward cell)
DIRS=np.array([[0,0],[1,0],[-1,0],[0,1],[0,-1]])
WALL_ROW=N//2          # central wall band at row N//2
HARD=2                 # digs needed to breach a wall cell (hardness)

class TerraformWorld:
    def __init__(self, B, seed, n_ag=2, blind=False, revert=False, noterrain=False, memwipe=False):
        self.B=B; self.n=n_ag; self.blind=blind; self.revert=revert
        self.noterrain=noterrain; self.memwipe=memwipe
        self.rng=np.random.default_rng(seed); self.reset()
    @property
    def odim(self): return 2+1 + N + 2 + (self.n-1)*2  # own pos(2)+facing(1) | wall-row passability(N) | goal rel(2) | others rel
    def reset(self):
        B,n=self.B,self.n
        # wall: row WALL_ROW, each cell has 'health' HARD (impassable until dug to 0)
        self.wall=np.full((B,N),HARD,np.int16)   # health per column of the wall row
        self.wall0=self.wall.copy()
        # agents start on near side (rows < WALL_ROW)
        self.ap=np.stack([self.rng.integers(0,WALL_ROW,size=(B,n)),self.rng.integers(0,N,size=(B,n))],axis=2).astype(int)
        self.face=self.rng.integers(1,5,size=(B,n))   # facing dir for dig
        # goal starts on FAR side; alternates sides each time it's reached (forces REPEATED crossing => hole REUSE)
        self.goal=np.stack([self.rng.integers(WALL_ROW+1,N,size=B),self.rng.integers(0,N,size=B)],axis=1)
        self.goal_far=np.ones(B,bool)   # True => goal currently on far side
        self._t=0; self._goals=np.zeros(B); return self.obs()
    def _passable(self, r, c, b):
        # a cell is passable unless it's the wall row with health>0
        wall_here = (r==WALL_ROW)
        return ~(wall_here & (self.wall[b,c]>0))
    def obs(self):
        B,n=self.B,self.n; out=np.zeros((B,n,self.odim),np.float32)
        wallview = (np.zeros_like(self.wall) if self.noterrain else self.wall)
        for i in range(n):
            o=0
            out[:,i,o:o+2]=self.ap[:,i]/(N-1); o+=2
            out[:,i,o]=self.face[:,i]/4.0; o+=1
            out[:,i,o:o+N]=wallview/HARD; o+=N
            out[:,i,o:o+2]=(self.goal-self.ap[:,i])/(N-1); o+=2
            for j in range(n):
                if j==i: continue
                rel=np.zeros((B,2),np.float32) if self.blind else (self.ap[:,j]-self.ap[:,i])/(N-1)
                out[:,i,o:o+2]=rel; o+=2
        return out
    def step(self, act):
        B,n=self.B,self.n; self._t+=1
        for i in range(n):
            a=act[:,i]
            mv=(a>=1)&(a<=4)
            self.face[:,i]=np.where(mv,a,self.face[:,i])
            # attempt move
            d=DIRS[np.where(mv,a,0)]
            nr=np.clip(self.ap[:,i,0]+d[:,0],0,N-1); nc=np.clip(self.ap[:,i,1]+d[:,1],0,N-1)
            ok=self._passable(nr,nc,np.arange(B))
            self.ap[:,i,0]=np.where(mv&ok,nr,self.ap[:,i,0]); self.ap[:,i,1]=np.where(mv&ok,nc,self.ap[:,i,1])
            # dig: reduce health of the cell the agent faces, if it's the wall row
            dig=(a==5)
            fd=DIRS[self.face[:,i]]
            tr=np.clip(self.ap[:,i,0]+fd[:,0],0,N-1); tc=np.clip(self.ap[:,i,1]+fd[:,1],0,N-1)
            hit=dig&(tr==WALL_ROW)&(self.wall[np.arange(B),tc]>0)
            self.wall[np.arange(B),tc]=np.where(hit,self.wall[np.arange(B),tc]-1,self.wall[np.arange(B),tc])
        if self.revert: self.wall=self.wall0.copy()   # modifications don't persist
        # reward: reached goal? (any agent on goal cell). respawn goal on far side.
        rew=np.zeros(B,np.float32)
        for i in range(n):
            on=np.all(self.ap[:,i]==self.goal,axis=1)
            rew+=on.astype(np.float32)
            if on.any():
                self._goals+=on
                # flip side for reached worlds: far->near, near->far (forces recrossing the wall)
                self.goal_far[on]=~self.goal_far[on]
                far=self.goal_far
                ng_r=np.where(far,self.rng.integers(WALL_ROW+1,N,size=B),self.rng.integers(0,WALL_ROW,size=B))
                ng_c=self.rng.integers(0,N,size=B)
                self.goal[on]=np.stack([ng_r,ng_c],axis=1)[on]
        return np.repeat(rew[:,None],n,axis=1)
    def score(self):  # mean goals reached per episode (per world), normalized by a soft cap
        return float(self._goals.mean())

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
            O[:,i,t]=obs[:,i]; A[:,i,t]=a.numpy(); LP[:,i,t]=d.log_prob(a).numpy(); V[:,i,t]=v[:,0].numpy()
            act[:,i]=a.numpy()
        r=env.step(act); R[:,:,t]=r; obs=env.obs()
    return O,A,LP,V,R

def gae(R,V,gamma=0.99,lam=0.95):
    B,Tt=R.shape; adv=np.zeros((B,Tt),np.float32); last=np.zeros(B,np.float32)
    for t in reversed(range(Tt)):
        nv=V[:,t+1] if t+1<Tt else np.zeros(B,np.float32)
        delta=R[:,t]+gamma*nv-V[:,t]; last=delta+gamma*lam*last; adv[:,t]=last
    return adv, adv+V

def train(iters, seed=0, n_ag=2, B=128):
    torch.manual_seed(seed); od=TerraformWorld(1,0,n_ag=n_ag).odim
    net=Pol(od); opt=torch.optim.Adam(net.parameters(),lr=3e-3)
    for it in range(iters):
        env=TerraformWorld(B,1000+it,n_ag=n_ag)
        O,A,LP,V,R=rollout(net,env)
        od_=od
        O2=O.reshape(B*n_ag,T,od_); A2=A.reshape(B*n_ag,T); LP2=LP.reshape(B*n_ag,T)
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
            ev=TerraformWorld(400,7000+it,n_ag=n_ag); rollout(net,ev,greedy=True)
            print(f"  iter {it:>3}: goals/episode {ev.score():.2f}",flush=True)
    return net

def _score(net, seed, B=400, **kw):
    env=TerraformWorld(B,seed,**kw); rollout(net,env,greedy=True); return env.score()
