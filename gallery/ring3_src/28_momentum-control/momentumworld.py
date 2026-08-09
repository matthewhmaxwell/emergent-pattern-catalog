"""MOMENTUM WORLD -- continuous-physics multi-agent substrate (structural axis: INERTIA/dynamics).

First substrate to leave the discrete grid. Agents are point masses with velocity + momentum in a 2D box.
Actions apply a bounded FORCE (thrust), not a teleport. Velocity persists (inertia), with light drag.

Task: N agents must together bring a heavy PUCK (also a point mass, heavier) to a goal region. A single agent's
max thrust is too weak to move the puck efficiently against drag -- but agents can BUILD UP MOMENTUM (accelerate
over several steps) and transfer it to the puck on contact (elastic-ish collision). So the *possible* emergent
mechanism is momentum management: wind up, strike, recover -- a behavior the grid family structurally cannot express
(no velocity state). Crucially the reward is ONLY puck-in-goal proximity; nothing rewards momentum, striking, or
coordination directly -> a physics-exploiting mechanism must be DISCOVERED, not built in (degenerate-novelty guard).

Ablation channels (for the battery, content-source axis):
  blind        : zero other agents' relative state (positions+velocities)
  freezeother  : freeze others' state to t=0 (removes their DYNAMICS, keeps presence)
  novel        : zero own VELOCITY in the observation (agent can't perceive its own momentum) -> if load-bearing,
                 the behavior conditions on OWN-DYNAMICS, a content source no catalog behavior loads
  freezeenv    : freeze the puck+goal state to t=0 (environment dynamics)
  memwipe      : reset recurrent hidden each step
"""
import numpy as np, torch, torch.nn as nn

DT=0.2; DRAG=0.06; FMAX=2.0; VMAX=4.0; T=50
NMOVE=5   # 0=none,1=+x,2=-x,3=+y,4=-y thrust directions
DIRS=np.array([[0,0],[1,0],[-1,0],[0,1],[0,-1]],np.float32)
BOX=10.0; PUCK_MASS=2.0; GOAL_R=2.0; H=96; CONTACT=1.2

def _clip_v(v): 
    sp=np.linalg.norm(v,axis=-1,keepdims=True); f=np.minimum(1.0,VMAX/np.maximum(sp,1e-6)); return v*f

class MomentumWorld:
    def __init__(self, B, seed, n_ag=2, blind=False, freezeother=False, novel=False,
                 freezeenv=False, memwipe=False, spawn_near=False):
        self.B=B; self.n=n_ag; self.blind=blind; self.freezeother=freezeother
        self.novel=novel; self.freezeenv=freezeenv; self.memwipe=memwipe
        self.spawn_near=spawn_near
        self.rng=np.random.default_rng(seed); self.reset()
    # obs per agent: own pos(2)+own vel(2) | puck rel(2)+puck vel(2)+goal rel(2) | others rel(2)+vel(2) each
    @property
    def odim(self): return 2+2 + 2+2+2 + (self.n-1)*4
    def reset(self):
        B,n=self.B,self.n
        self.pp=self.rng.uniform(3,BOX-3,size=(B,2)).astype(np.float32)
        self.pv=np.zeros((B,2),np.float32)
        self.goal=self.rng.uniform(1,BOX-1,size=(B,2)).astype(np.float32)
        if getattr(self,'spawn_near',False):
            self.ap=(self.pp[:,None,:]+self.rng.uniform(-2,2,size=(B,n,2))).clip(0,BOX).astype(np.float32)
        else:
            self.ap=self.rng.uniform(1,BOX-1,size=(B,n,2)).astype(np.float32)
        self.av=np.zeros((B,n,2),np.float32)
        self.ap0=self.ap.copy(); self.av0=self.av.copy(); self.pp0=self.pp.copy(); self.pv0=self.pv.copy()
        self._t=0; return self.obs()
    def obs(self):
        B,n=self.B,self.n; out=np.zeros((B,n,self.odim),np.float32)
        pp = self.pp0 if self.freezeenv else self.pp
        pv = self.pv0 if self.freezeenv else self.pv
        for i in range(n):
            o=0
            out[:,i,o:o+2]=self.ap[:,i]/BOX; o+=2
            out[:,i,o:o+2]=(0*self.av[:,i] if self.novel else self.av[:,i])/VMAX; o+=2  # own velocity (novel zeros it)
            out[:,i,o:o+2]=(pp-self.ap[:,i])/BOX; o+=2
            out[:,i,o:o+2]=pv/VMAX; o+=2
            out[:,i,o:o+2]=(self.goal-self.ap[:,i])/BOX; o+=2
            for j in range(n):
                if j==i: continue
                if self.blind:
                    rel=np.zeros((B,2),np.float32); vel=np.zeros((B,2),np.float32)
                else:
                    apj=self.ap0[:,j] if self.freezeother else self.ap[:,j]
                    avj=self.av0[:,j] if self.freezeother else self.av[:,j]
                    rel=(apj-self.ap[:,i])/BOX; vel=avj/VMAX
                out[:,i,o:o+2]=rel; o+=2; out[:,i,o:o+2]=vel; o+=2
        return out
    def step(self, mv):
        B,n=self.B,self.n; self._t+=1
        F=DIRS[mv]*FMAX                         # (B,n,2) thrust
        self.av=_clip_v((self.av+F*DT)*(1-DRAG))
        self.ap=np.clip(self.ap+self.av*DT,0,BOX)
        # agent-puck elastic-ish momentum transfer on contact
        for i in range(n):
            d=self.pp-self.ap[:,i]; dist=np.linalg.norm(d,axis=1)
            hit=dist<CONTACT
            if hit.any():
                nrm=d/np.maximum(dist[:,None],1e-6)
                # transfer component of agent velocity along contact normal into puck
                vproj=np.sum(self.av[:,i]*nrm,axis=1,keepdims=True)
                give=np.maximum(vproj,0)*nrm
                self.pv=np.where(hit[:,None], self.pv+give*(1.0/PUCK_MASS)*n, self.pv)
                self.av[:,i]=np.where(hit[:,None], self.av[:,i]-give*0.5, self.av[:,i])
        self.pv=_clip_v(self.pv*(1-DRAG*0.5))
        self.pp=np.clip(self.pp+self.pv*DT,0,BOX)
        # reward = ONLY reduction in puck-goal distance (dense shaping toward the ONLY objective)
        dgoal=np.linalg.norm(self.pp-self.goal,axis=1)
        # dense reward for puck VELOCITY toward goal (moving the puck goalward, achievable only via momentum
        # transfer; a slow push earns it too, so the striking MECHANISM is not itself rewarded).
        to_goal=(self.goal-self.pp); to_goal=to_goal/np.maximum(np.linalg.norm(to_goal,axis=1,keepdims=True),1e-6)
        puck_prog=np.sum(self.pv*to_goal,axis=1)                 # puck speed toward goal
        dag=np.stack([np.linalg.norm(self.pp-self.ap[:,i],axis=1) for i in range(n)],1)
        rew_team=puck_prog*0.3 + (dgoal<GOAL_R).astype(np.float32)*2.0
        rew=np.repeat(rew_team[:,None],n,axis=1) - dag*0.02      # keep mild engagement pull
        return rew
    def score(self):
        # fraction of final states with puck in goal
        return float((np.linalg.norm(self.pp-self.goal,axis=1)<GOAL_R).mean())

class Pol(nn.Module):
    def __init__(self, odim):
        super().__init__(); self.recurrent=True
        self.gru=nn.GRU(odim,H,batch_first=True); self.pi=nn.Linear(H,NMOVE); self.v=nn.Linear(H,1)
    def forward(self,x,h=None):
        z,h=self.gru(x,h); return self.pi(z),self.v(z).squeeze(-1),h

def rollout(net, env, greedy=False):
    B,n=env.B,env.n; obs=env.obs(); h=[None]*n; memwipe=getattr(env,"memwipe",False)
    O=np.zeros((B,n,T,env.odim),np.float32); MV=np.zeros((B,n,T),int)
    LP=np.zeros((B,n,T),np.float32); V=np.zeros((B,n,T),np.float32); R=np.zeros((B,n,T),np.float32)
    for t in range(T):
        mv=np.zeros((B,n),int)
        for i in range(n):
            hin=None if memwipe else h[i]
            with torch.no_grad(): lg,v,h[i]=net(torch.from_numpy(obs[:,i])[:,None,:],hin)
            d=torch.distributions.Categorical(logits=lg[:,0]); a=lg[:,0].argmax(1) if greedy else d.sample()
            O[:,i,t]=obs[:,i]; MV[:,i,t]=a.numpy(); LP[:,i,t]=d.log_prob(a).numpy(); V[:,i,t]=v[:,0].numpy()
            mv[:,i]=a.numpy()
        r=env.step(mv); R[:,:,t]=r; obs=env.obs()
    return O,MV,LP,V,R

def gae(R,V,gamma=0.97,lam=0.95):
    B,Tt=R.shape; adv=np.zeros((B,Tt),np.float32); last=np.zeros(B,np.float32)
    for t in reversed(range(Tt)):
        nv=V[:,t+1] if t+1<Tt else np.zeros(B,np.float32)
        delta=R[:,t]+gamma*nv-V[:,t]; last=delta+gamma*lam*last; adv[:,t]=last
    return adv, adv+V

def train(iters, seed=0, n_ag=2, B=128):
    torch.manual_seed(seed); od=MomentumWorld(1,0,n_ag=n_ag).odim
    net=Pol(od); opt=torch.optim.Adam(net.parameters(),lr=3e-3)
    for it in range(iters):
        # curriculum: first 40% of training spawns agents near the puck (learn contact), then full box
        near = it < int(0.4*iters)
        env=MomentumWorld(B,1000+it,n_ag=n_ag,spawn_near=near)
        O,MV,LP,V,R=rollout(net,env)
        O2=O.reshape(B*n_ag,T,od); MV2=MV.reshape(B*n_ag,T); LP2=LP.reshape(B*n_ag,T)
        V2=V.reshape(B*n_ag,T); R2=R.reshape(B*n_ag,T)
        adv,ret=gae(R2,V2); adv=(adv-adv.mean())/(adv.std()+1e-8)
        Ot,MVt,LPt=torch.from_numpy(O2),torch.from_numpy(MV2),torch.from_numpy(LP2)
        advt,rett=torch.from_numpy(adv),torch.from_numpy(ret)
        for _ in range(4):
            lg,v,_=net(Ot); d=torch.distributions.Categorical(logits=lg)
            ratio=torch.exp(d.log_prob(MVt)-LPt); s1=ratio*advt; s2=torch.clamp(ratio,0.8,1.2)*advt
            loss=-torch.min(s1,s2).mean()+0.5*((v-rett)**2).mean()-0.03*d.entropy().mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(),0.5); opt.step()
        if it%100==0 or it==iters-1:
            ev=MomentumWorld(400,7000+it,n_ag=n_ag); rollout(net,ev,greedy=True)
            print(f"  iter {it:>3}: puck-in-goal {ev.score():.2f}",flush=True)
    return net

def _score(net, seed, B=400, **kw):
    env=MomentumWorld(B,seed,**kw); rollout(net,env,greedy=True); return env.score()
