import numpy as np, csv, json
from scipy.special import erf, erfinv
from pathlib import Path
N=1024; J=10; R=100; KS=[1,3,5,10]
OUT=Path(__file__).parent

def wavmat(h):
 g=(h*(-1.)**np.arange(1,len(h)+1))[::-1];h=np.r_[h,np.zeros(N)];g=np.r_[g,np.zeros(N)];W=np.eye(8)
 for k in range(7,0,-1):
  m=2**(J-k);ix=(N+np.arange(1,2*m+1)[:,None]-2*np.arange(1,m+1)[None,:]+2)%(2*m);ix=np.where(ix==0,2*m,ix)-1;W=np.vstack((W@h[ix].T,g[ix].T))
 return W

def gate(x,fi):
 if fi==0:
  g=np.tanh(x); gp=1-g*g
 elif fi==1:g=erf(x/np.sqrt(2));gp=np.sqrt(2/np.pi)*np.exp(-x*x/2)
 elif fi==2:g=-np.expm1(-x);gp=np.exp(-x)
 elif fi==3:
  e=np.exp(-np.pi*x/2);g=1-4/np.pi*np.arctan(e);gp=2*e/(1+e*e)
 elif fi==4:g=np.minimum(x,1);gp=(x<1).astype(float)
 else:g=2/np.pi*np.arctan(x);gp=2/np.pi/(1+x*x)
 return g,gp

def quantile(k,fi):
 r=2**(-1/k)
 return [np.arctanh(r),np.sqrt(2)*erfinv(r),-np.log1p(-r),2/np.pi*np.arcsinh(np.tan(np.pi*r/2)),r,np.tan(np.pi*r/2)][fi]

# Selection function receives NO truth. Integer grid indexes avoid interpolation.
# alpha=0:.05:1.5; log2(c)=-3:.05:10.
def select(d,sig,fi):
 aa=np.arange(31); cc=np.arange(-60,201); ai,ci=np.meshgrid(aa,cc,indexing='ij');ai=ai.ravel();ci=ci.ravel()
 bs=np.arange(-330,201); rows=[]
 for k in KS:
  risk=np.zeros((len(ai),R));q=quantile(k,fi)
  for j in range(3,10):
   z=d[2**j:2**(j+1)];x0=q*abs(z)/sig;curve=np.empty((len(bs),R))
   for bi,b in enumerate(bs):
    x=x0*2**(b/20);g,gp=gate(x,fi);att=g**k;der=att+k*x*g**(k-1)*gp
    curve[bi]=np.sum((z*(att-1))**2+2*sig**2*der-sig**2,axis=0)
   risk+=curve[ci-ai*j+330]
  choice=np.argmin(risk,axis=0);rows.append((risk[choice,np.arange(R)],ai[choice]/20,ci[choice]/20,k))
 best=np.argmin(np.array([r[0] for r in rows]),axis=0)
 fixed=rows[KS.index(5)]
 chosen=(np.array([rows[b][0][i] for i,b in enumerate(best)]),np.array([rows[b][1][i] for i,b in enumerate(best)]),np.array([rows[b][2][i] for i,b in enumerate(best)]),np.array([rows[b][3] for b in best]))
 return fixed,chosen

def estimate(d,sig,fi,selection):
 _,a,lc,k=selection;k=np.broadcast_to(k,(R,));q=np.array([quantile(v,fi) for v in k]);est=d.copy()
 for j in range(3,10):
  sl=slice(2**j,2**(j+1));x=abs(d[sl])/sig*q*2**(lc-a*j);g,_=gate(x,fi);est[sl]*=g**k
 return est

if __name__=='__main__':
 t=np.arange(1,N+1)/N;p=np.array([.1,.13,.15,.23,.25,.4,.44,.65,.76,.78,.81]);bl=((1+np.sign(t[:,None]-p))/2)@np.array([4,-5,3,-4,5,-4.2,2.1,4.3,-3.1,2.1,-4.2]);bu=np.sum(np.array([4,5,3,4,5,4.2,2.1,4.3,3.1,5.1,4.2])/(1+abs((t[:,None]-p)/np.array([.005,.005,.006,.01,.01,.03,.01,.01,.005,.008,.005])))**4,axis=1);he=4*np.sin(4*np.pi*t)-np.sign(t-.3)-np.sign(.72-t);dop=np.sqrt(t*(1-t))*np.sin(2*np.pi*1.05/(t+.05))
 haar=np.array([1,1])/np.sqrt(2);db=np.array([.3326705529500826,.8068915093110928,.4598775021184915,-.1350110200102546,-.0854412738820267,.0352262918857095]);sy=np.array([-.0757657147892733,-.0296355276459985,.4976186676320155,.8037387518059161,.2978577956052774,-.0992195435768472,-.0126039672620378,.0322231006040427]);eps=np.random.default_rng(20261024).normal(size=(N,R));raw=[];summary=[]
 for sn,f,h in zip(['Blocks','Bumps','HeaviSine','Doppler'],[bl,bu,he,dop],[haar,db,sy,sy]):
  f=f*np.sqrt(5/np.mean(f*f));W=wavmat(h);assert np.max(abs(W@W.T-np.eye(N)))<1e-10;d=W@(f[:,None]+eps);sig=np.median(abs(d[512:]),axis=0)/.6744897501960817
  for fi,fn in enumerate(['Logistic','Normal','Laplace','Sech','Uniform','Cauchy']):
   choices=select(d,sig,fi)
   for mode,s in zip(['k5','select_k'],choices):
    est=estimate(d,sig,fi,s);ms=np.mean((W.T@est-f[:,None])**2,axis=0);k=np.broadcast_to(s[3],(R,));boundary=(s[1]==0)|(s[1]==1.5)|(s[2]==-3)|(s[2]==10)
    row=dict(signal=sn,F=fn,mode=mode,AMSE=float(ms.mean()),MCSE=float(ms.std(ddof=1)/10),median_alpha=float(np.median(s[1])),median_c=float(np.median(2**s[2])),boundary_count=int(boundary.sum()),k_counts=str(dict(zip(*np.unique(k,return_counts=True)))))
    summary.append(row);print(json.dumps(row),flush=True)
    for r in range(R):raw.append(dict(signal=sn,F=fn,mode=mode,rep=r+1,MSE=ms[r],alpha=s[1][r],c=2**s[2][r],k=k[r],sigmahat=sig[r]))
  for name,data in [('SCOPE_SURE_summary.csv',summary),('SCOPE_SURE_replications.csv',raw)]:
   with open(OUT/name,'w',newline='') as out:
    w=csv.DictWriter(out,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
