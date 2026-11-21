function [Summary, Replications] = SCOPE_SURE_test()
% Plug-in SURE selection; no truth enters select_sure.
% Uses Wavmat(h,N,J-3,2) if supplied; otherwise the embedded Wavmat algorithm.
% No wavelet or statistics toolbox required. MATLAB R2016b or later.
% Seed matches Python numerically, but random streams differ across languages.
N=1024; J=10; R=100; rng(20261024,'twister'); E=randn(N,R);
t=(1:N)'/N; p=[.1 .13 .15 .23 .25 .4 .44 .65 .76 .78 .81];
bl=((1+sign(t-p))/2)*[4 -5 3 -4 5 -4.2 2.1 4.3 -3.1 2.1 -4.2]';
bu=sum([4 5 3 4 5 4.2 2.1 4.3 3.1 5.1 4.2]./(1+abs((t-p)./[.005 .005 .006 .01 .01 .03 .01 .01 .005 .008 .005])).^4,2);
he=4*sin(4*pi*t)-sign(t-.3)-sign(.72-t);
dop=sqrt(t.*(1-t)).*sin(2*pi*1.05./(t+.05));
H={[1 1]/sqrt(2),[.3326705529500826 .8068915093110928 .4598775021184915 -.1350110200102546 -.0854412738820267 .0352262918857095],[-.0757657147892733 -.0296355276459985 .4976186676320155 .8037387518059161 .2978577956052774 -.0992195435768472 -.0126039672620378 .0322231006040427]};
H{4}=H{3}; signals={bl,bu,he,dop}; sn={'Blocks','Bumps','HeaviSine','Doppler'};
fn={'Logistic','Normal','Laplace','Sech','Uniform','Cauchy'}; modes={'k5','select_k'};
srows={}; rrows={};
for s=1:4
 f=signals{s};f=f*sqrt(5/mean(f.^2));
 if exist('Wavmat','file')==2, W=Wavmat(H{s},N,J-3,2); else, W=wavmat_embedded(H{s},N,J-3,2); end
 assert(max(max(abs(W*W'-eye(N))))<1e-10,'Transform is not orthogonal');
 d=W*(f+E); sigmahat=median(abs(d(N/2+1:end,:)),1)/.6744897501960817;
 for fi=1:6
  choices=select_sure(d,sigmahat,fi);
  for mode=1:2
   pars=choices{mode}; alpha=pars(1,:); lc=pars(2,:); k=pars(3,:);
   q=arrayfun(@(kk) qvalue(kk,fi),k); est=d;
   for j=3:J-1
    ix=2^j+1:2^(j+1); x=abs(d(ix,:))./sigmahat.*q.*2.^(lc-alpha*j);
    [g,~]=gate(x,fi);est(ix,:)=d(ix,:).*g.^k;
   end
   mse=mean((W'*est-f).^2,1); boundary=sum(alpha==0 | alpha==1.5 | lc==-3 | lc==10);
   srows(end+1,:)={sn{s},fn{fi},modes{mode},mean(mse),std(mse)/sqrt(R),median(alpha),median(2.^lc),boundary}; %#ok<AGROW>
   for r=1:R
    rrows(end+1,:)={sn{s},fn{fi},modes{mode},r,mse(r),alpha(r),2^lc(r),k(r),sigmahat(r)}; %#ok<AGROW>
   end
   fprintf('%s %-8s %-8s AMSE %.6f SE %.6f\n',sn{s},fn{fi},modes{mode},mean(mse),std(mse)/sqrt(R));
  end
 end
 Summary=cell2table(srows,'VariableNames',{'Signal','F','Mode','AMSE','MCSE','MedianAlpha','MedianC','BoundaryCount'});
 Replications=cell2table(rrows,'VariableNames',{'Signal','F','Mode','Replication','MSE','Alpha','c','k','SigmaHat'});
 writetable(Summary,'SCOPE_SURE_MATLAB_summary.csv');writetable(Replications,'SCOPE_SURE_MATLAB_replications.csv');
 save('SCOPE_SURE_MATLAB_results.mat','Summary','Replications','N','R');
end
end

function choices=select_sure(d,sig,fi)
% Parameters selected using only noisy coefficients and sigma estimate.
% Exact grid: alpha=0:.05:1.5; log2(c)=-3:.05:10; k in [1 3 5 10].
R=size(d,2);[ai,ci]=ndgrid(0:30,-60:200);ai=ai(:);ci=ci(:);bs=-330:200;
ks=[1 3 5 10];scores=zeros(4,R);pars=cell(4,1);
for ki=1:4
 k=ks(ki);q=qvalue(k,fi);risk=zeros(numel(ai),R);
 for j=3:9
  ix=2^j+1:2^(j+1);z=d(ix,:);x0=q*abs(z)./sig;curve=zeros(numel(bs),R);
  for bi=1:numel(bs)
   x=x0*2^(bs(bi)/20);[g,gp]=gate(x,fi);att=g.^k;der=att+k*x.*g.^(k-1).*gp;
   curve(bi,:)=sum((z.*(att-1)).^2+2*sig.^2.*der-sig.^2,1);
  end
  risk=risk+curve(ci-ai*j+331,:);
 end
 [scores(ki,:),idx]=min(risk,[],1);
 pars{ki}=[reshape(ai(idx),1,R)/20;reshape(ci(idx),1,R)/20;k*ones(1,R)];
end
[~,best]=min(scores,[],1);chosen=zeros(3,R);
for r=1:R, chosen(:,r)=pars{best(r)}(:,r); end
choices={pars{3},chosen};
end

function q=qvalue(k,fi)
r=2^(-1/k);qs=[atanh(r),sqrt(2)*erfinv(r),-log1p(-r),2/pi*asinh(tan(pi*r/2)),r,tan(pi*r/2)];q=qs(fi);
end

function [g,gp]=gate(x,fi)
switch fi
 case 1, g=tanh(x);gp=1-g.^2;
 case 2, g=erf(x/sqrt(2));gp=sqrt(2/pi)*exp(-x.^2/2);
 case 3, g=-expm1(-x);gp=exp(-x);
 case 4, e=exp(-pi*x/2);g=1-4/pi*atan(e);gp=2*e./(1+e.^2);
 case 5, g=min(x,1);gp=double(x<1);
 case 6, g=2/pi*atan(x);gp=2/pi./(1+x.^2);
end
end

function W=wavmat_embedded(h,N,k0,shift)
% Algorithm from Brani Vidakovic's Wavmat, 2004, retaining coarse dimension 8.
J=log2(N);h=h(:)';g=fliplr(conj(h).*(-1).^(1:length(h)));
h=[h zeros(1,N)];g=[g zeros(1,N)];oldmat=eye(2^(J-k0));
for k=k0:-1:1
 m=2^(J-k);[ii,jj]=ndgrid(1:2*m,1:m);ix=mod(N+ii-2*jj+shift,2*m);ix(ix==0)=2*m;
 hmat=h(ix);gmat=g(ix);W=[oldmat*hmat';gmat'];oldmat=W;
end
end
