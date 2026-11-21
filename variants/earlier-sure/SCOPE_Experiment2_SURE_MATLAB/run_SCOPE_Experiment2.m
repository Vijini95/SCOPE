function run_SCOPE_Experiment2(mode)
% RUN_SCOPE_EXPERIMENT2   SCOPE SURE and N-SCOPE SURE, Section 6 Experiment II.
%
% Full manuscript experiment (4 DJ signals x 3 SNRs x M=1000):
%       run_SCOPE_Experiment2
% Quick functionality check, NOT publishable:
%       run_SCOPE_Experiment2('smoke')
%
% SELF-CONTAINED: Needs only Wavmat.m and MakeONFilterExt.m included here.
% No MATLAB Wavelet Toolbox, Statistics Toolbox, WaveLab, or Optimization
% Toolbox required. MATLAB with local functions, tables, and export graphics
% via print() is sufficient (R2018a or newer recommended).
%
% DESIGN: n=1024, t_i=i/n, Gaussian noise sigma=1, power SNR={3,5,7}.
% Periodic orthonormal DWT with 8 unshrunk scaling coefficients:
%  Doppler, HeaviSine: Symmlet-4 (8 taps)
%  Bumps:              Daubechies length 6
%  Blocks:             Haar
% Details j=3,...,9. SigmaHat=median(abs(finest details))/0.674489750196.
%
% SIX CDFs: Normal, Logistic, Laplace, Sech, Uniform, Cauchy.
% ADAPTIVE GRIDS: k={3,5,7}; r={0,1,2,3}; b={0.65,0.8,1,1.2,1.5}.
% lambda_{j,r}(F,k,b) = b*a_{F,r,k}(eta_j) / SigmaHat,
% where eta_j=0.01*2^(9-j), and a solves the *exact* Gaussian null-risk
% equation with U ~ ChiSquare(2r+1):
%   E[(U/(2r+1)) * Fstar(a*sqrt(U/(2r+1)))^(2*k)] = eta_j.
%
% For each CDF and replicate, SCOPE_SURE MINIMIZES the summed plug-in SURE
% over all (k,b) at r=0; NSCOPE_SURE minimizes over all (k,r,b).
% All tuning selections are GLOBAL across levels j=3:9, not levelwise.
% The score contains Stein's diagonal derivative, INCLUDING the dependence
% of central RMS on the central coefficient; scale/SigmaHat held fixed in
% the derivative. The minimum SURE is not an unbiased estimate of the risk
% after selection; we evaluate squared error against truth only AFTER choice.
%
% Two FIXED baselines use k=5,b=1,r=0 or r=1 and the same CDF calibration.
% Hard/soft universal threshold reference methods are also evaluated.
% Matched noise streams are seeded as in Experiment I, facilitating checks.
%
% RESUME: a saved partial checkpoint after every batch permits restarting
% the current signal-SNR combination without repeating completed batches.
% A separate completed checkpoint is saved after each signal-SNR condition.
%
% OUTPUT: outputs_E2_M1000 or outputs_SMOKE_NOT_FOR_PUBLICATION. Consult
% README.txt for CSV, MAT, and plotted results.

if nargin<1 || isempty(mode), mode='full'; end
if ~ischar(mode), mode=char(mode); end
smoke=strcmpi(mode,'smoke');
assert(smoke || strcmpi(mode,'full'), 'Mode must be ''smoke'' or ''full''.');
rootDir=fileparts(mfilename('fullpath'));
addpath(rootDir);

cfg.N=1024; cfg.J=10; cfg.j0=3;
cfg.kGrid=[3 5 7];
cfg.rGrid=0:3;
cfg.bGrid=[0.65 0.80 1.00 1.20 1.50];
cfg.eta=0.01;
cfg.generators={'Normal','Logistic','Laplace','Sech','Uniform','Cauchy'};
cfg.signals={'Doppler','HeaviSine','Bumps','Blocks'};
cfg.snrs=[3 5 7];
cfg.M=1000; cfg.batchSize=40;
cfg.seed=20261026; cfg.shift=2;
cfg.makePlots=true; cfg.resume=true;
if smoke
    cfg.generators={'Normal','Logistic'};
    cfg.signals={'Doppler','Blocks'};
    cfg.snrs=5;
    cfg.M=6; cfg.batchSize=3;
    outDir=fullfile(rootDir,'outputs_SMOKE_NOT_FOR_PUBLICATION');
    calibFile=fullfile(rootDir,'cached_E2_calibration_SMOKE.mat');
else
    outDir=fullfile(rootDir,'outputs_E2_M1000');
    calibFile=fullfile(rootDir,'cached_E2_calibration_FULL.mat');
end
if ~exist(outDir,'dir'), mkdir(outDir); end
plotDir=fullfile(outDir,'plots');
if ~exist(plotDir,'dir'), mkdir(plotDir); end

cand=make_candidates(cfg);
validate_divergence(); % finite-difference check for all six CDF derivatives
names=method_names(cfg.generators);
fCount=numel(cfg.generators); cCount=numel(cand.k);
assert(cCount==numel(cfg.kGrid)*numel(cfg.rGrid)*numel(cfg.bGrid));
fprintf('Experiment II (%s): M=%d, %d CDFs, %d candidates per CDF.\n',...
    mode,cfg.M,fCount,cCount);

% Exactly the same quadrature convention as Experiment I, except k varies.
% Full calibration: 6 x 3 k x 4 radii x 7 levels = 504 constants.
calibMeta=struct('generators',{cfg.generators},'kGrid',cfg.kGrid,...
    'rGrid',cfg.rGrid,'eta',cfg.eta,'j0',cfg.j0,'J',cfg.J);
if exist(calibFile,'file')
    z=load(calibFile,'a','calibMeta');
    assert(isfield(z,'a') && isfield(z,'calibMeta') && ...
       isequaln(z.calibMeta,calibMeta),...
       'The cached calibration belongs to different settings: %s',calibFile);
    a=z.a;
    fprintf('Loaded null-risk constants from %s\n',calibFile);
else
    fprintf('Calibrating Gaussian null-risk constants (one-time task)...\n');
    a=compute_constants(cfg);
    save(calibFile,'a','calibMeta','-v7');
end
assert(isequal(size(a),[fCount,numel(cfg.kGrid),numel(cfg.rGrid),7]));
check_k5_constants(a,cfg);
write_constants(cfg,a,fullfile(outDir,'Experiment2_ScaleConstants.csv'));

% Build the supplied periodic orthonormal wavelet matrices ONCE.
Wbank=cell(1,3);
Wbank{1}=wavelet_matrix('Symmlet',4,cfg);
Wbank{2}=wavelet_matrix('Daubechies',6,cfg);
Wbank{3}=wavelet_matrix('Haar',0,cfg);
results=cell(numel(cfg.signals),numel(cfg.snrs));
totalTimer=tic;
for is=1:numel(cfg.signals)
    signalName=cfg.signals{is};
    [f0,filterName,whichW]=make_DJ_signal(signalName,cfg.N);
    W=Wbank{whichW};
    for iq=1:numel(cfg.snrs)
        snr=cfg.snrs(iq);
        fn=fullfile(outDir,sprintf('checkpoint_%s_SNR%d.mat',signalName,snr));
        partialFn=fullfile(outDir,sprintf('PARTIAL_%s_SNR%d.mat',signalName,snr));
        if cfg.resume && exist(fn,'file')
            tmp=load(fn,'result');
            result=tmp.result;
            assert(isequal(size(result.mse),[cfg.M,numel(names)]) && ...
                isfield(result,'studySettings') && ...
                isequaln(result.studySettings,study_settings(cfg)),...
                'Existing completed checkpoint is incompatible with the settings.');
            fprintf('RESUMED completed: %s SNR=%g\n',signalName,snr);
        else
            f=f0*sqrt(snr/mean(f0.^2));
            fprintf('\nStarting %s | SNR=%g | M=%d\n',signalName,snr,cfg.M);
            rng(cfg.seed+100*is+iq,'twister');
            result=simulate_condition(f,W,a,cfg,cand,names,...
                signalName,filterName,snr,partialFn);
            save(fn,'result','-v7');
            if exist(partialFn,'file'), delete(partialFn); end
        end
        results{is,iq}=result;
        write_mse_csv(result,fullfile(outDir,sprintf(...
            'MSE_replications_%s_SNR%d.csv',signalName,snr)));
        write_selection_csv(result,fullfile(outDir,sprintf(...
            'SelectedParams_%s_SNR%d.csv',signalName,snr)));
        if cfg.makePlots
            condition_plots(result,cfg,plotDir);
        end
        fprintf('Done: %s SNR=%g; total elapsed %.1f sec\n',...
            signalName,snr,toc(totalTimer));
    end
end

[summary,paired,selFreq,tripleFreq,ranking]=compile_tables(results,cfg,names);
writetable(summary,fullfile(outDir,'Experiment2_AMSE_all.csv'));
writetable(paired,fullfile(outDir,'Experiment2_PairedDifferences.csv'));
writetable(selFreq,fullfile(outDir,'Experiment2_SelectionFrequencies.csv'));
writetable(tripleFreq,fullfile(outDir,'Experiment2_SelectedTriples.csv'));
writetable(ranking,fullfile(outDir,'Experiment2_AggregateRanking.csv'));
if cfg.makePlots
    aggregate_plot(ranking,cfg,plotDir);
end
save(fullfile(outDir,'Experiment2_Results.mat'),...
    'cfg','a','cand','results','-v7.3');
write_results_readme(outDir,cfg);
fprintf('\nExperiment II finished: %s\n',outDir);
fprintf('AMSE: Experiment2_AMSE_all.csv | Paired: Experiment2_PairedDifferences.csv\n');
fprintf('Total elapsed %.1f seconds.\n',toc(totalTimer));
end

% ========================================================================
function cand=make_candidates(cfg)
% CANDIDATE ORDER: radius outer, then k, then b. r=0 candidates are the
% exact SCOPE SURE candidate set and are REUSED inside N-SCOPE SURE.
nC=numel(cfg.rGrid)*numel(cfg.kGrid)*numel(cfg.bGrid);
cand.r=zeros(nC,1); cand.k=zeros(nC,1);
cand.b=zeros(nC,1); cand.ir=zeros(nC,1); cand.ik=zeros(nC,1);
t=0;
for ir=1:numel(cfg.rGrid)
 for ik=1:numel(cfg.kGrid)
  for ib=1:numel(cfg.bGrid)
   t=t+1; cand.r(t)=cfg.rGrid(ir); cand.k(t)=cfg.kGrid(ik);
   cand.b(t)=cfg.bGrid(ib); cand.ir(t)=ir; cand.ik(t)=ik;
  end
 end
end
cand.scopeIndices=find(cand.r==0);
cand.fixedScope=find(cand.r==0 & cand.k==5 & abs(cand.b-1)<1e-12);
cand.fixedNScope=find(cand.r==1 & cand.k==5 & abs(cand.b-1)<1e-12);
assert(isscalar(cand.fixedScope) && isscalar(cand.fixedNScope));
end

% ========================================================================
function W=wavelet_matrix(filterName,par,cfg)
h=MakeONFilterExt(filterName,par);
assert(ismember(numel(h),[2 6 8]));
W=Wavmat(h,cfg.N,cfg.J-cfg.j0,cfg.shift);
assert(isequal(size(W),[cfg.N cfg.N]));
v=sin((1:cfg.N)'*1.2345);
err=norm(W'*(W*v)-v)/norm(v);
assert(err<3e-8,'Wavmat not orthogonal to numerical precision.');
fprintf('Wavelet %-11s %d taps; reconstruction check %.3g\n',...
    filterName,numel(h),err);
end

% ========================================================================
function [f0,label,whichW]=make_DJ_signal(name,N)
% MATCHES Experiment I *exactly*, including sampling grid and DJ heights.
t=(1:N)'/N;
pos=[0.10 0.13 0.15 0.23 0.25 0.40 0.44 0.65 0.76 0.78 0.81];
htB=[4 5 3 4 5 4.2 2.1 4.3 3.1 5.1 4.2];
htK=[4 -5 3 -4 5 -4.2 2.1 4.3 -3.1 2.1 -4.2];
w=[.005 .005 .006 .010 .010 .030 .010 .010 .005 .008 .005];
switch lower(name)
 case 'doppler'
  f0=sqrt(t.*(1-t)).*sin(2*pi*1.05./(t+0.05));
  label='Symmlet4_8tap'; whichW=1;
 case 'heavisine'
  f0=4*sin(4*pi*t)-sign(t-.30)-sign(.72-t);
  label='Symmlet4_8tap'; whichW=1;
 case 'bumps'
  f0=zeros(size(t));
  for ii=1:numel(pos), f0=f0+htB(ii)./(1+abs((t-pos(ii))/w(ii))).^4; end
  label='Daubechies3_6tap'; whichW=2;
 case 'blocks'
  f0=zeros(size(t));
  for ii=1:numel(pos), f0=f0+htK(ii)*(t>=pos(ii)); end
  label='Haar_2tap'; whichW=3;
 otherwise
  error('Unknown DJ signal %s',name);
end
end

% ========================================================================
function a=compute_constants(cfg)
% Gaussian-null per-coefficient risk; no truth/signal enters calibration.
% Original Experiment I uses integral(...,0,Inf) and fzero. Here repeated
% for k={3,5,7} to ensure identical k=5 constants.
nF=numel(cfg.generators); nK=numel(cfg.kGrid);
nR=numel(cfg.rGrid); nJ=cfg.J-cfg.j0;
a=zeros(nF,nK,nR,nJ);
for iF=1:nF
 for ik=1:nK
  for ir=1:nR
   for jj=1:nJ
    j=cfg.j0+jj-1;
    eta_j=cfg.eta*2^(cfg.J-1-j);
    fun=@(x)nullrisk(x,cfg.generators{iF},cfg.rGrid(ir),cfg.kGrid(ik));
    upper=1;
    while fun(upper)<eta_j
       upper=upper*2;
       assert(upper<1e6,'Null-risk root failed to bracket.');
    end
    a(iF,ik,ir,jj)=fzero(@(x)fun(x)-eta_j,[0 upper],...
          optimset('TolX',1e-10,'Display','off'));
   end
  end
 end
 fprintf('  Calibrated %-9s | k={3,5,7}, r=0:3\n',cfg.generators{iF});
end
end

function risk=nullrisk(a,F,r,k)
if a<=0, risk=0; return; end
m=2*r+1; nu=m/2;
logC=-nu*log(2)-gammaln(nu)-log(m);
fun=@(u)exp(nu.*log(max(u,realmin))-.5*u+logC) ...
    .*cdf_centered(F,a.*sqrt(u/m)).^(2*k);
risk=integral(fun,0,Inf,'AbsTol',2e-10,'RelTol',1e-8);
risk=max(0,min(1,risk));
end

function [g,gprime]=cdf_centered(F,u)
% Positive arguments only (local RMS >= 0). Conventions are EXACTLY those
% in SCOPE6, Table 1: logistic is TANH(u), not tanh(u/2).
u=max(0,u);
switch lower(F)
 case 'normal'
    g=erf(u/sqrt(2));
    if nargout>1, gprime=sqrt(2/pi)*exp(-u.^2/2); end
 case 'logistic'
    g=tanh(u);
    if nargout>1, gprime=1-g.^2; end
 case 'laplace'
    g=-expm1(-u);
    if nargout>1, gprime=exp(-u); end
 case 'sech'
    g=(2/pi)*atan(sinh(pi*u/2));
    if nargout>1, gprime=1./cosh(pi*u/2); end
 case 'uniform'
    g=min(u,1);
    if nargout>1, gprime=double(u<1); end % derivative a.e.
 case 'cauchy'
    g=(2/pi)*atan(u);
    if nargout>1, gprime=(2/pi)./(1+u.^2); end
 otherwise
    error('Unknown F %s',F);
end
g=max(0,min(g,1));
end

% ========================================================================
function validate_divergence()
% Numerically verify the *diagonal* Jacobian term used in N-SCOPE SURE.
% The other observed coefficients, sigma, and lambda are held fixed.
% Uniform kink locations are avoided. This is a deterministic unit test.
Fs={'Normal','Logistic','Laplace','Sech','Uniform','Cauchy'};
d=[0.24; -0.74; 0.51; 1.11; -0.42; 0.66; -0.35; 0.83];
center=4; lambda=0.57; epsFD=1e-6;
for f=1:numel(Fs)
 for rr=[0 1 3]
  mm=2*rr+1;
  ids=mod((center-rr:center+rr)-1,numel(d))+1;
  S=sqrt(sum(d(ids).^2)/mm);
  u=lambda*S;
  [g,gp]=cdf_centered(Fs{f},u);
  for k=[3 5 7]
   analytic=g^k+(k*lambda*d(center)^2/(mm*S))*g^(k-1)*gp;
   plus=d; plus(center)=plus(center)+epsFD;
   minus=d; minus(center)=minus(center)-epsFD;
   sp=sqrt(sum(plus(ids).^2)/mm);
   sm=sqrt(sum(minus(ids).^2)/mm);
   fp=plus(center)*cdf_centered(Fs{f},lambda*sp)^k;
   fm=minus(center)*cdf_centered(Fs{f},lambda*sm)^k;
   numerical=(fp-fm)/(2*epsFD);
   assert(abs(numerical-analytic)<5e-6,...
      'SURE derivative validation failed: F=%s, k=%d, r=%d',Fs{f},k,rr);
  end
 end
end
fprintf('Passed finite-difference checks of N-SCOPE diagonal divergence.\n');
end

% ========================================================================
function check_k5_constants(a,cfg)
% Same 24 k=5, j=9 null-risk constants as Experiment I/SCOPE6.
expected=[.4071 .5865 .6646 .7103; ...
          .3443 .4932 .5571 .5941; ...
          .4826 .6840 .7678 .8153; ...
          .3627 .5175 .5832 .6209; ...
          .2545 .3766 .4341 .4693; ...
          .7936 1.1018 1.2232 1.2900];
idx5=find(cfg.kGrid==5);
for iF=1:numel(cfg.generators)
 if iF>size(expected,1), continue; end
 allF={'Normal','Logistic','Laplace','Sech','Uniform','Cauchy'};
 pos=find(strcmp(cfg.generators{iF},allF));
 if isempty(pos), continue; end
 x=reshape(a(iF,idx5,:,end),1,[]);
 d=max(abs(x-expected(pos,:)));
 assert(d<.002,'k=5 null-risk constants disagree with Experiment I.');
end
fprintf('Validated k=5, j=9 calibration against Experiment I.\n');
end

function write_constants(cfg,a,fileName)
F={}; K=[]; R=[]; J=[]; ETA=[]; A=[];
for f=1:numel(cfg.generators)
 for ik=1:numel(cfg.kGrid)
  for ir=1:numel(cfg.rGrid)
   for jj=1:(cfg.J-cfg.j0)
    j=cfg.j0+jj-1;
    F{end+1,1}=cfg.generators{f}; %#ok<AGROW>
    K(end+1,1)=cfg.kGrid(ik); %#ok<AGROW>
    R(end+1,1)=cfg.rGrid(ir); %#ok<AGROW>
    J(end+1,1)=j; %#ok<AGROW>
    ETA(end+1,1)=cfg.eta*2^(cfg.J-1-j); %#ok<AGROW>
    A(end+1,1)=a(f,ik,ir,jj); %#ok<AGROW>
   end
  end
 end
end
T=table(F,K,R,J,ETA,A,'VariableNames',...
   {'Generator','K','Radius','Level_j','Eta_j','a_NullRisk'});
writetable(T,fileName);
end

% ========================================================================
function S=study_settings(cfg)
% Protect checkpoints from accidental reuse after tuning-grid changes.
S=struct('N',cfg.N,'M',cfg.M,'J',cfg.J,'j0',cfg.j0,...
  'eta',cfg.eta,'kGrid',cfg.kGrid,'rGrid',cfg.rGrid,...
  'bGrid',cfg.bGrid,'generators',{cfg.generators},...
  'batchSize',cfg.batchSize,'seed',cfg.seed,'shift',cfg.shift);
end

function result=simulate_condition(f,W,a,cfg,cand,names,tag,filterName,snr,partialFn)
N=cfg.N; M=cfg.M; nbMethods=numel(names);
nF=numel(cfg.generators); nbCandidates=numel(cand.k);

% State contains all in-progress values, including the exact RNG state.
% This makes restarting mid-condition reproducible and deterministic.
if cfg.resume && exist(partialFn,'file')
    z=load(partialFn,'state'); state=z.state;
    assert(state.M==M && state.nF==nF && state.nbCandidates==nbCandidates && ...
        isfield(state,'studySettings') && ...
        isequaln(state.studySettings,study_settings(cfg)),...
        'Incompatible partial checkpoint. Move/delete %s',partialFn);
    rng(state.rngState);
    fprintf('  Resuming PARTIAL at replicate %d/%d.\n',state.completed,M);
else
    state=struct(); state.M=M; state.nF=nF;
    state.nbCandidates=nbCandidates; state.completed=0;
    state.studySettings=study_settings(cfg);
    state.mse=NaN(M,nbMethods);
    state.sigmahat=NaN(M,1);
    state.scopeIdx=zeros(M,nF);
    state.nscopeIdx=zeros(M,nF);
    state.firstReplication=[];
    state.rngState=rng;
end
trueD=W*f;
finest=(N/2+1):N;
for istart=(state.completed+1):cfg.batchSize:M
    iend=min(M,istart+cfg.batchSize-1);
    nb=iend-istart+1;
    Y=bsxfun(@plus,f,randn(N,nb));
    D=W*Y;
    sigmahat=median(abs(D(finest,:)),1)/0.6744897501960817;
    assert(all(sigmahat>0 & isfinite(sigmahat)));
    state.sigmahat(istart:iend)=sigmahat';
    if istart==1
        state.firstReplication=struct('y',Y(:,1),'d',D(:,1),...
           'sigmahat',sigmahat(1));
    end
    % Store seven wavelet detail blocks and 7x4 overlapping RMS arrays.
    % Neighborhoods are PERIODIC within each resolution level.
    db=cell(7,1); tr=cell(7,1); S=cell(7,4);
    for jj=1:7
        j=cfg.j0+jj-1;
        ind=(2^j+1):2^(j+1);
        db{jj}=D(ind,:); tr{jj}=trueD(ind);
        d2=db{jj}.^2;
        for ir=1:4
            rad=cfg.rGrid(ir);
            local=d2;
            for lag=1:rad
                local=local+circshift(d2,[lag 0])+circshift(d2,[-lag 0]);
            end
            S{jj,ir}=sqrt(local/(2*rad+1));
        end
    end
    errScaling=sum(bsxfun(@minus,D(1:2^cfg.j0,:),...
        trueD(1:2^cfg.j0)).^2,1);
    mseBatch=NaN(nb,nbMethods);

    for iF=1:nF
        F=cfg.generators{iF};
        sureScores=zeros(nbCandidates,nb);
        % Evaluate EVERY candidate, summing the SURE score over ALL 7 levels.
        % Scale is held fixed with respect to d_jl in the derivative.
        for c=1:nbCandidates
            ir=cand.ir(c); ik=cand.ik(c);
            r=cand.r(c); k=cand.k(c); b=cand.b(c);
            mm=2*r+1;
            for jj=1:7
                dj=db{jj}; Sj=S{jj,ir};
                lam=b*a(iF,ik,ir,jj)./sigmahat; % 1 x nb
                u=bsxfun(@times,Sj,lam);
                [g,gp]=cdf_centered(F,u);
                w=g.^k;
                % Derivative of each *central* estimate wrt its own datum:
                %   delta' = w + (k*lambda*d^2/(m*S)) * g^(k-1)*g'
                % d=0 and S=0 -> correction interpreted as 0 by continuity.
                correction=(k/mm)*((dj.^2)./max(Sj,realmin)) ...
                   .*g.^(k-1).*gp;
                deriv=w+bsxfun(@times,correction,lam);
                residual=dj.*(w-1);
                blockRisk=sum(residual.^2,1) + ...
                   2*(sigmahat.^2).*sum(deriv,1) - ...
                   size(dj,1)*(sigmahat.^2);
                sureScores(c,:)=sureScores(c,:)+blockRisk;
            end
        end
        % SCOPE SURE is nested exactly inside N-SCOPE SURE, using r=0.
        % The N-SCOPE MINIMUM *estimated* risk must be no larger, but its
        % true realized MSE need not be lower (SURE selection can overfit).
        assert(all(isfinite(sureScores(:))), ...
            'Non-finite plug-in SURE encountered in candidate search.');
        [minScope,local]=min(sureScores(cand.scopeIndices,:),[],1);
        idxScope=cand.scopeIndices(local);
        [minNScope,idxNScope]=min(sureScores,[],1);
        assert(all(minNScope<=minScope+1e-7.*max(1,abs(minScope))), ...
            'Nested SURE search consistency check failed.');
        state.scopeIdx(istart:iend,iF)=idxScope(:);
        state.nscopeIdx(istart:iend,iF)=idxNScope(:);
        % Four per-CDF methods in the stated output column order.
        choices={repmat(cand.fixedScope,1,nb),...
                 repmat(cand.fixedNScope,1,nb),idxScope,idxNScope};
        for meth=1:4
            loss=evaluate_selected_MSE(db,tr,S,errScaling,...
                  sigmahat,F,a,iF,cand,choices{meth},N);
            mseBatch(:,(iF-1)*4+meth)=loss(:);
        end
    end
    % Common universal-threshold references; same data, same sigmahat.
    hardLoss=errScaling; softLoss=errScaling;
    threshold=sigmahat*sqrt(2*log(N));
    for jj=1:7
        dj=db{jj}; tar=tr{jj};
        hard=bsxfun(@times,dj,bsxfun(@gt,abs(dj),threshold));
        soft=sign(dj).*max(bsxfun(@minus,abs(dj),threshold),0);
        hardLoss=hardLoss+sum(bsxfun(@minus,hard,tar).^2,1);
        softLoss=softLoss+sum(bsxfun(@minus,soft,tar).^2,1);
    end
    mseBatch(:,end-1)=hardLoss(:)/N;
    mseBatch(:,end)=softLoss(:)/N;
    assert(all(isfinite(mseBatch(:)) & mseBatch(:)>=0));
    state.mse(istart:iend,:)=mseBatch;
    state.completed=iend;
    state.rngState=rng;
    save(partialFn,'state','-v7');
    fprintf('   %s / SNR %g: %4d/%4d (%5.1f%%)\n',...
        tag,snr,iend,M,100*iend/M);
end
assert(state.completed==M && all(isfinite(state.mse(:))));
result=struct('signal',tag,'snr',snr,'filter',filterName,'M',M,...
   'names',{names},'mse',state.mse,'sigmahat',state.sigmahat,...
   'scopeIdx',state.scopeIdx,'nscopeIdx',state.nscopeIdx,...
   'cand',cand,'firstReplication',state.firstReplication,...
   'waveletShift',cfg.shift,'studySettings',study_settings(cfg));
end

function loss=evaluate_selected_MSE(db,tr,S,errScaling,sigmahat,F,a,iF,cand,idx,N)
% Squared error computed only for the *selected* candidate. Crucially, the
% TRUE coefficient tr is never accessed during SURE candidate selection.
nb=numel(idx); loss=errScaling;
for jj=1:7
    dj=db{jj}; target=tr{jj}; theta=zeros(size(dj));
    % Group replications with the same candidate; vectorize within group.
    distinct=unique(idx);
    for ii=1:numel(distinct)
        c=distinct(ii); cols=find(idx==c);
        ir=cand.ir(c); ik=cand.ik(c);
        lam=cand.b(c)*a(iF,ik,ir,jj)./sigmahat(cols);
        u=bsxfun(@times,S{jj,ir}(:,cols),lam);
        theta(:,cols)=dj(:,cols).*cdf_centered(F,u).^cand.k(c);
    end
    loss=loss+sum(bsxfun(@minus,theta,target).^2,1);
end
loss=loss/N;
end

% ========================================================================
function names=method_names(Fs)
nF=numel(Fs); names=cell(nF*4+2,1);
suffix={'FixedSCOPE','FixedNSCOPE','SCOPE_SURE','NSCOPE_SURE'};
for f=1:nF
 for meth=1:4
    names{(f-1)*4+meth}=sprintf('%s_%s',Fs{f},suffix{meth});
 end
end
names{end-1}='UniversalHard'; names{end}='UniversalSoft';
end

function write_mse_csv(R,fileName)
T=table((1:R.M)',R.sigmahat,...
     'VariableNames',{'Replication','SigmaHat'});
for c=1:numel(R.names)
    T.(R.names{c})=R.mse(:,c);
end
writetable(T,fileName);
end

function write_selection_csv(R,fileName)
T=table((1:R.M)',R.sigmahat,'VariableNames',{'Replication','SigmaHat'});
for fi=1:size(R.scopeIdx,2)
    % Derive generator tag from fixed-method name; no filename collisions.
    fullName=R.names{4*(fi-1)+1};
    tag=fullName(1:end-numel('_FixedSCOPE'));
    for itype=1:2
        if itype==1, ids=R.scopeIdx(:,fi); method='SCOPE';
        else, ids=R.nscopeIdx(:,fi); method='NSCOPE'; end
        T.([tag '_' method '_k'])=R.cand.k(ids);
        T.([tag '_' method '_r'])=R.cand.r(ids);
        T.([tag '_' method '_b'])=R.cand.b(ids);
    end
end
writetable(T,fileName);
end

% ========================================================================
function [overall,paired,selFreq,tripleFreq,ranking]=compile_tables(results,cfg,names)
methodTypes={'FixedSCOPE','FixedNSCOPE','SCOPE_SURE','NSCOPE_SURE'};
Overall={}; Paired={}; Sel={}; Triple={};
nF=numel(cfg.generators); nM=numel(names);
relative=zeros(numel(results),nM); casenum=0;
for is=1:size(results,1)
 for iq=1:size(results,2)
    R=results{is,iq}; M=R.M;
    means=mean(R.mse,1); sd=std(R.mse,0,1);
    casenum=casenum+1;
    relative(casenum,:)=means/means(end-1); % common HARD benchmark
    for fi=1:nF
        ind=(fi-1)*4+(1:4);
        for meth=1:4
            ii=ind(meth);
            gain=100*(means(ind(1))-means(ii))/means(ind(1));
            Overall(end+1,:)={R.signal,R.snr,R.filter,cfg.generators{fi},...
                 methodTypes{meth},R.names{ii},M,means(ii),sd(ii),...
                 sd(ii)/sqrt(M),gain}; %#ok<AGROW>
        end
        comparisons=[4 3;3 1;4 2;4 1;2 1];
        for pair=1:size(comparisons,1)
            iA=ind(comparisons(pair,1)); iB=ind(comparisons(pair,2));
            % Delta>0 indicates first method A has LOWER AMSE than B.
            delta=R.mse(:,iB)-R.mse(:,iA);
            dmean=mean(delta); se=std(delta)/sqrt(M);
            Paired(end+1,:)={R.signal,R.snr,cfg.generators{fi},...
                R.names{iA},R.names{iB},dmean,se,dmean-1.96*se,...
                dmean+1.96*se,100*dmean/means(iB)}; %#ok<AGROW>
        end
        for kind=1:2
            if kind==1, ids=R.scopeIdx(:,fi); mt='SCOPE_SURE';
            else, ids=R.nscopeIdx(:,fi); mt='NSCOPE_SURE'; end
            for which=1:3
                switch which
                 case 1, variable='k'; vals=cfg.kGrid; cv=R.cand.k(ids);
                 case 2, variable='r'; vals=cfg.rGrid; cv=R.cand.r(ids);
                 otherwise, variable='b'; vals=cfg.bGrid; cv=R.cand.b(ids);
                end
                for v=vals
                    n=sum(abs(cv-v)<1e-9);
                    Sel(end+1,:)={R.signal,R.snr,cfg.generators{fi},mt,...
                        variable,v,n,n/M}; %#ok<AGROW>
                end
            end
            for c=1:numel(R.cand.k)
                if kind==1 && R.cand.r(c)~=0, continue; end
                n=sum(ids==c);
                if n>0
                    Triple(end+1,:)={R.signal,R.snr,cfg.generators{fi},mt,...
                       R.cand.k(c),R.cand.r(c),R.cand.b(c),n,n/M}; %#ok<AGROW>
                end
            end
        end
    end
    for ref=1:2
       ii=nF*4+ref;
       Overall(end+1,:)={R.signal,R.snr,R.filter,'Universal',...
           'Reference',R.names{ii},M,means(ii),sd(ii),sd(ii)/sqrt(M),NaN}; %#ok<AGROW>
    end
 end
end
overall=cell2table(Overall,'VariableNames',{'Signal','SNR','Filter',...
 'Generator','Method','Procedure','M','AMSE','SD_MSE','SE_AMSE','GainPct_vs_FixedSCOPE'});
paired=cell2table(Paired,'VariableNames',{'Signal','SNR','Generator',...
 'Procedure_A','Procedure_B','Delta_BminusA','PairedSE','CI95_Lower','CI95_Upper',...
 'GainPct_A_over_B'});
selFreq=cell2table(Sel,'VariableNames',{'Signal','SNR','Generator',...
 'Method','Parameter','Value','Count','Frequency'});
tripleFreq=cell2table(Triple,'VariableNames',{'Signal','SNR','Generator',...
 'Method','K','R','B','Count','Frequency'});
meanRelative=mean(relative,1)';
[~,order]=sort(meanRelative,'ascend');
ranking=table(names(order),meanRelative(order),...
    'VariableNames',{'Procedure','MeanRelativeAMSE_to_UniversalHard'});
end

% ========================================================================
function condition_plots(R,cfg,plotDir)
nF=numel(cfg.generators); mu=mean(R.mse,1);
bars=zeros(nF,4); gain=zeros(nF,3); rFreq=zeros(nF,numel(cfg.rGrid));
for fi=1:nF
    bars(fi,:)=mu((fi-1)*4+(1:4));
    gain(fi,:)=100*(bars(fi,1)-bars(fi,2:4))/bars(fi,1);
    idx=R.nscopeIdx(:,fi);
    for ir=1:numel(cfg.rGrid)
        rFreq(fi,ir)=mean(R.cand.r(idx)==cfg.rGrid(ir));
    end
end
fig=figure('Visible','off','Color','w','Position',[100 100 1200 640]);
bar(bars,'grouped');
set(gca,'XTick',1:nF,'XTickLabel',cfg.generators,...
    'FontSize',11,'XTickLabelRotation',18);
ylabel('Average MSE');
title(sprintf('%s: fixed and SURE SCOPE/N-SCOPE, SNR %g, M=%d',...
   R.signal,R.snr,R.M));
legend({'SCOPE Fixed','N-SCOPE Fixed r=1','SCOPE SURE','N-SCOPE SURE'},...
   'Location','northoutside','Orientation','horizontal');
grid on;
base=fullfile(plotDir,sprintf('E2_AMSE_%s_SNR%d',R.signal,R.snr));
export_figure(fig,base);

fig=figure('Visible','off','Color','w','Position',[100 100 1200 640]);
bar(gain,'grouped');
set(gca,'XTick',1:nF,'XTickLabel',cfg.generators,...
    'FontSize',11,'XTickLabelRotation',18);
ylabel('AMSE reduction vs fixed SCOPE for same CDF (%)');
title(sprintf('%s: adaptive and neighborhood gains, SNR %g',R.signal,R.snr));
legend({'Fixed N-SCOPE r=1','SCOPE SURE','N-SCOPE SURE'},...
    'Location','northoutside','Orientation','horizontal');
grid on; hold on; lim=xlim; plot(lim,[0 0],'k-'); hold off;
base=fullfile(plotDir,sprintf('E2_Gains_%s_SNR%d',R.signal,R.snr));
export_figure(fig,base);

fig=figure('Visible','off','Color','w','Position',[100 100 1200 640]);
bar(100*rFreq,'stacked');
set(gca,'XTick',1:nF,'XTickLabel',cfg.generators,...
   'FontSize',11,'XTickLabelRotation',18);
ylabel('N-SCOPE SURE selected radius frequency (%)');
title(sprintf('%s: SURE radius selection, SNR %g, M=%d',R.signal,R.snr,R.M));
legend({'r=0','r=1','r=2','r=3'},...
    'Location','northoutside','Orientation','horizontal');
ylim([0 100]); grid on;
base=fullfile(plotDir,sprintf('E2_RadiusSelection_%s_SNR%d',R.signal,R.snr));
export_figure(fig,base);
end

function aggregate_plot(ranking,cfg,plotDir)
nF=numel(cfg.generators); vals=zeros(nF,4);
for fi=1:nF
 for meth=1:4
    nm=method_names(cfg.generators(fi));
    id=strcmp(ranking.Procedure,nm{meth});
    assert(sum(id)==1);
    vals(fi,meth)=ranking.MeanRelativeAMSE_to_UniversalHard(id);
 end
end
fig=figure('Visible','off','Color','w','Position',[100 100 1150 640]);
bar(vals,'grouped');
set(gca,'XTick',1:nF,'XTickLabel',cfg.generators,...
    'FontSize',11,'XTickLabelRotation',18);
ylabel('Average (AMSE / universal hard AMSE) over all signal-SNR cases');
title('Experiment II: aggregate non-oracle performance');
legend({'Fixed SCOPE','Fixed N-SCOPE','SCOPE SURE','N-SCOPE SURE'},...
    'Location','northoutside','Orientation','horizontal');
grid on;
export_figure(fig,fullfile(plotDir,'E2_AggregateRelativeAMSE'));
end

function export_figure(fig,basename)
print(fig,[basename '.png'],'-dpng','-r180');
print(fig,[basename '.pdf'],'-dpdf','-painters');
close(fig);
end

function write_results_readme(outDir,cfg)
fid=fopen(fullfile(outDir,'Experiment2_README.txt'),'w');
assert(fid>0);
fprintf(fid,'SCOPE Experiment II, October 2026 | n=1024, M=%d, sigma=1\n',cfg.M);
fprintf(fid,'DJ signals %s | power SNR %s\n',strjoin(cfg.signals,', '),mat2str(cfg.snrs));
fprintf(fid,'CDFs %s\n',strjoin(cfg.generators,', '));
fprintf(fid,'kGrid=%s; rGrid=%s; bGrid=%s\n',...
    mat2str(cfg.kGrid),mat2str(cfg.rGrid),mat2str(cfg.bGrid));
fprintf(fid,'Global tuning across all 7 shrunk detail levels j=3,...,9.\n');
fprintf(fid,'lambda_j=b*a_F,r,k(eta_j)/sigmahat, eta_j=.01*2^(9-j).\n');
fprintf(fid,'Gaussian null-risk constants via chi-square integral and fzero.\n');
fprintf(fid,'SURE uses the diagonal derivative holding sigmahat fixed.\n');
fprintf(fid,'Fixed k=5,b=1,r=0 and r=1 are matched reference estimators.\n');
fprintf(fid,'Error is measured after SURE selection, using wavelet-domain Parseval.\n');
fprintf(fid,'Same wavelet filters/signals and signal-SNR seeds as Experiment I.\n\n');
fprintf(fid,'Experiment2_AMSE_all.csv: mean MSE, SD, SE, gains for 4x6 methods.\n');
fprintf(fid,'Experiment2_PairedDifferences.csv: paired comparisons and 95%% MC CIs.\n');
fprintf(fid,'Experiment2_SelectionFrequencies.csv: k,r,b marginal frequencies.\n');
fprintf(fid,'Experiment2_SelectedTriples.csv: selected k,r,b joint frequencies.\n');
fprintf(fid,'Experiment2_ScaleConstants.csv: all F,k,r,level calibrated constants.\n');
fprintf(fid,'Experiment2_AggregateRanking.csv: 12-case offline relative AMSE.\n');
fprintf(fid,'MSE_replications_*.csv: paired per-run MSE, all methods.\n');
fprintf(fid,'SelectedParams_*.csv: per-replication selected k,r,b for both SURE.\n');
fprintf(fid,'checkpoint_*.mat: complete condition, PARTIAL_*.mat: resume.\n');
fprintf(fid,'plots/*: AMSE grouped bars, gains, radius frequencies, aggregate.\n');
fprintf(fid,'Experiment2_Results.mat: complete results and configuration.\n\n');
fprintf(fid,'IMPORTANT: plug-in SURE with estimated sigma is not exactly unbiased.\n');
fprintf(fid,'No candidate is optimized on clean signal; no oracle tuning.\n');
fprintf(fid,'Do not use these same replications as final independent Experiment III.\n');
fclose(fid);
end
