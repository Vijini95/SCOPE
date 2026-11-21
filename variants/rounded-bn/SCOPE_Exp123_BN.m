function SCOPE_Exp123_BN(mode, experiment, Ns_input, snrs_input, M_input)
% SCOPE_Exp123_BN | Section 4 rounded-B_N calibration (2026-10-30).
% Experiments I, II, III and PNG/PDF figures. No external MATLAB functions.
% Only built-in MATLAB functionality (R2018b or newer); no toolboxes.
%
% SIMULATION: N=[512 1024 2048], ROOT SNR=std(f,1)/sigma=[3 5 7],
% sigma=1, M=1000 paired replications per condition, four DJ signals.
% IMPORTANT: root SNR is an SD RATIO, NOT E(f^2)/sigma^2.
% Fixed Gaussian null-risk SCOPE with direct B_N/2^j targets, plug-in SURE, independent E3 benchmark.
% CALIBRATION: B_512=3.5, B_1024=5, B_2048=7.5 (also 128 and 4096).
% The original generator conventions, RNG seeds, SNR definition, filters,
% MSE and SURE choices are otherwise unchanged.
%
% EXPERIMENT III: exactly 26 methods, in the established table order.
% #5 BAMS reproduces the user's original author MATLAB settings:
%     Cscale=2, gamma=2, depth L=J-3, tau_j=tau_u*sqrt(j).
% #23 Chang--Yu--Vetterli context: 1D orthonormal DWT adaptation,
%     four lateral neighbours plus parent, LS context, grouped soft threshold.
% BAMSO (an alternative variant) and old Context_1D_proxy are NOT included.
% Some other published-method-inspired 1D adaptations remain provisional;
% see METHOD_IMPLEMENTATION_NOTES.csv produced by Experiment III.
%
% Original results and checkpoints are NEVER overwritten or migrated:
% output goes to SCOPE_Exp123_results_BNrounded. New SCOPE results are
% recomputed using the same seeded noise streams.
%
% EXAMPLES:
%   SCOPE_Exp123_BN('smoke','all')
%   SCOPE_Exp123_BN('full','Exp1',1024,[3 5 7],1000)
%   SCOPE_Exp123_BN('full','Exp2',1024,[3 5 7],1000)
%   SCOPE_Exp123_BN('full','Exp3',1024,[3 5 7],1000)
%   SCOPE_Exp123_BN('full','all',[512 1024 2048],[3 5 7],1000)
% Defaults: SCOPE_Exp123_BN('full','all') uses the complete design.
% All figures are written to plots subdirectories as .png and .pdf.
%
if nargin<1 || isempty(mode), mode='full'; end
if nargin<2 || isempty(experiment), experiment='all'; end
if nargin<3, Ns_input=[]; end
if nargin<4, snrs_input=[]; end
if nargin<5, M_input=[]; end
if isstring(mode), mode=char(mode); end
if isstring(experiment), experiment=char(experiment); end
assert(any(strcmpi(mode,{'full','smoke'})),'Mode must be full/smoke');
assert(any(strcmpi(experiment,{'all','Exp1','Exp2','Exp3'})), ...
 'Experiment must be all/Exp1/Exp2/Exp3');
if strcmpi(mode,'smoke')
 Ns=512; M=2; snrs=5;
else
 Ns=[512 1024 2048]; M=1000; snrs=[3 5 7];
 if ~isempty(Ns_input), Ns=Ns_input; end
 if ~isempty(snrs_input), snrs=snrs_input; end
 if ~isempty(M_input), M=M_input; end
end
Ns=unique(Ns(:).','stable');snrs=unique(snrs(:).','stable');
assert(isnumeric(Ns) && all(isfinite(Ns)) && ...
 all(Ns>=16 & log2(Ns)==floor(log2(Ns))), ...
 'Ns must contain powers of two >= 16.');
assert(isnumeric(snrs) && ~isempty(snrs) && all(isfinite(snrs) & snrs>0),...
 'SNRs must be positive finite numbers.');
assert(isnumeric(M) && isscalar(M) && isfinite(M) && M>=1 && M==floor(M),...
 'M must be a positive integer.');
% Isolate single-SNR runs in separate result directories; no CSV overwrite.
if isequal(snrs,[3 5 7])
 snrTag='';
else
 snrTag=['_SNR' strjoin(arrayfun(@(z)strrep(num2str(z,'%.8g'),'.','p'),...
                  snrs,'UniformOutput',false),'_')];
end
outBase=fullfile(fileparts(mfilename('fullpath')),'SCOPE_Exp123_results_BNrounded');
if ~exist(outBase,'dir'),mkdir(outBase);end
logLabel=sprintf('%s_N%s_M%d%s',upper(experiment),...
 strrep(strrep(mat2str(Ns),' ', '_'),'[',''),M,snrTag);
logLabel=strrep(logLabel,']','');
logPath=fullfile(outBase,['RUNLOG_' logLabel '.txt']);
diary(logPath); diary on;
diaryGuard=onCleanup(@()diary('off')); %#ok<NASGU>
fprintf('SAFE run: experiment=%s, N=%s, SNR=%s, M=%d.\n',...
 experiment,mat2str(Ns),mat2str(snrs),M);
fprintf('Log: %s\n',logPath);
fprintf('Noise: N(0,1); NOTE this code uses ROOT SNR = std(signal)/1.\n');
try
for N=Ns
 fprintf('\n================ N = %d ================\n',N);
 if any(strcmpi(experiment,{'all','Exp1'}))
   fprintf('BEGIN Exp1, N=%d\n',N);
   e1_run_SCOPE_Experiment1(mode,N,M,outBase,snrs,snrTag);
   fprintf('END Exp1, N=%d\n',N);
 end
 if any(strcmpi(experiment,{'all','Exp2'}))
   fprintf('BEGIN Exp2, N=%d\n',N);
   e2_run_SCOPE_Experiment2(mode,N,M,outBase,snrs,snrTag);
   fprintf('END Exp2, N=%d\n',N);
 end
 if any(strcmpi(experiment,{'all','Exp3'}))
   fprintf('BEGIN Exp3, N=%d\n',N);
   e3_run_SCOPE_Experiment3(mode,N,M,outBase,snrs,snrTag);
   fprintf('END Exp3, N=%d\n',N);
 end
end
if strcmpi(mode,'full')
 scope_consolidate_N(outBase,Ns,M,experiment,snrTag);
  if any(strcmpi(experiment,{'all','Exp1'}))
   scope_plot_N_SNR_sensitivity(outBase,Ns,snrs,M,snrTag);
  end
end
fprintf('\nCompleted. Results folder: %s\n',outBase);
catch ME
 fprintf(2,'\nMATLAB ERROR (not native crash):\n%s\n', ...
         getReport(ME,'extended','hyperlinks','off'));
 fprintf(2,'Run log: %s\n',logPath);
 rethrow(ME);
end
end
 
function [B_N,j0]=scope_B_for_N(N)
% Section 4: prespecified retained standardized null-noise energy PER LEVEL.
% N=2^J; each detail level has 2^j coefficients. Directly target
% R_F,r,k(a_{j})=B_N/(2^j). NOT an empirical-risk or oracle calibration.
switch N
 case 128,  B_N=1.75; j0=3;
 case 512,  B_N=3.5;  j0=3;
 case 1024, B_N=5;    j0=3;
 case 2048, B_N=7.5;  j0=3;
 case 4096, B_N=11.5; j0=4;
 otherwise
  error('Section 4 rounded B_N is specified only for N = 128,512,1024,2048,4096 (received %g).',N);
end
assert(B_N>0 && B_N<2^j0 && 2*max(0:3)+1<=2^j0, ...
 'Invalid null-risk budget or coarsest shrunk wavelet level.');
fprintf('Section 4 calibration: N=%d, B_N=%g, j0=%d, finest target=%g.\n', ...
        N,B_N,j0,B_N/(N/2));
end

function f=scope_scale_to_SNR(f0,SNR)
% Preserve sampled signal mean and waveform; population SD uses divisor N.
f0=f0(:); assert(isfinite(SNR)&&SNR>0 && std(f0,1)>0);
f=f0*(SNR/std(f0,1));
assert(abs(std(f,1)-SNR)<1e-9*max(1,SNR));
end
 
function h=MakeONFilterExt(Type,Par)
% Exact tap coefficients from the supplied MakeONFilterExt, not toolbox calls.
switch lower(Type)
 case 'haar'
  h=[1 1]/sqrt(2);
 case 'daubechies'
  assert(Par==6,'Use six-tap Daubechies-3 for Bumps');
  h=[0.332670552950082615998512 0.806891509311092576494494 ...
     0.459877502118491570095152 -0.135011020010254588696390 ...
    -0.085441273882026661692819 0.0352262918857095366027407];
 case 'symmlet'
  assert(Par==4,'Use Symmlet-4 for Doppler/HeaviSine');
  h=[-0.07576571478950221 -0.029635527646002493 ...
   0.497618667632775 0.8037387518051321 ...
   0.29785779560530606 -0.09921954357663353 ...
   -0.012603967262031304 0.032223100604051466];
 otherwise, error('Unrecognized filter %s',Type);
end
h=h/norm(h);
end
 
function W=Wavmat(h,N,k0,shift)
% SPARSE equivalent of the user's Wavmat(h,N,k0,shift) with shift=2.
% Output coefficient ordering is [2^j0 scaling, j0 detail, ... J-1 detail].
% Equivalent to computing hmat(ii,jj)=h(mod(N+ii-2jj+shift,nin)),
% gmat(ii,jj)=g(same index) and multiplying [oldmat*hmat';gmat'].
if nargin<4,shift=2;end
J=log2(N); assert(J==floor(J) && k0>=1 && k0<=J);
h=h(:)'; g=fliplr(conj(h).*(-1).^(1:numel(h)));
old=speye(2^(J-k0));
for k=k0:-1:1
  nout=2^(J-k); nin=2*nout; nh=numel(h);
  rows=repmat((1:nout)',nh,1);
  taps=kron((1:nh)',ones(nout,1));
  cols=mod(taps+2*rows-shift-1,nin)+1;
  H=sparse(rows,cols,repelem(h(:),nout),nout,nin);
  G=sparse(rows,cols,repelem(g(:),nout),nout,nin);
  old=[old*H;G];
end
W=old;
end
 
function z=scope_BAMS_author_coefficients(d,J,L,Cscale,gam)
% BAMS(20261028) author-defined variant in wavelet coordinates.
% Reuses the same W*y=d; avoids BAMS's second wavelet matrix multiplication.
% Original signal-domain BAMS returns W'*swt; applying W returns swt.
sw=d(:);n=numel(sw);coarsest=J-L;
assert(2^J==n && L>=1 && coarsest>=0 && coarsest<J, ...
    'BAMS requires a valid decomposition DEPTH L, not the coarsest level');
swt=sw; vals=sort(sw(n/2+1:n)); mm=numel(vals);
lo=floor(mm/4); hi=floor(3*mm/4);
pseudos=abs(vals(hi)-vals(lo))/Cscale;
mu=1/max(pseudos^2,realmin);
tauu=sqrt(max(std(sw).^2-1/mu,0.0001));
for j=J-1:-1:coarsest
 epsj=1-1/(j-coarsest+1)^gam;
 tauj=tauu*sqrt(j);ind=2^j+1:2^(j+1);
 for ii=ind
   [~,swt(ii)]=bams_bayesrule(sw(ii),mu,tauj,epsj);
 end
end
z=swt;
end
 
function scope_consolidate_N(outBase,Ns,M,experiment,snrTag)
% Combined files group all sizes; do not pool AMSE over different N.
for ie=1:3
 tag=sprintf('Exp%d',ie);
 if ~any(strcmpi(experiment,{'all',tag})),continue;end
 allT=table();
 for N=Ns
   folder=fullfile(outBase,sprintf('%s_N%d_M%d%s',tag,N,M,snrTag));
   file=fullfile(folder,sprintf('Experiment%d_AMSE_all.csv',ie));
   if exist(file,'file')
     tab=readtable(file);if isempty(allT), allT=tab; else, allT=[allT;tab]; end %#ok<AGROW>
   end
 end
 if ~isempty(allT)
   writetable(allT,fullfile(outBase,sprintf('AllN_%s_AMSE_RMSE%s.csv',tag,snrTag)));
 end
end
end
 
 
function scope_plot_N_SNR_sensitivity(outBase,Ns,snrs,M,snrTag)
% Main-paper sensitivity panel: fixed Normal generator r=1 vs r=0.
% Read actual simulation CSVs only; never fabricate or interpolate data.
if numel(Ns)<2 || numel(snrs)<2,return;end
signals={'Doppler','HeaviSine','Bumps','Blocks'};
gain=NaN(4,numel(Ns),numel(snrs));
for in=1:numel(Ns)
 folder=fullfile(outBase,sprintf('Exp1_N%d_M%d%s',Ns(in),M,snrTag));
 fn=fullfile(folder,'Experiment1_AMSE_all.csv');
 if ~exist(fn,'file'),return;end
 T=readtable(fn);
 for si=1:4
  for ki=1:numel(snrs)
   aa=strcmp(T.Signal,signals{si}) & T.SNR==snrs(ki) & ...
        strcmp(T.Procedure,'Normal_r0');
   bb=strcmp(T.Signal,signals{si}) & T.SNR==snrs(ki) & ...
        strcmp(T.Procedure,'Normal_r1');
   if sum(aa)~=1 || sum(bb)~=1,return;end
   gain(si,in,ki)=100*(T.AMSE(aa)-T.AMSE(bb))/T.AMSE(aa);
  end
 end
end
fig=figure('Visible','off','Color','w','Position',[50 50 1180 850]);
for si=1:4
 subplot(2,2,si);hold on;
 for ki=1:numel(snrs)
  plot(Ns,squeeze(gain(si,:,ki)),'-o','LineWidth',1.5,'MarkerSize',6);
 end
 plot(xlim,[0 0],'k:','LineWidth',.8);
 set(gca,'XTick',Ns);grid on;box on;
 xlabel('Sample size N');ylabel('AMSE reduction (%)');
 title(signals{si});
 if si==1
  leg=arrayfun(@(x)sprintf('SNR %g',x),snrs,'UniformOutput',false);
  legend(leg,'Location','best');
 end
end
plotDir=fullfile(outBase,'plots');if ~exist(plotDir,'dir'),mkdir(plotDir);end
base=fullfile(plotDir,sprintf('E1_Figure2_Sensitivity_N_SNR_M%d%s',M,snrTag));
print(fig,[base '.png'],'-dpng','-r180');
print(fig,[base '.pdf'],'-dpdf','-painters');close(fig);
end

function e1_run_SCOPE_Experiment1(mode,N_requested,M_requested,outputBase,snrsUser,snrTag)
% e1_run_SCOPE_Experiment1 -- Fixed SCOPE/N-SCOPE, Experiment I (October 2026)
%
% Completely self-contained MATLAB Monte Carlo using the supplied
% Wavmat.m and MakeONFilterExt.m (both included in this folder).
% NO Wavelet Toolbox, WaveLab, Statistics Toolbox or external files needed.
%
% Run the full experiment (4 signals x 3 SNRs x M=1000):
%     e1_run_SCOPE_Experiment1
% or a small executable check (2 signals x SNR 5 x M=12):
%     e1_run_SCOPE_Experiment1('smoke')
%
% SCIENTIFIC SPECIFICATION (matches Section 5 and Experiment I in Section 6)
%  * n=512,1024,2048; t_i=i/n; Gaussian noise variance one.
%  * SNR=std(f,1)/std(epsilon), the ROOT (SD-ratio) SNR; SNR={3,5,7}.
%  * Gaussian epsilon ~ N(0,1), std(epsilon)=1; scale f to std(f,1)=SNR.
%  * Population SD std(f,1) uses divisor n; do NOT rescale noise draws.
%  * Symmlet-4 (8 taps): Doppler, HeaviSine
%    Daubechies length-6: Bumps; Haar: Blocks.
%  * Periodic orthonormal DWT, 8 scaling coefficients UNTOUCHED,
%    levels j=3:(J-1) shrunk, finest level has N/2 coefficients.
%  * Six centered CDFs: Normal, Logistic, Laplace, hyperbolic secant,
%    Uniform, Cauchy. Fixed k=5, radii r=0:3.
%  * Section 4 null-risk target at each level j: B_N/2^j.
%  * Solve the EXACT one-dimensional chi-square integral for a(F,r,j).
%    lambda_j(F,r) = a(F,r,j) / sigmahat, separately for each radius.
%  * sigmahat=median(abs(finest detail coefficients))/0.6744897501960817.
%  * No oracle parameter fitting, no SURE, no CV in Experiment I.
%  * Paired comparisons: identical noise for ALL 24 fixed procedures.
%  * MSE computed IN THE ORTHONORMAL WAVELET DOMAIN, equivalent by
%    Parseval's identity to signal-domain reconstruction MSE.
%  * Universal hard and soft rules included as reference methods.
%
% OUTFILES (created alongside this file in outputs_E1_SNRstd_M1000):
%  Experiment1_AMSE_all.csv; Experiment1_PairedGains.csv;
%  Experiment1_AggregateRanking.csv; Experiment1_Constants.csv;
%  MSE_replications_*.csv; checkpoint_*.mat; PNG/PDF bar plots;
%  Experiment1_Results.mat; Experiment1_README.txt.
%
% This function is intended as the reproducible production program.
% The demonstration ('smoke') run does NOT replace M=1000 evidence.
 
if nargin < 1 || isempty(mode), mode = 'full'; end
if ~ischar(mode), mode=char(mode); end
smoke = strcmpi(mode,'smoke');
if ~(smoke || strcmpi(mode,'full'))
    error('Mode must be ''full'' or ''smoke''.');
end
 
rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
 
cfg.N = N_requested;
cfg.J = round(log2(cfg.N));
[cfg.B_N,cfg.j0] = scope_B_for_N(cfg.N);
cfg.k = 5;
cfg.radii = 0:3;
cfg.generators = {'Normal','Logistic','Laplace','Sech','Uniform','Cauchy'};
cfg.signals = {'Doppler','HeaviSine','Bumps','Blocks'};
cfg.snrs = snrsUser;
cfg.snrDefinition = 'std_signal_over_std_noise_iidN01';
cfg.M = M_requested;
cfg.batchSize = 100;
cfg.baseSeed = 20261026 + cfg.N;
cfg.makePlots = true;
cfg.resume = true;
cfg.shift = 2;  % original supplied Wavmat periodic phase convention
 
if smoke
    cfg.signals = {'Doppler','Blocks'};
    cfg.snrs = 5;
    cfg.M = M_requested;
    cfg.batchSize = min(6,cfg.M);
    cfg.makePlots = true;
    outDir = fullfile(outputBase,sprintf('Exp1_N%d_SMOKE',cfg.N));
else
    outDir = fullfile(outputBase,sprintf('Exp1_N%d_M%d%s',cfg.N,cfg.M,snrTag));
end
if ~exist(outDir,'dir'), mkdir(outDir); end
plotDir = fullfile(outDir,'plots');
if ~exist(plotDir,'dir'), mkdir(plotDir); end
 
names = e1_method_names(cfg.generators,cfg.radii);
assert(numel(names)==26); % 6 CDFs x 4 radii + universal hard/soft
 
% Calibrate once; this depends on the noise model, never on the signal.
constPath = fullfile(outDir,'calibration_E1.mat');
if exist(constPath,'file')
    d = load(constPath,'a','cfgCal');
    cfgCal = struct('k',cfg.k,'B_N',cfg.B_N,'j0',cfg.j0,...
      'N',cfg.N,'generators',{cfg.generators},'radii',cfg.radii);
    assert(isfield(d,'cfgCal') && isequaln(d.cfgCal,cfgCal),...
       'Cached calibration metadata mismatch. Delete cache and rerun.');
    a = d.a;
else
    fprintf('Computing Gaussian null-risk calibration (6 CDF x 4 radii x %d levels)...\n',cfg.J-cfg.j0);
    a = e1_compute_constants(cfg);
    cfgCal = struct('k',cfg.k,'B_N',cfg.B_N,'j0',cfg.j0,...
      'N',cfg.N,'generators',{cfg.generators},'radii',cfg.radii);
    save(constPath,'a','cfgCal');
end
assert(isequal(size(a),[6 4 cfg.J-cfg.j0]),'Incorrect calibration array size.');
e1_write_constant_csv(cfg,a,fullfile(outDir,'Experiment1_Constants.csv'));
e1_check_constants(a,cfg); % diagnostic consistency check with current paper
 
% Build wavelet transform matrices once, not inside Monte Carlo loops.
% Wavmat order: [8 scaling, 8 details(j=3), ..., N/2 details(j=J-1)].
N = cfg.N;
Wbank = cell(1,3);
Wbank{1} = e1_wavelet_matrix('Symmlet',4,cfg);  % 8 taps
Wbank{2} = e1_wavelet_matrix('Daubechies',6,cfg); % 6 taps
Wbank{3} = e1_wavelet_matrix('Haar',0,cfg);
 
results = cell(numel(cfg.signals),numel(cfg.snrs));
tic;
for is = 1:numel(cfg.signals)
    tag = cfg.signals{is};
    [f0,filterTag,whichW] = e1_make_DJ_signal(tag,cfg.N);
    W = Wbank{whichW};
    fprintf('\n================ %s (%s) ================\n',tag,filterTag);
    for iq = 1:numel(cfg.snrs)
        snr = cfg.snrs(iq);
        checkpoint = fullfile(outDir,sprintf('checkpoint_%s_SNR%d.mat',tag,snr));
        if cfg.resume && exist(checkpoint,'file')
            loadStruct = load(checkpoint,'result');
            if isfield(loadStruct,'result') && ...
                    isfield(loadStruct.result,'snrDefinition') && ...
                    strcmp(loadStruct.result.snrDefinition,cfg.snrDefinition) && ...
                    size(loadStruct.result.mse,1)==cfg.M && ...
                    size(loadStruct.result.mse,2)==numel(names) && ...
                    isfield(loadStruct.result,'calibration_BN') && ...
                    loadStruct.result.calibration_BN==cfg.B_N && ...
                    isfield(loadStruct.result,'j0') && ...
                    loadStruct.result.j0==cfg.j0
                result = loadStruct.result;
                fprintf('RESUMED %s / SNR=%g from checkpoint.\n',tag,snr);
            else
                error('Checkpoint incompatible with current M or methods: %s',checkpoint);
            end
        else
            fprintf('RUNNING %s, SNR=std(f)/std(noise)=%g, M=%d...\n',tag,snr,cfg.M);
            f = scope_scale_to_SNR(f0,snr); % std(f,1)=SNR; epsilon~N(0,1)
            fprintf('  Verified signal SD = %.8f; population noise SD = 1.\n',std(f,1));
            rng(cfg.baseSeed + 100*is + iq,'twister'); % deterministic per condition
            partialFn = fullfile(outDir,sprintf('PARTIAL_%s_SNR%d.mat',tag,snr));
            result = e1_simulate_condition(f,W,a,cfg,names,tag,snr,filterTag,partialFn);
            save(checkpoint,'result','-v7'); % recover after interruption
            if exist(partialFn,'file'), delete(partialFn); end
        end
        results{is,iq} = result;
        e1_write_mse_csv(result,fullfile(outDir,sprintf('MSE_replications_%s_SNR%d.csv',tag,snr)));
        if cfg.makePlots
            e1_make_barplots(result,cfg,plotDir);
        end
        fprintf('FINISHED %s / SNR=%g (elapsed %.1f s).\n',tag,snr,toc);
    end
end
 
[amsetab,gainstab,ranktab] = e1_compile_tables(results,cfg,names);
amsetab.N = repmat(cfg.N,height(amsetab),1);
amsetab.RMSE = sqrt(amsetab.AMSE);
writetable(amsetab,fullfile(outDir,'Experiment1_AMSE_all.csv'));
writetable(gainstab,fullfile(outDir,'Experiment1_PairedGains.csv'));
writetable(ranktab,fullfile(outDir,'Experiment1_AggregateRanking.csv'));
save(fullfile(outDir,'Experiment1_Results.mat'),'cfg','a','results','-v7.3');
e1_write_readme(outDir,cfg);
 
fprintf('\nFULL EXPERIMENT FINISHED. Output: %s\n',outDir);
fprintf('Inspect Experiment1_AMSE_all.csv and Experiment1_PairedGains.csv.\n');
fprintf('Total elapsed: %.1f seconds.\n',toc);
end
 
% -------------------------------------------------------------------------
function W = e1_wavelet_matrix(filterName,par,cfg)
% Wavmat is the original matrix representation provided by the author.
% Calling with k0=J-j0 retains 2^j0=8 scaling coefficients.
h = MakeONFilterExt(filterName,par);
assert(numel(h)==8 || numel(h)==6 || numel(h)==2, ...
    'Unexpected wavelet filter length.');
W = Wavmat(h,cfg.N,cfg.J-cfg.j0,cfg.shift);
assert(isequal(size(W),[cfg.N cfg.N]),'Wavelet matrix dimension mismatch.');
v = sin((1:cfg.N)'*1.2345); % deterministic vector; no RNG consumed
err = norm(W'*(W*v)-v)/norm(v);
assert(err<3e-8,'Wavmat failed the invertibility/orthogonality check.');
fprintf('Wavelet: %-12s length=%d, reconstruction check %.3g\n',...
    filterName,numel(h),err);
end
 
% -------------------------------------------------------------------------
function [f0,label,whichW] = e1_make_DJ_signal(name,N)
% Classical Donoho-Johnstone definitions. Section 6 uses t_i=i/N.
t = (1:N)'/N;
pos = [0.10 0.13 0.15 0.23 0.25 0.40 0.44 0.65 0.76 0.78 0.81];
htB = [4.0 5.0 3.0 4.0 5.0 4.2 2.1 4.3 3.1 5.1 4.2];
htK = [4.0 -5.0 3.0 -4.0 5.0 -4.2 2.1 4.3 -3.1 2.1 -4.2];
w = [0.005 0.005 0.006 0.010 0.010 0.030 0.010 0.010 0.005 0.008 0.005];
switch lower(name)
    case 'doppler'
        f0 = sqrt(t.*(1-t)).*sin(2*pi*1.05./(t+0.05));
        label='Symmlet4_8tap'; whichW=1;
    case 'heavisine'
        f0 = 4*sin(4*pi*t)-sign(t-0.30)-sign(0.72-t);
        label='Symmlet4_8tap'; whichW=1;
    case 'bumps'
        f0 = zeros(size(t));
        for i=1:numel(pos)
            f0 = f0+htB(i)./(1+abs((t-pos(i))/w(i))).^4;
        end
        label='Daubechies3_6tap'; whichW=2;
    case 'blocks'
        f0 = zeros(size(t));
        for i=1:numel(pos)
            f0 = f0+htK(i)*(t>=pos(i));
        end
        label='Haar_2tap'; whichW=3;
    otherwise
        error('Unknown signal: %s',name);
end
assert(iscolumn(f0) && all(isfinite(f0)) && mean(f0.^2)>0);
end
 
% -------------------------------------------------------------------------
function a = e1_compute_constants(cfg)
% Exact null-risk calibration from Proposition 4 (Eq. 68):
%   U ~ chi-square_m, m=2*r+1.
%   R_F,r(a) = E[ (U/m)*Fstar(a*sqrt(U/m))^(2*k) ].
% For j=j0,...,J-1, solve R(a)=B_N/2^j (direct Section 4 rule).
% Crucially, the constants are SPECIFIC to r and F, not just quantiles.
nF=numel(cfg.generators); nR=numel(cfg.radii); nJ=cfg.J-cfg.j0;
a = zeros(nF,nR,nJ);
for iF=1:nF
    for ir=1:nR
        for jj=1:nJ
            j=cfg.j0+jj-1;
            target=cfg.B_N/2^j;
            fun=@(z) e1_nullrisk(z,cfg.generators{iF},cfg.radii(ir),cfg.k);
            upper=1;
            while fun(upper)<target
                upper=2*upper;
                assert(upper<1e6,'Failed to bracket null-risk equation.');
            end
            a(iF,ir,jj)=fzero(@(z)fun(z)-target,[0,upper],...
                optimset('TolX',1e-10,'Display','off'));
        end
        fprintf('  %-9s r=%d, a_j9=%.4f\n',cfg.generators{iF},cfg.radii(ir),a(iF,ir,end));
    end
end
end
 
% -------------------------------------------------------------------------
function val = e1_nullrisk(a,F,r,k)
% MATLAB-only chi-square expectation, no Statistics Toolbox chi2pdf.
if a==0, val=0; return; end
m=2*r+1; nu=m/2;
logC=-nu*log(2)-gammaln(nu)-log(m);
% chi-square density times U/m. u=0 is treated by continuous limit.
fun=@(u) exp(nu.*log(max(u,realmin))-.5*u+logC) ...
    .* e1_cdf_centered(F,a.*sqrt(u/m)).^(2*k);
val=integral(fun,0,Inf,'AbsTol',2e-10,'RelTol',1e-8);
val=min(1,max(0,val));
end
 
% -------------------------------------------------------------------------
function z = e1_cdf_centered(F,u)
% F*(u)=2F(u)-1 for u>=0, with EXACT paper parameter conventions.
u=max(0,u); % always local RMS argument; nonnegative
switch lower(F)
    case 'normal'
        z=erf(u/sqrt(2));               % 2*Phi(u)-1
    case 'logistic'
        z=tanh(u);                      % manuscript absorbs factor 1/2
    case 'laplace'
        z=-expm1(-u);                   % 1-exp(-u)
    case 'sech'
        z=(2/pi).*atan(sinh(pi*u/2));   % Hyperbolic-secant CDF
    case 'uniform'
        z=min(u,1);
    case 'cauchy'
        z=(2/pi).*atan(u);
    otherwise
        error('Unrecognized centered CDF: %s',F);
end
z=max(0,min(1,z));
end
 
% -------------------------------------------------------------------------
function e1_check_constants(a,cfg)
% Check actual finest-level Gaussian null risks instead of a fixed N=1024 table.
for iF=1:numel(cfg.generators)
 for ir=1:numel(cfg.radii)
  r=e1_nullrisk(a(iF,ir,end),cfg.generators{iF},cfg.radii(ir),cfg.k);
  assert(abs(r-cfg.B_N/2^(cfg.J-1))<1e-5,'E1 calibration root failed');
 end
end
end
function result = e1_simulate_condition(f,W,a,cfg,names,tag,snr,filterTag,partialFn)
N=cfg.N; M=cfg.M; nMeth=numel(names);
trueD = W*f;
finest=(N/2+1):N; % j=9, indices 513..1024
% E1 is now restartable within a signal/SNR condition, preserving RNG.
if cfg.resume && exist(partialFn,'file')
    tmp=load(partialFn,'state'); state=tmp.state;
    assert(state.N==N && state.M==M && ...
       isfield(state,'calibration_BN') && ...
       state.calibration_BN==cfg.B_N && ...
       isfield(state,'j0') && state.j0==cfg.j0 && ...
       strcmp(state.signal,tag) && state.snr==snr && ...
       isequal(size(state.mse),[M nMeth]), ...
       'Incompatible E1 partial checkpoint: %s',partialFn);
    mse=state.mse; sigmahats=state.sigmahats;
    repFirst=state.repFirst; completed=state.completed;
    rng(state.rngState);
    fprintf('  E1 resumed partial %s SNR %g at %d/%d\n',tag,snr,completed,M);
else
    mse=zeros(M,nMeth);sigmahats=zeros(M,1);
    repFirst=[];completed=0;
end
for istart=(completed+1):cfg.batchSize:M
    iend=min(M,istart+cfg.batchSize-1);
    nb=iend-istart+1;
    Y=bsxfun(@plus,f,randn(N,nb));
    Dall=W*Y; % one transform for ALL methods
    sigmahat=median(abs(Dall(finest,:)),1)/0.6744897501960817;
    assert(all(sigmahat>0) && all(isfinite(sigmahat)));
    sigmahats(istart:iend)=sigmahat';
    if istart==1
        repFirst=struct('y',Y(:,1),'D',Dall(:,1),'sigmahat',sigmahat(1));
    end
    % Contribution of unshrunk scaling coefficients common to all methods.
    errScaling=sum(bsxfun(@minus,Dall(1:2^cfg.j0,:),trueD(1:2^cfg.j0)).^2,1);
    loss=repmat(errScaling,nMeth,1); % [method x replication]
    for jj=1:(cfg.J-cfg.j0)
        j=cfg.j0+jj-1;
        i1=2^j+1; i2=2^(j+1);  % 1-based wavelet ordering
        dj=Dall(i1:i2,:);
        target=trueD(i1:i2);
        rmsvals=cell(1,4);
        rmsvals{1}=abs(dj);
        dsq=dj.^2;
        for ir=2:4
            r=cfg.radii(ir);
            e=dsq;
            for lag=1:r
                e=e+circshift(dsq,[lag 0])+circshift(dsq,[-lag 0]);
            end
            rmsvals{ir}=sqrt(e/(2*r+1));
        end
        for iF=1:numel(cfg.generators)
            for ir=1:4
                im=(iF-1)*4+ir;
                lam=a(iF,ir,jj)./sigmahat; % row vector, one value per MC trial
                u=bsxfun(@times,rmsvals{ir},lam);
                attenuation=e1_cdf_centered(cfg.generators{iF},u).^cfg.k;
                theta=dj.*attenuation;
                loss(im,:)=loss(im,:)+sum(bsxfun(@minus,theta,target).^2,1);
            end
        end
        % Prespecified universal hard and soft comparator references.
        threshold=sigmahat*sqrt(2*log(N));
        thrmat=repmat(threshold,size(dj,1),1);
        hard=dj.*(abs(dj)>thrmat);
        soft=sign(dj).*max(abs(dj)-thrmat,0);
        loss(25,:)=loss(25,:)+sum(bsxfun(@minus,hard,target).^2,1);
        loss(26,:)=loss(26,:)+sum(bsxfun(@minus,soft,target).^2,1);
    end
    mse(istart:iend,:)=loss.'/N;
    state=struct('N',N,'M',M,'signal',tag,'snr',snr,...
       'calibration_BN',cfg.B_N,'j0',cfg.j0,...
       'mse',mse,'sigmahats',sigmahats,'repFirst',repFirst,...
       'completed',iend,'rngState',rng);
    save(partialFn,'state','-v7');
    fprintf('  %s / SNR %g : %4d/%4d (%5.1f%%)\n',...
        tag,snr,iend,M,100*iend/M);
end
result=struct('signal',tag,'snr',snr,'filter',filterTag,'M',M,...
   'names',{names},'mse',mse,'sigmahat',sigmahats,'f',f,...
   'firstReplication',repFirst,'waveletShift',cfg.shift,...
   'snrDefinition',cfg.snrDefinition,'calibration_BN',cfg.B_N,'j0',cfg.j0);
end
 
% -------------------------------------------------------------------------
function names = e1_method_names(Fs,rs)
names=cell(numel(Fs)*numel(rs)+2,1);
i=0;
for f=1:numel(Fs)
    for ir=1:numel(rs)
        i=i+1;
        names{i}=sprintf('%s_r%d',Fs{f},rs(ir));
    end
end
names{end-1}='UniversalHard'; names{end}='UniversalSoft';
end
 
% -------------------------------------------------------------------------
function e1_write_constant_csv(cfg,a,filePath)
G={}; R=[]; J=[]; riskTarget=[]; budget=[]; aa=[];
for fi=1:numel(cfg.generators)
    for ir=1:numel(cfg.radii)
        for jj=1:(cfg.J-cfg.j0)
            j=cfg.j0+jj-1;
            G{end+1,1}=cfg.generators{fi}; %#ok<AGROW>
            R(end+1,1)=cfg.radii(ir); %#ok<AGROW>
            J(end+1,1)=j; %#ok<AGROW>
            riskTarget(end+1,1)=cfg.B_N/2^j; %#ok<AGROW>
            budget(end+1,1)=cfg.B_N; %#ok<AGROW>
            aa(end+1,1)=a(fi,ir,jj); %#ok<AGROW>
        end
    end
end
T=table(G,R,J,budget,riskTarget,aa,'VariableNames',...
    {'Generator','Radius','Level_j','Budget_BN','TargetNullRisk','a_NullRisk'});
writetable(T,filePath);
end
 
% -------------------------------------------------------------------------
function e1_write_mse_csv(R,path)
M=size(R.mse,1);
T=table((1:M)',R.sigmahat,'VariableNames',{'Replication','SigmaHat'});
for i=1:numel(R.names)
    T.(R.names{i})=R.mse(:,i);
end
writetable(T,path);
end
 
% -------------------------------------------------------------------------
function [overall,gains,ranking] = e1_compile_tables(results,cfg,names)
% One row per condition and procedure; paired comparisons on replication MSE.
Tcell={}; Gcell={}; scores=zeros(numel(results),numel(names));
caseidx=0;
for is=1:size(results,1)
    for iq=1:size(results,2)
        r=results{is,iq};
        caseidx=caseidx+1;
        m=mean(r.mse,1);
        sd=std(r.mse,0,1);
        M=r.M;
        scores(caseidx,:)=m./m(1); % shared Normal-r0 reference, per condition
        for i=1:numel(names)
            if i<=24
                fi=ceil(i/4); rad=mod(i-1,4);
                family=cfg.generators{fi};
                if rad==0, method='SCOPE_Fixed'; else, method='NSCOPE_Fixed'; end
                base=mean(r.mse(:,(fi-1)*4+1));
                gain=100*(base-m(i))/base;
            else
                fi=0; rad=-1; family='Universal'; method='Reference'; gain=NaN;
            end
            Tcell(end+1,:)={r.signal,r.snr,r.filter,names{i},method,...
                family,rad,M,m(i),sd(i),sd(i)/sqrt(M),gain}; %#ok<AGROW>
            if i<=24 && rad>0
                u=r.mse(:,(fi-1)*4+1)-r.mse(:,i);
                du=mean(u);
                se=std(u)/sqrt(M);
                Gcell(end+1,:)={r.signal,r.snr,family,rad,du,se,...
                    du-1.96*se,du+1.96*se,gain}; %#ok<AGROW>
            end
        end
    end
end
overall=cell2table(Tcell,'VariableNames',{'Signal','SNR','Filter','Procedure',...
    'Class','Generator','Radius','M','AMSE','SD_MSE','SE_AMSE','GainPct_vs_sameF_r0'});
gains=cell2table(Gcell,'VariableNames',{'Signal','SNR','Generator','Radius',...
    'PairedDeltaAMSE','PairedSE_Delta','CI95_Lower','CI95_Upper','GainPercent'});
% Descriptive OFFLINE ranking; not in-sample tuning of any reconstruction.
relative=mean(scores,1)';
[~,order]=sort(relative,'ascend');
ranking=table(names(order),relative(order),'VariableNames',...
    {'Procedure','MeanRelativeAMSE_to_Normal_r0'});
end
 
% -------------------------------------------------------------------------
function e1_make_barplots(R,cfg,plotDir)
% TWO manuscript-friendly bar figures for each signal and SNR:
% (1) AMSE: six CDF x {r0,r1,r2,r3}.
% (2) Percent gain vs SAME generator with r0 (positive = improvement).
mu=mean(R.mse,1);
labels=cfg.generators;
bars=zeros(6,4); gains=zeros(6,3);
for fi=1:6
    bars(fi,:)=mu((fi-1)*4+(1:4));
    gains(fi,:)=100*(bars(fi,1)-bars(fi,2:4))/bars(fi,1);
end
fig=figure('Visible','off','Color','w','Position',[100 100 1150 620]);
bar(bars,'grouped');
set(gca,'XTick',1:6,'XTickLabel',labels,'FontSize',11,'XTickLabelRotation',18);
ylabel('Average MSE (MSE per signal point)');
title(sprintf('%s: fixed SCOPE/N-SCOPE, SNR (SD ratio) %g, M=%d',R.signal,R.snr,R.M));
legend({'r=0 (SCOPE)','r=1','r=2','r=3'},'Location','northoutside',...
    'Orientation','horizontal');
grid on;
base=fullfile(plotDir,sprintf('E1_AMSE_%s_SNR%d',R.signal,R.snr));
print(fig,[base,'.png'],'-dpng','-r180');
print(fig,[base,'.pdf'],'-dpdf','-painters');
close(fig);
 
fig=figure('Visible','off','Color','w','Position',[100 100 1150 620]);
bar(gains,'grouped');
set(gca,'XTick',1:6,'XTickLabel',labels,'FontSize',11,'XTickLabelRotation',18);
ylabel('AMSE reduction vs same-CDF SCOPE r=0 (%)');
title(sprintf('%s: neighborhood gains, SNR (SD ratio) %g, M=%d',R.signal,R.snr,R.M));
legend({'r=1','r=2','r=3'},'Location','northoutside',...
    'Orientation','horizontal');
grid on; hold on; xl=xlim; plot(xl,[0 0],'k-','LineWidth',0.8); hold off;
base=fullfile(plotDir,sprintf('E1_GainPercent_%s_SNR%d',R.signal,R.snr));
print(fig,[base,'.png'],'-dpng','-r180');
print(fig,[base,'.pdf'],'-dpdf','-painters');
close(fig);
end
 
% -------------------------------------------------------------------------
function e1_write_readme(outDir,cfg)
fid=fopen(fullfile(outDir,'Experiment1_README.txt'),'w');
assert(fid>0,'Could not create output README.');
fprintf(fid,'SCOPE Experiment I | fixed k=%d B_N=%.5g | n=%d M=%d\n',...
    cfg.k,cfg.B_N,cfg.N,cfg.M);
fprintf(fid,'Donoho-Johnstone: Doppler, HeaviSine, Bumps, Blocks.\n');
fprintf(fid,'Time grid t_i=i/n; SNR=std(f,1)/std(epsilon)=%s; epsilon iid N(0,1).\n',mat2str(cfg.snrs));
fprintf(fid,'std(epsilon)=1 theoretically; f scaled to std(f,1)=SNR (divisor n).\n');
fprintf(fid,'Symmlet-4 (8 taps), Daubechies-3 (6 taps), Haar, periodic DWT.\n');
fprintf(fid,'SCOPE = r0; N-SCOPE = r1,r2,r3; 6 centered CDFs.\n');
fprintf(fid,'Detail levels j=%d,...,%d; %d unshrunk scaling coefficients.\n', ...
    cfg.j0,cfg.J-1,2^cfg.j0);
fprintf(fid,'Calibration: R_F,r,k(a_j)=B_N/2^j; B_N=%.8g; j0=%d; J=%d.\n', ...
    cfg.B_N,cfg.j0,cfg.J);
fprintf(fid,'Per replication sigmaHat=median(abs(finest (j=%d) details))/0.6744897501960817.\n', ...
    cfg.J-1);
fprintf(fid,'Methodwise MSE computed in wavelet domain via Parseval, exactly reconstruction MSE.\n');
fprintf(fid,'Within signal/SNR all methods use the identical simulated Gaussian noise.\n');
fprintf(fid,'Each signal/SNR uses a fixed, reproducible MT19937 seed.\n\n');
fprintf(fid,'FILES\n  Experiment1_AMSE_all.csv -- all 24 fixed rules + 2 references.\n');
fprintf(fid,'  Experiment1_PairedGains.csv -- paired r0 minus r1/r2/r3 differences.\n');
fprintf(fid,'  Experiment1_AggregateRanking.csv -- descriptive offline aggregate.\n');
fprintf(fid,'  Experiment1_Constants.csv -- all N-dependent calibrated scale constants.\n');
fprintf(fid,'  MSE_replications_*.csv -- per-replication MSE for each method.\n');
fprintf(fid,'  checkpoint_*.mat -- resumable condition-level results.\n');
fprintf(fid,'  plots/E1_AMSE_*.{png,pdf} -- bar plots by signal/SNR.\n');
fprintf(fid,'  plots/E1_GainPercent_*.{png,pdf} -- gains relative r=0.\n');
fprintf(fid,'  Experiment1_Results.mat -- complete numeric arrays / configuration.\n\n');
fprintf(fid,'IMPORTANT: Fixed calibration only. No SURE, no oracle.\n');
fprintf(fid,'The offline ranking must NOT be interpreted as per-realization tuning.\n');
fprintf(fid,'Validation: direct B_N/2^j finest null-risk target checked numerically.\n');
fclose(fid);
end
 
 
function e2_run_SCOPE_Experiment2(mode,N_requested,M_requested,outputBase,snrsUser,snrTag)
% RUN_SCOPE_EXPERIMENT2   SCOPE SURE and N-SCOPE SURE, Section 6 Experiment II.
%
% Full manuscript experiment (4 DJ signals x 3 SNRs x M=1000):
%       call SCOPE_Exp123(mode,'Exp2') from the top-level entry point
% Quick functionality check, NOT publishable:
%       call SCOPE_Exp123(mode,'Exp2') from the top-level entry point('smoke')
%
% SELF-CONTAINED: Needs only Wavmat.m and MakeONFilterExt.m included here.
% No MATLAB Wavelet Toolbox, Statistics Toolbox, WaveLab, or Optimization
% Toolbox required. MATLAB with local functions, tables, and export graphics
% via print() is sufficient (R2018a or newer recommended).
%
% DESIGN: n=512,1024,2048; t_i=i/n, iid Gaussian noise N(0,1).
% SNR=std(f,1)/std(epsilon)={3,5,7}; std(epsilon)=1, so std(f,1)=SNR.
% std(f,1) uses divisor n. Never normalize individual noise realizations.
% Periodic orthonormal DWT with 8 unshrunk scaling coefficients:
%  Doppler, HeaviSine: Symmlet-4 (8 taps)
%  Bumps:              Daubechies length 6
%  Blocks:             Haar
% Details j=3,...,J-1. SigmaHat=median(abs(finest details))/0.674489750196.
%
% SIX CDFs: Normal, Logistic, Laplace, Sech, Uniform, Cauchy.
% ADAPTIVE GRIDS: k={3,5,7}; r={0,1,2,3}; b={0.65,0.8,1,1.2,1.5}.
% lambda_{j,r}(F,k,b) = b*R^{-1}_{F,r,k}(B_N/2^j) / SigmaHat,
% where a solves the *exact* Gaussian null-risk
% equation with U ~ ChiSquare(2r+1):
%   E[(U/(2r+1)) * Fstar(a*sqrt(U/(2r+1)))^(2*k)] = B_N/2^j.
%
% For each CDF and replicate, SCOPE_SURE MINIMIZES the summed plug-in SURE
% over all (k,b) at r=0; NSCOPE_SURE minimizes over all (k,r,b).
% All tuning selections are GLOBAL across levels j=3:(J-1), not levelwise.
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
% OUTPUT: SCOPE_Exp123_results/Exp2_N<N>_M1000 (or smoke).
% Consult Experiment2_README.txt for generated CSV/MAT records.
 
if nargin<1 || isempty(mode), mode='full'; end
if ~ischar(mode), mode=char(mode); end
smoke=strcmpi(mode,'smoke');
assert(smoke || strcmpi(mode,'full'), 'Mode must be ''smoke'' or ''full''.');
rootDir=fileparts(mfilename('fullpath'));
addpath(rootDir);
 
cfg.N=N_requested; cfg.J=round(log2(cfg.N));
[cfg.B_N,cfg.j0]=scope_B_for_N(cfg.N);
cfg.kGrid=[3 5 7];
cfg.rGrid=0:3;
cfg.bGrid=[0.65 0.80 1.00 1.20 1.50];
cfg.generators={'Normal','Logistic','Laplace','Sech','Uniform','Cauchy'};
cfg.signals={'Doppler','HeaviSine','Bumps','Blocks'};
cfg.snrs=snrsUser;
cfg.snrDefinition='std_signal_over_std_noise_iidN01';
cfg.M=M_requested; cfg.batchSize=20;
cfg.seed=20261026+cfg.N; cfg.shift=2;
cfg.makePlots=true; cfg.resume=true;
if smoke
    cfg.generators={'Normal','Logistic'};
    cfg.signals={'Doppler','Blocks'};
    cfg.snrs=5;
    cfg.M=M_requested; cfg.batchSize=min(3,cfg.M);
    outDir=fullfile(outputBase,sprintf('Exp2_N%d_SMOKE',cfg.N));
    calibFile=fullfile(outDir,'calibration_E2.mat');
else
    outDir=fullfile(outputBase,sprintf('Exp2_N%d_M%d%s',cfg.N,cfg.M,snrTag));
    calibFile=fullfile(outDir,'calibration_E2.mat');
end
if ~exist(outDir,'dir'), mkdir(outDir); end
plotDir=fullfile(outDir,'plots');
if ~exist(plotDir,'dir'), mkdir(plotDir); end
 
cand=e2_make_candidates(cfg);
e2_validate_divergence(); % finite-difference check for all six CDF derivatives
names=e2_method_names(cfg.generators);
fCount=numel(cfg.generators); cCount=numel(cand.k);
assert(cCount==numel(cfg.kGrid)*numel(cfg.rGrid)*numel(cfg.bGrid));
fprintf('Experiment II (%s): M=%d, %d CDFs, %d candidates per CDF.\n',...
    mode,cfg.M,fCount,cCount);
 
% Exactly the same quadrature convention as Experiment I, except k varies.
% Calibration: 6 x 3 k x 4 radii x (J-j0) levels.
calibMeta=struct('generators',{cfg.generators},'kGrid',cfg.kGrid,...
    'rGrid',cfg.rGrid,'B_N',cfg.B_N,'j0',cfg.j0,'J',cfg.J);
if exist(calibFile,'file')
    z=load(calibFile,'a','calibMeta');
    assert(isfield(z,'a') && isfield(z,'calibMeta') && ...
       isequaln(z.calibMeta,calibMeta),...
       'The cached calibration belongs to different settings: %s',calibFile);
    a=z.a;
    fprintf('Loaded null-risk constants from %s\n',calibFile);
else
    fprintf('Calibrating Gaussian null-risk constants (one-time task)...\n');
    a=e2_compute_constants(cfg);
    save(calibFile,'a','calibMeta','-v7');
end
assert(isequal(size(a),[fCount,numel(cfg.kGrid),numel(cfg.rGrid),cfg.J-cfg.j0]));
e2_check_k5_constants(a,cfg);
e2_write_constants(cfg,a,fullfile(outDir,'Experiment2_ScaleConstants.csv'));
 
% Build the supplied periodic orthonormal wavelet matrices ONCE.
Wbank=cell(1,3);
Wbank{1}=e2_wavelet_matrix('Symmlet',4,cfg);
Wbank{2}=e2_wavelet_matrix('Daubechies',6,cfg);
Wbank{3}=e2_wavelet_matrix('Haar',0,cfg);
results=cell(numel(cfg.signals),numel(cfg.snrs));
totalTimer=tic;
for is=1:numel(cfg.signals)
    signalName=cfg.signals{is};
    [f0,filterName,whichW]=e2_make_DJ_signal(signalName,cfg.N);
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
                isequaln(result.studySettings,e2_study_settings(cfg)),...
                'Existing completed checkpoint is incompatible with the settings.');
            fprintf('RESUMED completed: %s SNR=%g\n',signalName,snr);
        else
            f=scope_scale_to_SNR(f0,snr);
            fprintf('  Verified signal SD = %.8f; population noise SD = 1.\n',std(f,1));
            fprintf('\nStarting %s | SNR=%g | M=%d\n',signalName,snr,cfg.M);
            rng(cfg.seed+100*is+iq,'twister');
            result=e2_simulate_condition(f,W,a,cfg,cand,names,...
                signalName,filterName,snr,partialFn);
            save(fn,'result','-v7');
            if exist(partialFn,'file'), delete(partialFn); end
        end
        results{is,iq}=result;
        e2_write_mse_csv(result,fullfile(outDir,sprintf(...
            'MSE_replications_%s_SNR%d.csv',signalName,snr)));
        e2_write_selection_csv(result,fullfile(outDir,sprintf(...
            'SelectedParams_%s_SNR%d.csv',signalName,snr)));
        if cfg.makePlots
            e2_condition_plots(result,cfg,plotDir);
        end
        fprintf('Done: %s SNR=%g; total elapsed %.1f sec\n',...
            signalName,snr,toc(totalTimer));
    end
end
 
[summary,paired,selFreq,tripleFreq,ranking]=e2_compile_tables(results,cfg,names);
summary.N = repmat(cfg.N,height(summary),1);
summary.RMSE = sqrt(summary.AMSE);
writetable(summary,fullfile(outDir,'Experiment2_AMSE_all.csv'));
writetable(paired,fullfile(outDir,'Experiment2_PairedDifferences.csv'));
writetable(selFreq,fullfile(outDir,'Experiment2_SelectionFrequencies.csv'));
writetable(tripleFreq,fullfile(outDir,'Experiment2_SelectedTriples.csv'));
writetable(ranking,fullfile(outDir,'Experiment2_AggregateRanking.csv'));
if cfg.makePlots
    e2_aggregate_plot(ranking,cfg,plotDir);
end
save(fullfile(outDir,'Experiment2_Results.mat'),...
    'cfg','a','cand','results','-v7.3');
e2_write_results_readme(outDir,cfg);
fprintf('\nExperiment II finished: %s\n',outDir);
fprintf('AMSE: Experiment2_AMSE_all.csv | Paired: Experiment2_PairedDifferences.csv\n');
fprintf('Total elapsed %.1f seconds.\n',toc(totalTimer));
end
 
% ========================================================================
function cand=e2_make_candidates(cfg)
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
function W=e2_wavelet_matrix(filterName,par,cfg)
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
function [f0,label,whichW]=e2_make_DJ_signal(name,N)
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
function a=e2_compute_constants(cfg)
% Gaussian-null per-coefficient risk; no truth/signal enters calibration.
% Experiment I uses integral(...,0,Inf) and fzero. Here repeated
% for k={3,5,7} to ensure identical k=5 constants.
nF=numel(cfg.generators); nK=numel(cfg.kGrid);
nR=numel(cfg.rGrid); nJ=cfg.J-cfg.j0;
a=zeros(nF,nK,nR,nJ);
for iF=1:nF
 for ik=1:nK
  for ir=1:nR
   for jj=1:nJ
    j=cfg.j0+jj-1;
    target=cfg.B_N/2^j;
    fun=@(x)e2_nullrisk(x,cfg.generators{iF},cfg.rGrid(ir),cfg.kGrid(ik));
    upper=1;
    while fun(upper)<target
       upper=upper*2;
       assert(upper<1e6,'Null-risk root failed to bracket.');
    end
    a(iF,ik,ir,jj)=fzero(@(x)fun(x)-target,[0 upper],...
          optimset('TolX',1e-10,'Display','off'));
   end
  end
 end
 fprintf('  Calibrated %-9s | k={3,5,7}, r=0:3\n',cfg.generators{iF});
end
end
 
function risk=e2_nullrisk(a,F,r,k)
if a<=0, risk=0; return; end
m=2*r+1; nu=m/2;
logC=-nu*log(2)-gammaln(nu)-log(m);
fun=@(u)exp(nu.*log(max(u,realmin))-.5*u+logC) ...
    .*e2_cdf_centered(F,a.*sqrt(u/m)).^(2*k);
risk=integral(fun,0,Inf,'AbsTol',2e-10,'RelTol',1e-8);
risk=max(0,min(1,risk));
end
 
function [g,gprime]=e2_cdf_centered(F,u)
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
function e2_validate_divergence()
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
  [g,gp]=e2_cdf_centered(Fs{f},u);
  for k=[3 5 7]
   analytic=g^k+(k*lambda*d(center)^2/(mm*S))*g^(k-1)*gp;
   plus=d; plus(center)=plus(center)+epsFD;
   minus=d; minus(center)=minus(center)-epsFD;
   sp=sqrt(sum(plus(ids).^2)/mm);
   sm=sqrt(sum(minus(ids).^2)/mm);
   fp=plus(center)*e2_cdf_centered(Fs{f},lambda*sp)^k;
   fm=minus(center)*e2_cdf_centered(Fs{f},lambda*sm)^k;
   numerical=(fp-fm)/(2*epsFD);
   assert(abs(numerical-analytic)<5e-6,...
      'SURE derivative validation failed: F=%s, k=%d, r=%d',Fs{f},k,rr);
  end
 end
end
fprintf('Passed finite-difference checks of N-SCOPE diagonal divergence.\n');
end
 
% ========================================================================
function e2_check_k5_constants(a,cfg)
% Validate the finest-level null-risk root for k=5 at each generator/radius.
ik=find(cfg.kGrid==5,1);
for iF=1:numel(cfg.generators)
 for ir=1:numel(cfg.rGrid)
  got=e2_nullrisk(a(iF,ik,ir,end),cfg.generators{iF},cfg.rGrid(ir),5);
  assert(abs(got-cfg.B_N/2^(cfg.J-1))<1e-5,'E2 null-risk calibration mismatch');
 end
end
end
function e2_write_constants(cfg,a,fileName)
F={}; K=[]; R=[]; J=[]; TARGET=[]; B=[]; A=[];
for f=1:numel(cfg.generators)
 for ik=1:numel(cfg.kGrid)
  for ir=1:numel(cfg.rGrid)
   for jj=1:(cfg.J-cfg.j0)
    j=cfg.j0+jj-1;
    F{end+1,1}=cfg.generators{f}; %#ok<AGROW>
    K(end+1,1)=cfg.kGrid(ik); %#ok<AGROW>
    R(end+1,1)=cfg.rGrid(ir); %#ok<AGROW>
    J(end+1,1)=j; %#ok<AGROW>
    TARGET(end+1,1)=cfg.B_N/2^j; %#ok<AGROW>
    B(end+1,1)=cfg.B_N; %#ok<AGROW>
    A(end+1,1)=a(f,ik,ir,jj); %#ok<AGROW>
   end
  end
 end
end
T=table(F,K,R,J,B,TARGET,A,'VariableNames',...
   {'Generator','K','Radius','Level_j','Budget_BN','TargetNullRisk','a_NullRisk'});
writetable(T,fileName);
end
 
% ========================================================================
function S=e2_study_settings(cfg)
% Protect checkpoints from accidental reuse after tuning-grid changes.
S=struct('N',cfg.N,'M',cfg.M,'J',cfg.J,'j0',cfg.j0,...
  'B_N',cfg.B_N,'kGrid',cfg.kGrid,'rGrid',cfg.rGrid,...
  'bGrid',cfg.bGrid,'generators',{cfg.generators},...
  'batchSize',cfg.batchSize,'seed',cfg.seed,'shift',cfg.shift,...
  'snrDefinition',cfg.snrDefinition);
end
 
function result=e2_simulate_condition(f,W,a,cfg,cand,names,tag,filterName,snr,partialFn)
N=cfg.N; M=cfg.M; nbMethods=numel(names);
nF=numel(cfg.generators); nbCandidates=numel(cand.k);
 
% State contains all in-progress values, including the exact RNG state.
% This makes restarting mid-condition reproducible and deterministic.
if cfg.resume && exist(partialFn,'file')
    z=load(partialFn,'state'); state=z.state;
    assert(state.M==M && state.nF==nF && state.nbCandidates==nbCandidates && ...
        isfield(state,'studySettings') && ...
        isequaln(state.studySettings,e2_study_settings(cfg)),...
        'Incompatible partial checkpoint. Move/delete %s',partialFn);
    rng(state.rngState);
    fprintf('  Resuming PARTIAL at replicate %d/%d.\n',state.completed,M);
else
    state=struct(); state.M=M; state.nF=nF;
    state.nbCandidates=nbCandidates; state.completed=0;
    state.studySettings=e2_study_settings(cfg);
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
    db=cell(cfg.J-cfg.j0,1); tr=cell(cfg.J-cfg.j0,1); S=cell(cfg.J-cfg.j0,4);
    for jj=1:(cfg.J-cfg.j0)
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
        % Evaluate EVERY candidate, summing SURE over all J-j0 levels.
        % Scale is held fixed with respect to d_jl in the derivative.
        for c=1:nbCandidates
            ir=cand.ir(c); ik=cand.ik(c);
            r=cand.r(c); k=cand.k(c); b=cand.b(c);
            mm=2*r+1;
            for jj=1:(cfg.J-cfg.j0)
                dj=db{jj}; Sj=S{jj,ir};
                lam=b*a(iF,ik,ir,jj)./sigmahat; % 1 x nb
                u=bsxfun(@times,Sj,lam);
                [g,gp]=e2_cdf_centered(F,u);
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
            loss=e2_evaluate_selected_MSE(db,tr,S,errScaling,...
                  sigmahat,F,a,iF,cand,choices{meth},N);
            mseBatch(:,(iF-1)*4+meth)=loss(:);
        end
    end
    % Common universal-threshold references; same data, same sigmahat.
    hardLoss=errScaling; softLoss=errScaling;
    threshold=sigmahat*sqrt(2*log(N));
    for jj=1:(cfg.J-cfg.j0)
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
   'waveletShift',cfg.shift,'studySettings',e2_study_settings(cfg));
end
 
function loss=e2_evaluate_selected_MSE(db,tr,S,errScaling,sigmahat,F,a,iF,cand,idx,N)
% Squared error computed only for the *selected* candidate. Crucially, the
% TRUE coefficient tr is never accessed during SURE candidate selection.
nb=numel(idx); loss=errScaling;
% The coefficient blocks already determine the number of detail levels.
% Do not access cfg here: this helper receives db/tr/S, not cfg.
nLevels=numel(db);
assert(numel(tr)==nLevels && size(S,1)==nLevels, ...
    'Inconsistent wavelet-level arrays in Experiment II.');
assert(numel(loss)==nb && all(idx>=1 & idx<=numel(cand.k)), ...
    'Invalid candidate selection indices in Experiment II.');
for jj=1:nLevels
    dj=db{jj}; target=tr{jj}; theta=zeros(size(dj));
    % Group replications with the same candidate; vectorize within group.
    distinct=unique(idx);
    for ii=1:numel(distinct)
        c=distinct(ii); cols=find(idx==c);
        ir=cand.ir(c); ik=cand.ik(c);
        lam=cand.b(c)*a(iF,ik,ir,jj)./sigmahat(cols);
        u=bsxfun(@times,S{jj,ir}(:,cols),lam);
        theta(:,cols)=dj(:,cols).*e2_cdf_centered(F,u).^cand.k(c);
    end
    loss=loss+sum(bsxfun(@minus,theta,target).^2,1);
end
loss=loss/N;
end
 
% ========================================================================
function names=e2_method_names(Fs)
nF=numel(Fs); names=cell(nF*4+2,1);
suffix={'FixedSCOPE','FixedNSCOPE','SCOPE_SURE','NSCOPE_SURE'};
for f=1:nF
 for meth=1:4
    names{(f-1)*4+meth}=sprintf('%s_%s',Fs{f},suffix{meth});
 end
end
names{end-1}='UniversalHard'; names{end}='UniversalSoft';
end
 
function e2_write_mse_csv(R,fileName)
T=table((1:R.M)',R.sigmahat,...
     'VariableNames',{'Replication','SigmaHat'});
for c=1:numel(R.names)
    T.(R.names{c})=R.mse(:,c);
end
writetable(T,fileName);
end
 
function e2_write_selection_csv(R,fileName)
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
function [overall,paired,selFreq,tripleFreq,ranking]=e2_compile_tables(results,cfg,names)
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
function e2_condition_plots(R,cfg,plotDir)
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
e2_export_figure(fig,base);
 
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
e2_export_figure(fig,base);
 
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
e2_export_figure(fig,base);
end
 
function e2_aggregate_plot(ranking,cfg,plotDir)
nF=numel(cfg.generators); vals=zeros(nF,4);
for fi=1:nF
 for meth=1:4
    nm=e2_method_names(cfg.generators(fi));
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
e2_export_figure(fig,fullfile(plotDir,'E2_AggregateRelativeAMSE'));
end
 
function e2_export_figure(fig,basename)
print(fig,[basename '.png'],'-dpng','-r180');
print(fig,[basename '.pdf'],'-dpdf','-painters');
close(fig);
end
 
function e2_write_results_readme(outDir,cfg)
fid=fopen(fullfile(outDir,'Experiment2_README.txt'),'w');
assert(fid>0);
fprintf(fid,'SCOPE Experiment II, October 2026 | n=%d, M=%d, sigma=1\n',cfg.N,cfg.M);
fprintf(fid,'DJ signals %s | SNR=std(f,1)/std(epsilon) %s\n',strjoin(cfg.signals,', '),mat2str(cfg.snrs));
fprintf(fid,'std(epsilon)=1; std(f,1)=target SNR; Gaussian iid N(0,1), divisor n.\n');
fprintf(fid,'CDFs %s\n',strjoin(cfg.generators,', '));
fprintf(fid,'kGrid=%s; rGrid=%s; bGrid=%s\n',...
    mat2str(cfg.kGrid),mat2str(cfg.rGrid),mat2str(cfg.bGrid));
fprintf(fid,'Global tuning across %d detail levels j=%d,...,%d.\n',...
    cfg.J-cfg.j0,cfg.j0,cfg.J-1);
fprintf(fid,'lambda_j=b*R^{-1}_{F,r,k}(B_N/2^j)/sigmahat, B_N=%g.\n',cfg.B_N);
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
 
 
function e3_run_SCOPE_Experiment3(mode,N_requested,M_requested,outputBase,snrsUser,snrTag)
% RUN_SCOPE_EXPERIMENT3   Independent Monte Carlo, SCOPE against 20 benchmarks.
%
% RUN through the canonical SCOPE_Exp123_BN entry point only.
%
% WAVELET: N=512,1024,2048; orthonormal PERIODIC DWT using the user's
% Wavmat and MakeONFilterExt; detail levels j=3..J-1 (eight scaling terms
% untouched). Symmlet-4 for Doppler/HeaviSine, Daub-6 for Bumps,
% Haar for Blocks. Gaussian noise iid N(0,1); SNR=std(f,1)/std(epsilon)=3,5,7.
% Standard deviation of epsilon is 1 (population); std(f,1)=SNR.
% The name SNR is kept although this is root SNR (an SD ratio, not a power ratio).
%
% PRESELECTED four procedures (never reselect using these replications):
%  1. Fixed SCOPE: Normal F, k=5,r=0,b=1.
%  2. Fixed N-SCOPE: Normal F, k=5,r=1,b=1.
%  3. SCOPE SURE: Normal F; globally tune k={3,5,7}, r=0, b-grid.
%  4. N-SCOPE SURE: Normal F; globally tune k, r={0,1,2,3}, b-grid.
% These four recipes were selected using the FRESH RSNR Experiment I/II
% rankings; do NOT use independent Experiment III data to reselect them.
% Same B_N/2^j null-risk calibration as Experiment II.
%
% OTHER METHODS: The 20 listed in Section 6, Table competitors.
% Legacy 1D comparator routines were carried over from the previously
% supplied 17-/23-method MATLAB implementations. The additional
% procedures are implemented in E3_extended_methods below. IMPORTANT: several
% image-origin or complex methods are PROVISIONAL one-dimensional variants.
% Their results are to be treated as an IMPLEMENTATION PILOT, NOT as a
% vetted exact reproduction of the named published algorithm. Review
% METHOD_IMPLEMENTATION_NOTES.csv before citing outcomes in the manuscript.
%
% Independent test set (new seed distinct from Experiments I and II),
% paired noise realizations across all 26 estimators, no oracle parameter
% fitting. MSE assessed in the wavelet domain (Parseval identity).
%
% MATLAB R2018b+; no Wavelet Toolbox or Statistics Toolbox. Do not rely
% on prior output directories. Checkpoint files allow interrupted resumption.
if nargin<1 || isempty(mode),mode='full';end
if isstring(mode),mode=char(mode);end
smoke=strcmpi(mode,'smoke');assert(smoke || strcmpi(mode,'full'));
root=fileparts(mfilename('fullpath'));addpath(root,'-begin');
% Four Normal-generator SCOPE procedures were PRESELECTED from RSNR
% Experiments I and II (user-reported ranking, 2026-10-28), before E3.
% DO NOT require E1/E2 output CSVs in this directory: these were produced
% in the user's separate Experiment I/II folders. Do not re-select based
% on the independent E3 outcomes.
fprintf('Frozen E3: Normal fixed r=0, Normal fixed r=1, Normal SCOPE SURE, Normal N-SCOPE SURE.\n');
cfg=struct();cfg.N=N_requested;cfg.J=round(log2(cfg.N));cfg.M=M_requested;
[cfg.B_N,cfg.j0]=scope_B_for_N(cfg.N);
cfg.signals={'Doppler','HeaviSine','Bumps','Blocks'};
cfg.snrs=snrsUser;cfg.snrDefinition='std_signal_over_std_noise_iidN01';
cfg.kGrid=[3 5 7];cfg.rGrid=0:3;
cfg.bGrid=[.65 .80 1 1.20 1.50];
cfg.Fs={'Normal'};cfg.seed=20271126+cfg.N;cfg.shift=2;
assert(isequal(cfg.Fs,{'Normal'}),'Frozen generator specification has changed.');
cfg.resume=true;cfg.plot=true;cfg.saveEvery=50;
% Checkpoint signatures include the calibration budget. No migration of old
% E3 checkpoints is allowed: ALL FOUR SCOPE methods need recalculation.
cfg.codeRevision=sprintf('E3_BNrounded_Normal_N%d_M%d_B%.6g_j0%d_FINAL26_v1', ...
 cfg.N,cfg.M,cfg.B_N,cfg.j0);
if smoke
 cfg.M=2;cfg.signals={'Doppler','Blocks'};cfg.snrs=5;
 outDir=fullfile(outputBase,sprintf('Exp3_N%d_SMOKE',cfg.N));
else
 outDir=fullfile(outputBase,sprintf('Exp3_N%d_M%d%s',cfg.N,cfg.M,snrTag));
end
if ~exist(outDir,'dir'),mkdir(outDir);end
plotDir=fullfile(outDir,'plots');if ~exist(plotDir,'dir'),mkdir(plotDir);end
names=namesE3(); assert(numel(names)==26);
write_audit(fullfile(outDir,'METHOD_IMPLEMENTATION_NOTES.csv'),names);
assert(all(cellfun(@isvarname,names)), 'Method names must be table-variable safe');
% Calibration uses observed-data-independent Gaussian chi-square risk.
% IMPORTANT: the old E3_calibration.mat contains LOGISTIC/SECH constants.
% Do NOT reuse it: the new winner is Normal for all four SCOPE recipes.
% The reference cache here uses the exact Normal slice of Experiment II's
% independent Gaussian-null calibration (identical eta, k, r and j).
calibFile=fullfile(outDir,'calibration_E3_Normal.mat');
% Cache signature includes the generator, preventing accidental reuse of
% logistic/sech coefficients even when eta, k, r and levels coincide.
calibSignature=[cfg.kGrid cfg.rGrid cfg.B_N cfg.N cfg.J cfg.j0];
generatorName=cfg.Fs{1};
if exist(calibFile,'file')
 c=load(calibFile,'a','calibSignature','generatorName');
 assert(isfield(c,'a') && isfield(c,'calibSignature') && ...
    isfield(c,'generatorName') && ...
    isequal(c.calibSignature,calibSignature) && ...
    strcmp(c.generatorName,generatorName), ...
    'Wrong E3 Gaussian-null calibration cache; replace %s',calibFile);
 a=c.a;
else
 fprintf('Computing Normal Gaussian null-risk scales (once)...\n');
 a=calibrate_E3(cfg);
 save(calibFile,'a','calibSignature','generatorName','-v7');
end
assert(isequal(size(a),[1,numel(cfg.kGrid),numel(cfg.rGrid),cfg.J-cfg.j0]), ...
 'Normal E3 calibration array has incorrect dimensions.');
validate_a(a,cfg);
Wbank=cell(3,1); Hbank=cell(3,1);
filters={'Symmlet','Daubechies','Haar'}; par=[4 6 0];
for i=1:3
 Hbank{i}=MakeONFilterExt(filters{i},par(i));
 Wbank{i}=Wavmat(Hbank{i},cfg.N,cfg.J-cfg.j0,cfg.shift);
 v=sin((1:cfg.N)'*1.35);
 assert(norm(Wbank{i}'*(Wbank{i}*v)-v)/norm(v)<1e-8, ...
   'Wavelet matrix not orthonormal for %s',filters{i});
end
results=cell(numel(cfg.signals),numel(cfg.snrs)); start=tic;
for si=1:numel(cfg.signals)
 [f0,whichW]=signal_E3(cfg.signals{si},cfg.N);W=Wbank{whichW};h=Hbank{whichW};
 for iq=1:numel(cfg.snrs)
  snr=cfg.snrs(iq);f=scope_scale_to_SNR(f0,snr);
  fprintf('  Verified signal SD = %.8f; population noise SD = 1.\n',std(f,1));
  cp=fullfile(outDir,sprintf('checkpoint_%s_SNR%d.mat',cfg.signals{si},snr));
  fprintf('\n%s | SNR %g | M=%d | %d estimators\n',cfg.signals{si},snr,cfg.M,numel(names));
  if cfg.resume && exist(cp,'file')
   q=load(cp,'st'); assert(isfield(q,'st'),'Missing E3 state: %s',cp);
   st=q.st;
   assert(isfield(st,'codeRevision') && ...
      st.M==cfg.M && st.SNR==snr && ...
      strcmp(st.Signal,cfg.signals{si}) && ...
      isfield(st,'snrDefinition') && ...
      strcmp(st.snrDefinition,cfg.snrDefinition), ...
      'E3 checkpoint design does not match: %s',cp);
   assert(strcmp(st.codeRevision,cfg.codeRevision), ...
      'Old E3 checkpoint is incompatible with Section 4 rounded B_N: %s',cp);
   assert(isequal(st.names,names) && ...
      size(st.mse,2)==26 && strcmp(st.codeRevision,cfg.codeRevision), ...
      'E3 checkpoint method/revision mismatch: %s',cp);
   rng(st.rngState);
   fprintf('Resuming at replicate %d.\n',st.completed+1);
  else
   rng(cfg.seed+100*si+iq,'twister');
   st=struct('mse',NaN(cfg.M,numel(names)),'sigmaHat',NaN(cfg.M,1),...
       'kScope',NaN(cfg.M,1),'bScope',NaN(cfg.M,1),...
       'kNScope',NaN(cfg.M,1),'rNScope',NaN(cfg.M,1),...
       'bNScope',NaN(cfg.M,1),'names',{names},...
       'M',cfg.M,'SNR',snr,'Signal',cfg.signals{si},'completed',0,...
       'rngState',rng,'rep1',struct(),'codeRevision',cfg.codeRevision,...
       'snrDefinition',cfg.snrDefinition);
  end
  truth=W*f;
  for it=(st.completed+1):cfg.M
    y=f+randn(cfg.N,1);d=W*y;
    sigma=median(abs(d(cfg.N/2+1:end)))/.6744897501960817;
    assert(sigma>0 && isfinite(sigma));st.sigmaHat(it)=sigma;
    if it==1,st.rep1=struct('f',f,'y',y,'sigma',sigma);end
    S=neighborhoods_E3(d,cfg); % normalized local RMS for each (j,r)
    % Four pre-selected SCOPE procedures. Source data for SURE are d,S;
    % the clean signal 'truth' enters ONLY after the estimated rule chosen.
    sc=zeros(cfg.N,4);
    sc(:,1)=apply_E3(d,S,sigma,a,cfg,1,5,0,1);
    sc(:,2)=apply_E3(d,S,sigma,a,cfg,1,5,1,1);
    [k0,~,b0]=select_E3(d,S,sigma,a,cfg,1,0);
    sc(:,3)=apply_E3(d,S,sigma,a,cfg,1,k0,0,b0); % Normal is generator index 1
    [k1,r1,b1]=select_E3(d,S,sigma,a,cfg,1,cfg.rGrid);
    sc(:,4)=apply_E3(d,S,sigma,a,cfg,1,k1,r1,b1);
    st.kScope(it)=k0;st.bScope(it)=b0;
    st.kNScope(it)=k1;st.rNScope(it)=r1;st.bNScope(it)=b1;
    % Eight existing MATLAB comparator algorithms: table IDs 1,2,3,4,
    % 5,6,9,13; one-dimensional adaptations must be reported as such.
    leg=E3_legacy_methods(d,sigma,y,h,W);
    % Twelve further methods: table IDs 7,8,10,11,12,14,15,16,17,
    % 18,19,20. See audit CSV for provisional status.
    ext=E3_extended_methods(d,sigma);
    dh=zeros(cfg.N,numel(names));dh(:,1:4)=sc;
    legids=[1 2 3 4 5 6 9 13]+4;
    extids=[7 8 10 11 12 14 15 16 17 18 19 20]+4;
    dh(:,legids)=leg;dh(:,extids)=ext;
    % Universal hard and soft, same sigmahat, same orthogonal DWT.
    H=d;T=sigma*sqrt(2*log(cfg.N));
    H(2^cfg.j0+1:end)=d(2^cfg.j0+1:end).*(abs(d(2^cfg.j0+1:end))>T);
    Q=d;Q(2^cfg.j0+1:end)=sign(d(2^cfg.j0+1:end)).*max(abs(d(2^cfg.j0+1:end))-T,0);
    dh(:,25)=H;dh(:,26)=Q;
    assert(all(isfinite(dh(:))),'Non-finite estimator');
    st.mse(it,:)=sum(bsxfun(@minus,dh,truth).^2,1)/cfg.N;
    st.completed=it;
    if mod(it,cfg.saveEvery)==0 || it==cfg.M || it==1
      st.rngState=rng;save(cp,'st','-v7');
      fprintf('  %s SNR %d: %4d/%4d (%.1f%%); elapsed %.1f sec\n',...
        cfg.signals{si},snr,it,cfg.M,100*it/cfg.M,toc(start));
    end
  end
  results{si,iq}=st;
  write_rep_csv(st,outDir);
  if cfg.plot,plot_condition(st,plotDir);end
 end
end
write_summaries(results,cfg,names,outDir,plotDir);
save(fullfile(outDir,'Experiment3_Results.mat'),'cfg','results','a','names','-v7.3');
fprintf('\nExperiment III finished: %s; elapsed %.1f sec.\n',outDir,toc(start));
end
 
function names=namesE3()
% CANONICAL 26 columns, matching the N=1024 x 3-SNR paper table.
names={'SCOPE_Fixed_Normal','NSCOPE_Fixed_Normal_r1',...
'SCOPE_SURE_Normal','NSCOPE_SURE_Normal',...
'BAMS_Author_C2_Gamma2','EBayesMedian','RaisedCosineBayes','EpanechnikovBayes',...
'SmoothSCAD_stabilized_author','SureShrink_Hybrid','FDR_BH_q05','EBBlock_1D_proxy',...
'BayesShrink_1D','NeighCoeff','NeighShrink_SURE_1D',...
'SURELET_1D_orthonormal','BiShrink_1D','NeighBlock','SureBlock_hybrid_CZ09',...
'BlockJS','Jeffreys_ABE_FigueiredoNowak','LocalWiener_1D',...
'Context_CYV2000_1D','HMT_1D_EM_proxy','UniversalHard','UniversalSoft'};
end

function st=e3_upgrade_checkpoint_FINAL(st,f,W,cfg,si,iq,names)
% Repair ONLY #5 BAMS and #23 CYV, safely reusing the other 24 methods.
% Legacy 28-column checkpoints had BAMSO(#27), CYV(#28), an old proxy(#23).
% 26-column checkpoints had original BAMS(#5), old proxy(#23).
% We recompute BOTH columns from exact seeded observations in all cases.
N=cfg.N; J=cfg.J; j0=cfg.j0;
assert(size(st.mse,2)==26 || size(st.mse,2)==28, ...
   'Checkpoint must have 26 or 28 method columns.');
assert(isfield(st,'names') && numel(st.names)==size(st.mse,2) && ...
   isfield(st,'completed') && st.completed>=0 && ...
   st.completed<=cfg.M && st.completed==floor(st.completed) && ...
   isfield(st,'sigmaHat') && numel(st.sigmaHat)==cfg.M && ...
   isfield(st,'rngState'), 'Incomplete E3 checkpoint metadata.');
% All names but the deliberately replaced BAMS and CYV must coincide.
idx=[1:4 6:22 24:26];
oldNames=st.names;
assert(isequal(oldNames(idx),legacy_e3_names_excluding_replacements()) || ...
       isequal(oldNames(idx),names(idx)), ...
       'Legacy E3 checkpoint method ordering unexpected.');
newBAMS=NaN(cfg.M,1);newCYV=NaN(cfg.M,1);
oldRNG=rng; restoreRNG=onCleanup(@()rng(oldRNG)); %#ok<NASGU>
rng(cfg.seed+100*si+iq,'twister');
truth=W*f;
for it=1:st.completed
   y=f+randn(N,1); d=W*y;
   sigma=median(abs(d(N/2+1:N)))/0.6744897501960817;
   assert(abs(sigma-st.sigmaHat(it))<1e-8*max(1,abs(sigma)), ...
       'Noise replay discrepancy at replication %d; NO checkpoint update.',it);
   if it==1
      assert(isfield(st,'rep1') && isfield(st.rep1,'y') && ...
         max(abs(y-st.rep1.y))<1e-8, ...
         'First observation does not replay; NO checkpoint update.');
   end
   b=scope_BAMS_author_coefficients(d,J,J-j0,2,2);
   c=cyv_context_coefficients(d,sigma,J,j0);
   newBAMS(it)=sum((b-truth).^2)/N;
   newCYV(it)=sum((c-truth).^2)/N;
   assert(isfinite(newBAMS(it)) && isfinite(newCYV(it)), ...
      'Non-finite replay estimate at replication %d.',it);
end
% Do not mutate st until ALL replay checks have passed.
newMSE=st.mse(:,1:26);
newMSE(1:st.completed,5)=newBAMS(1:st.completed);
newMSE(1:st.completed,23)=newCYV(1:st.completed);
st.mse=newMSE;st.names=names;
if st.completed>0
 fprintf('  Verified %d trials; BAMS AMSE %.6f; CYV AMSE %.6f\n', ...
   st.completed,mean(newBAMS(1:st.completed)),mean(newCYV(1:st.completed)));
end
end

function names=legacy_e3_names_excluding_replacements()
% Original 26-column names excluding positions 5 and 23.
a={'SCOPE_Fixed_Normal','NSCOPE_Fixed_Normal_r1',...
'SCOPE_SURE_Normal','NSCOPE_SURE_Normal',...
'EBayesMedian','RaisedCosineBayes','EpanechnikovBayes',...
'SmoothSCAD_stabilized_author','SureShrink_Hybrid','FDR_BH_q05',...
'EBBlock_1D_proxy','BayesShrink_1D','NeighCoeff','NeighShrink_SURE_1D',...
'SURELET_1D_orthonormal','BiShrink_1D','NeighBlock',...
'SureBlock_hybrid_CZ09','BlockJS','Jeffreys_ABE_FigueiredoNowak',...
'LocalWiener_1D','HMT_1D_EM_proxy','UniversalHard','UniversalSoft'};
names=a;
end

function z=cyv_context_coefficients(d,sigma,J,j0)
% Chang, Yu & Vetterli (2000), IEEE Trans. Image Process. 9(9):1522-1531,
% "Spatially adaptive wavelet thresholding with context modeling",
% specialized to 1D orthonormal periodic DWT, levels j0..J-1.
% For each detail level (subband):
%  1. context Z_k = w'|u_k|, u_k = |neighbours| EXCLUDING d_k itself:
%     d_{k-1}, d_{k+1}, d_{k-2}, d_{k+2} (periodic) and the parent
%     d_{j-1,ceil(k/2)} (no parent at the coarsest level j0);
%     w = least-squares fit of |d_k| on u_k (no intercept);
%  2. sort coefficients by Z; for each, take the L coefficients nearest in
%     context order, L = max(50, 2% of the level size) (all if smaller);
%     sigmaY^2 = mean(d^2) over them; sigmaX = sqrt(max(sigmaY^2-sigma^2,0));
%  3. soft threshold with T_k = sigma^2/sigmaX_k (T=Inf, i.e. 0, if sigmaX=0).
% The 2D method uses 8 spatial neighbours + parent; the 1D neighbour set
% above is our analogue. Levels too small to fit the context regression
% (n < 2p+2) use BayesShrink, which is CYV with L = level size.
d=d(:); n=numel(d); assert(2^J==n && j0>=1 && j0<J);
z=d;
for j=j0:J-1
   ind=(2^j+1):2^(j+1);
   if j>j0, parent=d((2^(j-1)+1):2^j); else, parent=[]; end
   z(ind)=cyv_context_level(d(ind),parent,sigma);
end
end
 
function out=cyv_context_level(x,parent,sigma)
x=x(:); n=numel(x); ax=abs(x);
U=[circshift(ax,1),circshift(ax,-1),circshift(ax,2),circshift(ax,-2)];
if ~isempty(parent), U=[U,abs(repelem(parent(:),2))]; end
p=size(U,2);
if n<2*p+2
   out=bayesshrink_level(x,sigma); out=out(:); return;
end
G=U'*U;
if rcond(G)<1e-12, w=pinv(U)*ax; else, w=U\ax; end
Zc=U*w;
[~,order]=sort(Zc);
xs=x(order);
L=min(n,max(50,ceil(0.02*n)));
cs=[0;cumsum(xs.^2)];
i=(1:n)';
lo=min(max(i-floor(L/2),1),n-L+1); hi=lo+L-1;
sigmaY2=(cs(hi+1)-cs(lo))/L;
sigmaX=sqrt(max(sigmaY2-sigma^2,0));
T=Inf(n,1); pos=sigmaX>0; T(pos)=sigma^2./sigmaX(pos);
est=sign(xs).*max(abs(xs)-T,0);
out=zeros(n,1); out(order)=est;
end
 
function write_audit(p,names)
% Algorithm-level transparency: some 1D adaptations are not independently
% validated reproductions of the corresponding original research code.
status=repmat({'Published-method-inspired 1D adaptation; audit needed'},numel(names),1);
for i=1:4
 status{i}='SCOPE specification: fixed null risk / plug-in SURE';
end
status{5}='User-supplied BAMS.m posterior mean: Cscale=2,gamma=2,depth=J-3,level tau_j=tau_u sqrt(j)';
status{6}='Empirical Bayes median: implementation audit pending';
status{9}='Stabilized SmoothSCAD AUTHOR variant, not unrestricted published SURE';
status{10}='Hybrid SureShrink with sparse-level fallback';
for i=[12,24]
 status{i}='PROVISIONAL 1D proxy, not a verified original published algorithm';
end
status{16}='1D orthonormal SURE-LET adaptation, not original image/interscale algorithm';
status{20}='Cai BlockJS formula, 1D periodic boundary convention';
status{21}='Jeffreys/ABE 1D shrinkage formula';
status{23}=['Chang-Yu-Vetterli (2000) CONTEXT, 1D orthonormal: 4 neighbors+parent, ' ...
  'LS context, L=max(50,2%% level), adaptive soft threshold'];
status{25}='Universal hard threshold';status{26}='Universal soft threshold';
id=(1:numel(names))';
T=table(id,names(:),status,'VariableNames',{'Index','Procedure','ValidationStatus'});
writetable(T,p);
end

function [f,whichW]=signal_E3(name,N)
t=(1:N)'/N;
pos=[.10 .13 .15 .23 .25 .40 .44 .65 .76 .78 .81];
htB=[4 5 3 4 5 4.2 2.1 4.3 3.1 5.1 4.2];
htK=[4 -5 3 -4 5 -4.2 2.1 4.3 -3.1 2.1 -4.2];
w=[.005 .005 .006 .010 .010 .030 .010 .010 .005 .008 .005];
switch lower(name)
 case 'doppler', f=sqrt(t.*(1-t)).*sin(2*pi*1.05./(t+.05));whichW=1;
 case 'heavisine',f=4*sin(4*pi*t)-sign(t-.3)-sign(.72-t);whichW=1;
 case 'bumps'
  f=zeros(N,1); for q=1:numel(pos),f=f+htB(q)./(1+abs((t-pos(q))/w(q))).^4;end
  whichW=2;
 case 'blocks'
  f=zeros(N,1); for q=1:numel(pos),f=f+htK(q)*(t>=pos(q));end
  whichW=3;
 otherwise,error('Bad signal');
end
end
 
function a=calibrate_E3(cfg)
% B_N/2^j is the prescribed retained null risk per coefficient at level j.
a=zeros(numel(cfg.Fs),numel(cfg.kGrid),numel(cfg.rGrid),cfg.J-cfg.j0);
for f=1:numel(cfg.Fs)
 for ki=1:numel(cfg.kGrid)
  for ri=1:numel(cfg.rGrid)
   for jj=1:(cfg.J-cfg.j0)
    j=jj+cfg.j0-1; target=cfg.B_N/2^j;
    fn=@(v)nullrisk_E3(v,cfg.Fs{f},cfg.kGrid(ki),cfg.rGrid(ri));
    ub=1;while fn(ub)<target,ub=2*ub;assert(ub<1e6);end
    a(f,ki,ri,jj)=fzero(@(v)fn(v)-target,[0 ub],...
        optimset('TolX',1e-10,'Display','off'));
   end
  end
 end
 fprintf('Calibrated %s\n',cfg.Fs{f});
end
end
 
function R=nullrisk_E3(a,F,k,r)
if a<=0,R=0;return;end
m=2*r+1;nu=m/2;logC=-nu*log(2)-gammaln(nu)-log(m);
fn=@(u)exp(nu*log(max(u,realmin))-u/2+logC).*...
         (fstar_E3(F,a*sqrt(u/m))).^(2*k);
R=integral(fn,0,Inf,'AbsTol',2e-10,'RelTol',1e-8);
end
 
function g=fstar_E3(F,u)
u=max(u,0);
switch lower(F)
 case 'normal',g=erf(u/sqrt(2));
 case 'logistic',g=tanh(u);
 case 'sech',g=(2/pi)*atan(sinh(pi*u/2));
 otherwise,error('Unsupported generator %s',F);
end
g=min(max(g,0),1);
end
 
function [g,gp]=fstar_and_prime_E3(F,u)
g=fstar_E3(F,u);
switch lower(F)
 case 'normal',gp=sqrt(2/pi)*exp(-u.^2/2);
 case 'logistic',gp=1-g.^2;
 case 'sech',gp=1./cosh(pi*u/2);
end
end
 
function validate_a(a,cfg)
assert(all(isfinite(a(:)) & a(:)>0));
% Numerical null-risk verification on selected finest-level constants.
for f=1:numel(cfg.Fs)
 for ir=1:numel(cfg.rGrid)
  idx=find(cfg.kGrid==5);
  x=a(f,idx,ir,end);r=nullrisk_E3(x,cfg.Fs{f},5,cfg.rGrid(ir));
  assert(abs(r-cfg.B_N/2^(cfg.J-1))<1e-5,'Null-risk root verification failed');
 end
end
end
 
function S=neighborhoods_E3(d,cfg)
S=cell(cfg.J-cfg.j0,4);
for jj=1:(cfg.J-cfg.j0)
 j=cfg.j0+jj-1;ind=2^j+1:2^(j+1);x=d(ind);xx=x.^2;
 for ir=1:4
  r=ir-1;local=xx;
  for h=1:r,local=local+circshift(xx,h)+circshift(xx,-h);end
  S{jj,ir}=sqrt(local/(2*r+1));
 end
end
end
 
function out=apply_E3(d,S,sigma,a,cfg,iF,k,r,b)
ki=find(cfg.kGrid==k); ri=find(cfg.rGrid==r);assert(isscalar(ki)&&isscalar(ri));
out=d;
for jj=1:(cfg.J-cfg.j0)
 j=cfg.j0+jj-1;ind=(2^j+1):2^(j+1);
 u=(b*a(iF,ki,ri,jj)/sigma)*S{jj,ri};
 out(ind)=d(ind).*fstar_E3(cfg.Fs{iF},u).^k;
end
end
 
function [kbest,rbest,bbest]=select_E3(d,S,sigma,a,cfg,iF,rGrid)
best=Inf;kbest=NaN;rbest=NaN;bbest=NaN;sigma2=sigma^2;
for r=rGrid
 ri=find(cfg.rGrid==r);m=2*r+1;
 for k=cfg.kGrid
  ki=find(cfg.kGrid==k);
  for b=cfg.bGrid
   risk=0;
   for jj=1:(cfg.J-cfg.j0)
    j=cfg.j0+jj-1;ind=(2^j+1):2^(j+1);
    x=d(ind);s=S{jj,ri};lam=b*a(iF,ki,ri,jj)/sigma;
    u=lam*s;[g,gp]=fstar_and_prime_E3(cfg.Fs{iF},u);
    w=g.^k;
    correction=(k/m)*x.^2./max(s,realmin).*g.^(k-1).*gp*lam;
    der=w+correction;
    risk=risk+sum((x.*(w-1)).^2+2*sigma2*der-sigma2);
   end
   if risk<best,best=risk;kbest=k;rbest=r;bbest=b;end
  end
 end
end
end
 
function write_rep_csv(st,outDir)
T=table((1:st.M)',st.sigmaHat,st.kScope,st.bScope,...
 st.kNScope,st.rNScope,st.bNScope,...
 'VariableNames',{'Replicate','SigmaHat','ScopeK','ScopeB','NScopeK','NScopeR','NScopeB'});
for i=1:numel(st.names),T.(st.names{i})=st.mse(:,i);end
writetable(T,fullfile(outDir,sprintf('MSE_%s_SNR%d.csv',st.Signal,st.SNR)));
end
 
function write_summaries(results,cfg,names,outDir,plotDir)
rows={};paired={};scores=zeros(numel(results),numel(names));caseNo=0;
for si=1:size(results,1)
 for iq=1:size(results,2)
  st=results{si,iq};mu=mean(st.mse,1);sd=std(st.mse,0,1);
  caseNo=caseNo+1;scores(caseNo,:)=mu/mu(25);
  for m=1:numel(names)
   rows(end+1,:)={st.Signal,st.SNR,names{m},st.M,mu(m),sd(m),sd(m)/sqrt(st.M)}; %#ok<AGROW>
  end
  for m=5:24
   for a=1:4
    dif=st.mse(:,m)-st.mse(:,a);
    dm=mean(dif);se=std(dif)/sqrt(st.M);
    paired(end+1,:)={st.Signal,st.SNR,names{a},names{m},...
       dm,se,dm-1.96*se,dm+1.96*se,100*dm/mu(m)}; %#ok<AGROW>
   end
  end
 end
end
T=cell2table(rows,'VariableNames',{'Signal','SNR','Method','M','AMSE','SD_MSE','SE_AMSE'});
T.N=repmat(cfg.N,height(T),1); T.RMSE=sqrt(T.AMSE);
P=cell2table(paired,'VariableNames',{'Signal','SNR','SCOPE','Other','OtherMinusSCOPE','PairedSE','CI_L','CI_U','ImprovementPct'});
writetable(T,fullfile(outDir,'Experiment3_AMSE_all.csv'));
writetable(P,fullfile(outDir,'Experiment3_PairedDifferences.csv'));
Q=table(names(:),mean(scores,1)','VariableNames',{'Method','MeanRelativeAMSE_toHard'});
Q=sortrows(Q,'MeanRelativeAMSE_toHard');
writetable(Q,fullfile(outDir,'Experiment3_AggregateRanking.csv'));
if cfg.plot
 fig=figure('Visible','off','Color','w','Position',[100 100 1100 760]);
 ix=1:24;[v,order]=sort(Q.MeanRelativeAMSE_toHard(ismember(Q.Method,names(ix))));
 qnames=Q.Method(ismember(Q.Method,names(ix)));
 barh(v);set(gca,'YTick',1:numel(v),'YTickLabel',qnames(order),'YDir','reverse','FontSize',9);
 xlabel('Mean AMSE / universal hard AMSE across 12 conditions');
 title('Experiment III: aggregate relative performance (preliminary implementation audit)');
 grid on;print(fig,fullfile(plotDir,'E3_AggregateRanking.png'),'-dpng','-r180');
 print(fig,fullfile(plotDir,'E3_AggregateRanking.pdf'),'-dpdf','-painters');close(fig);
end
fid=fopen(fullfile(outDir,'Experiment3_README.txt'),'w');
fprintf(fid,'Independent comparison, M=%d per signal-SNR; 4 SCOPE+%d comparators+hard/soft.\n',cfg.M,numel(names)-6);
fprintf(fid,'26 methods: #5 BAMS author Cscale=2 gamma=2 depth J-3; #23 Chang-Yu-Vetterli context 1D.\n');
fprintf(fid,'Frozen F/k/r: Normal fixed r=0,1 k=5; Normal SCOPE SURE; Normal N-SCOPE SURE.\n');
fprintf(fid,'Section 4 calibration B_N=%g; first shrunk level j0=%d.\n',cfg.B_N,cfg.j0);
fprintf(fid,'Performance uses paired noisy realizations.\n');
fprintf(fid,'SNR=std(f,1)/std(epsilon) with iid epsilon N(0,1); std(epsilon)=1.\n');
fprintf(fid,'The sampled signal is scaled to std(f,1)=SNR, divisor n.\n');
fprintf(fid,'New rounded-B_N output directory; no old E3 checkpoint migration.\n');
fprintf(fid,'METHOD_IMPLEMENTATION_NOTES.csv identifies provisional 1D proxy methods.\n');
fprintf(fid,'These are not eligible to be described as exact original algorithms without further validation.\n');
fclose(fid);
end
 
function plot_condition(st,plotDir)
mu=mean(st.mse,1);idx=1:24;nI=numel(idx);
[val,ord]=sort(mu(idx));nm=st.names(idx);nm=nm(ord);
fig=figure('Visible','off','Color','w','Position',[100 100 1150 760]);
barh(val);set(gca,'YTick',1:nI,'YTickLabel',nm,'YDir','reverse','FontSize',9);
hold on;
plot([mu(25) mu(25)],[0.4 nI+0.6],'k--','LineWidth',1.3);
plot([mu(26) mu(26)],[0.4 nI+0.6],'k:','LineWidth',1.3);
xlim([0 1.08*max([mu(idx),mu(25),mu(26)])]); grid on;
xlabel('AMSE');title(sprintf('%s, SNR (SD ratio) %g, M=%d (unvalidated methods flagged in audit)',...
 st.Signal,st.SNR,st.M));
prefix=fullfile(plotDir,sprintf('E3_AMSE_%s_SNR%d',st.Signal,st.SNR));
print(fig,[prefix '.png'],'-dpng','-r160');
print(fig,[prefix '.pdf'],'-dpdf','-painters');close(fig);
% Paired-replication distribution of MSEs: four SCOPE and the six
% best-performing *plotted* benchmarks (figure selection only).
compIds=5:24;
[~,bestSix]=sort(mu(compIds)); ids=[1:4,compIds(bestSix(1:6))];
fig=figure('Visible','off','Color','w','Position',[100 100 1240 680]);
hold on;
for ii=1:numel(ids)
 v=sort(st.mse(:,ids(ii)));
 qq=empirical_quants(v,[.05 .25 .50 .75 .95]);
 plot([ii ii],[qq(1) qq(5)],'-','Color',[.3 .3 .3]);
 rectangle('Position',[ii-.29 qq(2) .58 max(qq(4)-qq(2),1e-9)], ...
    'EdgeColor',[.2 .2 .2],'FaceColor',[.88 .88 .88]);
 plot([ii-.29 ii+.29],[qq(3) qq(3)],'k-','LineWidth',1.4);
end
plot([.5 numel(ids)+.5],[mu(25) mu(25)],'k--');
plot([.5 numel(ids)+.5],[mu(26) mu(26)],'k:');
set(gca,'XTick',1:numel(ids),'XTickLabel',st.names(ids), ...
    'XTickLabelRotation',37,'FontSize',9);
ylabel('Replication MSE');title(sprintf('%s, SNR %g: 4 SCOPE plus 6 top plotted benchmarks',st.Signal,st.SNR));
box on;grid on;xlim([.5 numel(ids)+.5]);
base=fullfile(plotDir,sprintf('E3_MSEBoxes_%s_SNR%d',st.Signal,st.SNR));
print(fig,[base '.png'],'-dpng','-r160');
print(fig,[base '.pdf'],'-dpdf','-painters');close(fig);
end
 
function q=empirical_quants(x,p)
% No Statistics Toolbox. Empirical quantiles with linear interpolation.
x=sort(x(:));n=numel(x);z=1+(n-1)*p;
lo=floor(z);hi=ceil(z);ww=z-lo;
q=x(lo).*(1-ww)+x(hi).*ww;
end
 
 
function estimates=E3_legacy_methods(d,sigma,y,h,W)
% E3_LEGACY_METHODS: Eight established wavelet baseline implementations.
% Retains the computational routines/parameter choices supplied previously
% in compare_scope_four_signals.m. Columns correspond to published-table
% methods [1 2 3 4 5 6 9 13]. The routines for image-inspired methods
% are 1D wavelet versions and should be audited against the cited articles.
% Inputs: d=W*y, sigma MAD estimate, h low-pass filter, W orthonormal.
N=numel(d); J=round(log2(N)); L=3;
assert(numel(y)==N && 2^J==N);
D=cell(J+1,1);
for j=L:J-1
    D{j+1}=d((2^j+1):2^(j+1)).';
end
estimates=zeros(N,8);
% 1: BAMS (cached orthogonal matrix and user-supplied Bayes rule)
% BAMS.m interprets L as decomposition DEPTH, not coarsest j0.  For the
% common 8-scaling-coefficient design, L_BAMS=J-j0=J-L (L here is j0=3).
% Passing L=3 previously left detail levels 3,...,J-4 untouched.
estimates(:,1)=scope_BAMS_author_coefficients(d,J,J-L,2,2);
% 2: Johnstone-Silverman EBayes posterior median
z=d;
for j=L:J-1
  ind=(2^j+1):2^(j+1);
  z(ind)=ebayes_laplace_median_level(D{j+1},sigma,0.5,false,logspace(-1,log10(3),31));
end
estimates(:,2)=z;
% Quadrature used in original raised-cosine and Epanechnikov rules.
persistent glx glw
if isempty(glx), [glx,glw]=gauss_legendre_rule(64); end
% 3: Raised-cosine Bayes
z=d; tau=max(abs(d(9:end))); tau=max(tau,1e-8);
for j=L:J-1
  ind=(2^j+1):2^(j+1);
  alpha=min(max(1-1/(j-L+1)^2,0),1-1e-12);
  z(ind)=raised_cosine_bayes_level(D{j+1},sigma,alpha,tau,glx,glw);
end
estimates(:,3)=z;
% 4: Epanechnikov Bayes
z=d; sFine=max(std(D{J},0),1e-8);
lam=1/sFine^2+(1/2)*exp(-sFine/2);
for j=L:J-1
  ind=(2^j+1):2^(j+1);
  alpha=min(max(1-1/(j-L+1)^2,0),1-1e-12);
  beta=max(abs(D{j+1})); beta=max(beta,1e-8);
  z(ind)=epanechnikov_bayes_level(D{j+1},alpha,beta,lam,glx,glw);
end
estimates(:,4)=z;
% 5: Smooth SCAD -- AUTHOR VARIANT, not exact published algorithm.
% Levelwise joint SURE selection over a=[2,2.5,3,3.5], 121 lambda
% grid steps up to sigma*sqrt(2log N); postselection floor lambda>=2sigma.
% This floor is NOT SURE-minimizing, so do not label as plain Smooth SCAD.
z=d;
for j=L:J-1
  ind=(2^j+1):2^(j+1);
  z(ind)=smooth_scad_sure_level(D{j+1},sigma,[2 2.5 3 3.5],121,2,N);
end
estimates(:,5)=z;
% 6: Hybrid SureShrink (Donoho--Johnstone sparsity fallback)
z=d;
for j=L:J-1
  ind=(2^j+1):2^(j+1);
  z(ind)=sureshrink_level(D{j+1},sigma);
end
estimates(:,6)=z;
% 9: BayesShrink
z=d;
for j=L:J-1
  ind=(2^j+1):2^(j+1);
  z(ind)=bayesshrink_level(D{j+1},sigma);
end
estimates(:,7)=z;
% 13: Sendur-Selesnick bivariate shrinkage (1D parent-child adaptation).
% The original supplied routine uses j=L-1 as the coarsest detail;
% this retains j=3,...,9 with an imputed zero parent for the first level.
z=d;
for j=L:J-1
    ind=(2^j+1):2^(j+1); dj=d(ind).';
    if j==L, p=zeros(size(dj));
    else
       par=d((2^(j-1)+1):2^j).';
       p=repelem(par,2);
    end
    v=max(periodic_local_mean(dj.^2,7)-sigma^2,1e-10);
    s=sqrt(max(v,1e-10));
    q=sqrt(dj.^2+p.^2);
    T=sqrt(3)*sigma^2./s;
    z(ind)=dj.*max(1-T./max(q,realmin),0);
end
estimates(:,8)=z;
assert(all(isfinite(estimates(:))),'Legacy estimator yielded nonfinite value');
end
 
function [base,h,g]=benchmark_signal_filter(signal,N)
% From the four original MATLAB comparison files.
t=(1:N)/(N+1);
pos=[.10,.13,.15,.23,.25,.40,.44,.65,.76,.78,.81];
switch signal
    case 'Doppler'
        base=sqrt(t.*(1-t)).*sin(2*pi*1.05./(t+0.05));
        h=[-.0757657147892733 -.0296355276459985 .4976186676320150 ...
           .8037387518059160 .2978577956052770 -.0992195435768472 ...
           -.0126039672620378 .0322231006040427];
    case 'HeaviSine'
        base=4*sin(4*pi*t)-sign(t-.30)-sign(.72-t);
        h=[-.0757657147892733 -.0296355276459985 .4976186676320150 ...
           .8037387518059160 .2978577956052770 -.0992195435768472 ...
           -.0126039672620378 .0322231006040427];
    case 'Bumps'
        ht=[4 5 3 4 5 4.2 2.1 4.3 3.1 5.1 4.2];
        w=[.005 .005 .006 .010 .010 .030 .010 .010 .005 .008 .005];
        base=zeros(size(t));
        for z=1:numel(pos)
            base=base+ht(z)./(1+abs((t-pos(z))/w(z))).^4;
        end
        h=[.0352262918857095 -.0854412738820267 -.1350110200102550 ...
           .4598775021184915 .8068915093110928 .3326705529500826];
    case 'Blocks'
        ht=[4 -5 3 -4 5 -4.2 2.1 4.3 -3.1 2.1 -4.2];
        base=zeros(size(t));
        for z=1:numel(pos)
            base=base+ht(z).*(t>=pos(z));
        end
        h=[1 1]/sqrt(2);
    otherwise
        error('Unknown signal %s',signal);
end
% The original scripts use the following high-pass orientation for these
% three wavelet families. For sym4/db3 it differs by an overall sign from
% the common QMF convention, which makes no difference for odd shrinkers.
if strcmp(signal,'Blocks')
    g=[1 -1]/sqrt(2);
else
    g=-((-1).^(0:length(h)-1)).*fliplr(h);
end
end
 
function [Df,Ds]=scope_both(D,L,J,sigmaHat,fi,alpha,kFixed,kGrid,bGrid,q,tref)
% Both SCOPE calibrations on exactly the same observed coefficients.
Df=D; Ds=D;
% Fixed: b=1,k=5 and distribution-matched half-attenuation.
for j=L:J-1
    d=D{j+1};
    lam=q(kGrid==kFixed)/(sigmaHat*tref)*2^(alpha*(J-1-j));
    [H,~]=scope_Fstar_and_deriv(lam*abs(d),fi);
    Df{j+1}=d.*H.^kFixed;
end
best=Inf; bestB=1; bestK=kFixed;
for ki=1:numel(kGrid)
    k=kGrid(ki);
    for bi=1:numel(bGrid)
        b=bGrid(bi); risk=0;
        for j=L:J-1
            d=D{j+1};
            lam=q(ki)/(b*sigmaHat*tref)*2^(alpha*(J-1-j));
            u=lam*abs(d);
            [H,h]=scope_Fstar_and_deriv(u,fi);
            shrink=H.^k;
            delta=d.*shrink;
            div=shrink+k*u.*H.^(k-1).*h;
            risk=risk+sum((delta-d).^2+2*sigmaHat^2*div-sigmaHat^2);
        end
        if risk<best, best=risk; bestB=b; bestK=k; end
    end
end
qStar=q(kGrid==bestK);
for j=L:J-1
    d=D{j+1};
    lam=qStar/(bestB*sigmaHat*tref)*2^(alpha*(J-1-j));
    [H,~]=scope_Fstar_and_deriv(lam*abs(d),fi);
    Ds{j+1}=d.*H.^bestK;
end
end
 
function q=scope_half_quantile(fi,k)
r=2^(-1/k);
switch fi
    case 1, q=2*atanh(r);                  % logistic
    case 2, q=sqrt(2)*erfinv(r);          % normal
    case 3, q=-log1p(-r);                 % unit Laplace
    case 4, q=4/pi*atanh(tan(pi*r/4));    % hyperbolic secant
    case 5, q=r;                          % uniform on [-1,1]
    case 6, q=tan(pi*r/2);                % standard Cauchy
end
end
 
function [H,h]=scope_Fstar_and_deriv(u,fi)
switch fi
    case 1, H=tanh(u/2); h=(1-H.^2)/2;
    case 2, H=erf(u/sqrt(2)); h=sqrt(2/pi)*exp(-u.^2/2);
    case 3, H=-expm1(-u); h=exp(-u);
    case 4, H=(4/pi)*atan(tanh(pi*u/4)); h=1./cosh(pi*u/2);
    case 5, H=min(u,1); h=double(u<1);
    case 6, H=(2/pi)*atan(u); h=(2/pi)./(1+u.^2);
end
end
 
function v=caucau_mean_fast(d,sigma,j,C,alpha,tgrid)
% Same Gaussian CauCau posterior mean as original, accelerated by a
% precomputed interpolation table for standardized |d|/sigma at each level.
% tau/sigma is fixed by the prescribed schedule; all other parameters stay
% unchanged. 'pchip' preserves shape; direct exact integral for z>30.
persistent jCache uCache mCache
if isempty(jCache), jCache=[]; uCache={}; mCache={}; end
ind=find(jCache==j,1);
if isempty(ind)
    u=0:.02:30;
    out=zeros(size(u));
    tauStd=C*2^(-alpha*j);
    for b=1:100:numel(u)
        ix=b:min(b+99,numel(u));
        out(ix)=caucau_mean_vectorized(u(ix),tauStd,1,tgrid);
    end
    jCache(end+1)=j;
    uCache{end+1}=u;
    mCache{end+1}=out;
    ind=numel(jCache);
end
z=abs(d)/sigma;
w=interp1(uCache{ind},mCache{ind},min(z,30),'pchip');
if any(z(:)>30)
    ind2=z>30;
    w(ind2)=caucau_mean_vectorized(z(ind2),C*2^(-alpha*j),1,tgrid);
end
v=sigma*sign(d).*w;
end
 
function benchmark_boxplot(MSE,names,universal,signal,SNR,outdir)
% Tukey boxplots, no Statistics Toolbox; red dotted Universal reference.
[~,K]=size(MSE);
hf=figure('Visible','off','Color','w','Position',[60 80 1650 630]);
hold on;
for k=1:K
    xs=sort(MSE(:,k)); n=numel(xs);
    q1=lin_quantile(xs,.25); md=lin_quantile(xs,.5);
    q3=lin_quantile(xs,.75); iqr0=q3-q1;
    inside=xs(xs>=q1-1.5*iqr0 & xs<=q3+1.5*iqr0);
    if isempty(inside), lo=q1; hi=q3;
    else, lo=inside(1); hi=inside(end);
    end
    plot([k k],[lo q1],'k-',[k k],[q3 hi],'k-');
    plot([k-.1 k+.1],[lo lo],'k-',[k-.1 k+.1],[hi hi],'k-');
    rectangle('Position',[k-.25 q1 .5 max(q3-q1,eps)],'EdgeColor','k');
    plot([k-.25 k+.25],[md md],'k-');
    plot([k-.15 k+.15],[mean(xs) mean(xs)],'r-','LineWidth',1.3);
    out=xs(xs<q1-1.5*iqr0 | xs>q3+1.5*iqr0);
    if ~isempty(out), plot(k*ones(size(out)),out,'k.','MarkerSize',4); end
end
plot([.5 K+.5],[universal universal],'r:','LineWidth',2);
xlim([.5 K+.5]); xticks(1:K);
xticklabels(arrayfun(@num2str,1:K,'UniformOutput',false));
xlabel('Method number (see methods key CSV)'); ylabel('MSE');
title(sprintf('%s, SNR (SD ratio)=%g, M=%d',signal,SNR,size(MSE,1)));
set(gca,'FontSize',11,'Box','on'); hold off;
stem=sprintf('%s_SNR%g',signal,SNR);
print(hf,fullfile(outdir,[stem '_boxplot.png']),'-dpng','-r220');
print(hf,fullfile(outdir,[stem '_boxplot.pdf']),'-dpdf','-painters');
close(hf);
key=table((1:K)',names','VariableNames',{'Number','Method'});
writetable(key,fullfile(outdir,'method_key.csv'));
end
 
function q=lin_quantile(x,p)
n=numel(x); loc=1+(n-1)*p; a=floor(loc); b=ceil(loc);
q=x(a)+(loc-a)*(x(b)-x(a));
end
 
% MATLAB functions copied from the user-provided 17-method script.
function txt = format_duration(sec)
% Compact hh:mm:ss display for Monte Carlo progress.
sec = max(0,round(sec));
h = floor(sec/3600);
m = floor(mod(sec,3600)/60);
s = mod(sec,60);
txt = sprintf('%02d:%02d:%02d',h,m,s);
end
 
 
function sigmaHat = robust_sigma_gaussian(dFine)
 
x = dFine(:);
medx = median(x);
madx = median(abs(x-medx));
 
z075 = 0.6744897501960817;
 
sigmaHat = madx/z075;
sigmaHat = max(sigmaHat,1e-8);
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STANDARDIZED Cau/Cau SPoM BIRTH GEOMETRY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function q = empirical_quantile_sorted(xs,p)
% Quantile of an already sorted vector using linear interpolation at
% r = 1 + (n-1)p.  No Statistics Toolbox required.
 
xs = xs(:);
n = numel(xs);
 
if n==0
    q = NaN;
    return;
end
 
if p<=0
    q = xs(1);
    return;
elseif p>=1
    q = xs(end);
    return;
end
 
r = 1 + (n-1)*p;
i1 = floor(r);
i2 = ceil(r);
 
if i1==i2
    q = xs(i1);
else
    w = r-i1;
    q = (1-w)*xs(i1) + w*xs(i2);
end
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAMBDA-NEIGHBORHOOD WAVELET SHRINKAGE (LNWS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BAYESSHRINK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [dhat,T] = bayesshrink_level(d,sigma)
% 1-D levelwise specialization of Chang-Yu-Vetterli BayesShrink.
%
% sigma_Y^2 = empirical second moment of noisy coefficients.
% sigma_X   = sqrt(max(sigma_Y^2-sigma^2,0)).
% T         = sigma^2/sigma_X.
% Apply soft thresholding.
 
sigmaY2 = mean(d(:).^2);
 
sigmaX = sqrt(max(sigmaY2-sigma^2,0));
 
if sigmaX <= sqrt(eps)*max(1,sigma)
    T = Inf;
else
    T = sigma^2/sigmaX;
end
 
dhat = soft_threshold(d,T);
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SENDUR-SELESNICK BIVARIATE SHRINKAGE: 1-D ORTHOGONAL-DWT SPECIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SENDUR-SELESNICK BIVARIATE SHRINKAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [Dhat,Tmean,Smean] = sendur_selesnick_1d( ...
                                      D,L,J,sigma,windowSize)
% Orthogonal-DWT 1-D specialization of Sendur-Selesnick BiShrink with
% local variance estimation.
%
% For child d and its parent p:
%
%       R = sqrt(d^2+p^2),
%       d_hat = d (R-T)_+ / R,
%       T = sqrt(3) sigma^2 / sigma_X.
%
% sigma_X is estimated locally from neighboring squared CHILD
% coefficients.  The original separable-DWT image implementation uses
% window size 7.  Here we use the corresponding 1-D periodic window.
 
Dhat = D;
 
Tmean = NaN(1,J);
Smean = NaN(1,J);
 
for j = J-1:-1:L
 
    child  = D{j+1}(:).';
    parent = D{j}(:).';
 
    if numel(child) ~= 2*numel(parent)
        error('Sendur-Selesnick parent/child size mismatch at level j=%d.',j);
    end
 
    % Each coarser parent is associated with its two dyadic children.
    parentExpanded = repelem(parent,2);
 
    % Local observed energy using periodic neighborhood.
    Wsig = periodic_local_mean(child.^2,windowSize);
 
    % Local signal standard deviation.
    Ssig = sqrt(max(Wsig-sigma^2,eps));
 
    T = sqrt(3)*sigma^2 ./ Ssig;
 
    R = sqrt(child.^2 + parentExpanded.^2);
 
    Rminus = max(R-T,0);
 
    fac = zeros(size(R));
 
    nz = (R > realmin);
    fac(nz) = Rminus(nz)./R(nz);
 
    out = child .* fac;
 
    Dhat{j+1} = reshape(out,size(D{j+1}));
 
    Tmean(j+1) = mean(T);
    Smean(j+1) = mean(Ssig);
 
end
 
% D{L} and the scaling coefficients at L-1 are untouched.  Together they
% reconstruct the same coarse V_L component used by the other methods.
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function z = periodic_local_mean(x,w)
% Centered periodic moving average with odd window length w.
% No toolbox required.
 
x = x(:).';
 
r = (w-1)/2;
 
z = zeros(size(x));
 
for k = -r:r
    z = z + circshift(x,[0 k]);
end
 
z = z/w;
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERIODIC ORTHONORMAL DWT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function out = soft_threshold(x,T)
 
out = sign(x).*max(abs(x)-T,0);
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function p = Phi0(x)
% Standard normal CDF without Statistics Toolbox.
 
p = 0.5*erfc(-x/sqrt(2));
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function x = PhiInv0(p)
% Standard normal inverse CDF without Statistics Toolbox.
 
p = min(max(p,1e-15),1-1e-15);
 
x = -sqrt(2)*erfcinv(2*p);
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function y = phi0(x)
% Standard normal density.
 
y = exp(-0.5*x.^2)/sqrt(2*pi);
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EBayes-MED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [dhat,wHat,aHat] = ebayes_laplace_median_level( ...
              d,sigma,aDefault,estimateA,aGrid)
% Johnstone-Silverman empirical Bayes posterior median.
%
% Standardize z=d/sigma and use the Laplace slab
%
%   gamma(theta;a) = (a/2) exp(-a |theta|).
%
% w is the prior probability of a NONZERO coefficient.
%
% If estimateA=false, a=aDefault and w is MML.
% If estimateA=true, a is chosen by profile MML over aGrid, with w MML
% at every candidate.  The latter is an intentionally transparent
% numerical profile-MML implementation.
 
z = d(:).'/sigma;
 
if ~estimateA
 
    aHat = aDefault;
    [wHat,~] = ebayes_fit_w(z,aHat);
 
else
 
    bestLL = -Inf;
    wHat   = 0.5;
    aHat   = aDefault;
 
    for ia = 1:length(aGrid)
 
        aa = aGrid(ia);
 
        [ww,ll] = ebayes_fit_w(z,aa);
 
        if ll > bestLL
            bestLL = ll;
            wHat   = ww;
            aHat   = aa;
        end
 
    end
 
end
 
muhat = ebayes_postmed_laplace(z,wHat,aHat);
 
dhat = sigma*muhat;
dhat = reshape(dhat,size(d));
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [wHat,llHat] = ebayes_fit_w(z,a)
% Profile MML estimate of w for fixed Laplace rate a.
% Constraint: posterior-median threshold <= sqrt(2 log n).
 
z = z(:).';
n = length(z);
 
x = abs(z);
 
[ph,g,~,~] = laplace_marginal_parts(x,a);
 
% Lower admissible w determined by universal threshold
tU = sqrt(2*log(max(n,2)));
 
[phT,~,ApT,AmT] = laplace_marginal_parts(tU,a);
 
denT = phT + (ApT-AmT);
 
if denT <= 0 || ~isfinite(denT)
    wMin = 1e-8;
else
    wMin = phT/denT;
end
 
wMin = min(max(wMin,1e-8),1-1e-8);
 
% Score S(w) = sum (g-phi)/((1-w)phi+w g), decreasing in w.
score = @(w) sum((g-ph)./max((1-w)*ph+w*g,realmin));
 
sLo = score(wMin);
sHi = score(1-1e-10);
 
if sLo <= 0
 
    wHat = wMin;
 
elseif sHi >= 0
 
    wHat = 1-1e-10;
 
else
 
    lo = wMin;
    hi = 1-1e-10;
 
    for it = 1:60
 
        mid = 0.5*(lo+hi);
        sm  = score(mid);
 
        if sm > 0
            lo = mid;
        else
            hi = mid;
        end
 
    end
 
    wHat = 0.5*(lo+hi);
 
end
 
marg = (1-wHat)*ph + wHat*g;
 
llHat = sum(log(max(marg,realmin)));
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function mu = ebayes_postmed_laplace(z,w,a)
% Explicit posterior median under spike + Laplace slab.
 
sz = size(z);
 
z = z(:).';
 
sgn = sign(z);
x   = abs(z);
 
[ph,g,Ap,~] = laplace_marginal_parts(x,a);
 
marg = (1-w)*ph + w*g;
wpost = w*g./max(marg,realmin);
 
% Conditional slab posterior upper-tail probability at zero
Ftail0 = Ap./max(g,realmin);
 
isZero = (wpost.*Ftail0 <= 0.5);
 
mabs = zeros(size(x));
 
idx = find(~isZero & x>0);
 
if ~isempty(idx)
 
    xx = x(idx);
    wp = wpost(idx);
 
    % D = exp(-a x) Phi(x-a) + exp(a x) Phi(-x-a)
    D = exp(-a*xx).*Phi0(xx-a) + ...
        exp( a*xx).*Phi0(-xx-a);
 
    % From Ftilde_1(m|x) = 1/(2 wpost):
    %
    % exp(-a x) Phi(x-a-m)/D = 1/(2 wpost).
    u = D.*exp(a*xx)./(2*wp);
 
    u = min(max(u,1e-14),1-1e-14);
 
    mm = xx - a - PhiInv0(u);
 
    mabs(idx) = max(mm,0);
 
end
 
mu = sgn.*mabs;
mu = reshape(mu,sz);
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [ph,g,Ap,Am] = laplace_marginal_parts(x,a)
% For x>=0:
%
% Ap = integral_{theta>0} Laplace(theta;a) phi(x-theta) dtheta
% Am = integral_{theta<0} Laplace(theta;a) phi(x-theta) dtheta
% g  = Ap+Am.
 
x = abs(x);
 
ph = phi0(x);
 
Ap = 0.5*a*exp(0.5*a^2-a*x).*Phi0(x-a);
Am = 0.5*a*exp(0.5*a^2+a*x).*Phi0(-x-a);
 
% Numerical guard against Inf*0 in extreme tails
bad = ~isfinite(Ap);
if any(bad)
    Ap(bad) = 0;
end
 
bad = ~isfinite(Am);
if any(bad)
    Am(bad) = 0;
end
 
g = Ap+Am;
g = max(g,realmin);
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RAISED-COSINE BAYES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function dhat = raised_cosine_bayes_level( ...
                   d,sigma,alpha,tau,xGL,wGL)
% Posterior mean under
%
% pi(theta)=alpha delta_0 +(1-alpha) g_RC(theta;tau),
%
% g_RC(theta;tau) =
%   [1 + cos(pi theta/tau)]/(2 tau), |theta|<tau.
%
% Bounded-support integrals are evaluated by vectorized Gauss-Legendre
% quadrature; this is mathematically the same posterior mean as the paper's
% analytic-plus-integral representation.
 
sz = size(d);
d  = d(:).';
 
if tau <= 1e-12
    dhat = zeros(sz);
    return;
end
 
theta = tau*xGL(:);
qw    = tau*wGL(:);
 
gprior = (1 + cos(pi*theta/tau))/(2*tau);
 
% nq x nd matrix
E = exp(-0.5*((theta-d)/sigma).^2);
 
slabMarg = sum((qw.*gprior).*E,1);
slabNum  = sum((qw.*gprior.*theta).*E,1);
 
zeroLik = exp(-0.5*(d/sigma).^2);
 
den = alpha*zeroLik + (1-alpha)*slabMarg;
 
out = (1-alpha)*slabNum./max(den,realmin);
 
dhat = reshape(out,sz);
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPANECHNIKOV BAYES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function dhat = epanechnikov_bayes_level( ...
                   d,alpha,beta,lambda,xGL,wGL)
% Posterior mean under spike + Epanechnikov slab.
%
% The Epanechnikov paper places Exp(lambda) on sigma^2 and integrates it
% out. The resulting likelihood in d-theta is double exponential with rate
%
%        r = sqrt(2 lambda).
%
% Constants common to atom and slab cancel in the posterior mean.
 
sz = size(d);
d  = d(:).';
 
if beta <= 1e-12
    dhat = zeros(sz);
    return;
end
 
theta = beta*xGL(:);
qw    = beta*wGL(:);
 
gprior = (3/(4*beta))*(1-(theta/beta).^2);
 
rate = sqrt(2*lambda);
 
E = exp(-rate*abs(theta-d));
 
slabMarg = sum((qw.*gprior).*E,1);
slabNum  = sum((qw.*gprior.*theta).*E,1);
 
zeroLik = exp(-rate*abs(d));
 
den = alpha*zeroLik + (1-alpha)*slabMarg;
 
out = (1-alpha)*slabNum./max(den,realmin);
 
dhat = reshape(out,sz);
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMOOTH SCAD-SURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [dhat,lambdaHat,aHat,lambdaSURE] = ...
    smooth_scad_sure_level(d,sigma,aGrid,nLambda,cFloor,Nglobal)
% Adaptive/stabilized levelwise SURE for the smooth raised-cosine SCAD rule.
%
% Jointly select a and the raw SURE threshold over
%
%       a in aGrid,   0 <= lambda <= lambda_U,
%
% where
%
%       lambda_U = sigma*sqrt(2*log(Nglobal)).
%
% The selected threshold is stabilized as
%
%       lambdaHat = min(lambda_U,max(lambdaSURE,cFloor*sigma)).
%
% The benchmark settings use
%
%       aGrid  = [2 2.5 3 3.5],   cFloor = 2.
 
sz = size(d);
d  = d(:).';
 
n = length(d);
 
if n == 0
    dhat = reshape(d,sz);
    lambdaHat = 0;
    aHat = aGrid(1);
    lambdaSURE = 0;
    return;
end
 
if nargin < 6 || isempty(Nglobal)
    Nglobal = n;
end
 
if nargin < 5 || isempty(cFloor)
    cFloor = 2;
end
 
aGrid = aGrid(:).';
 
if any(aGrid <= 1)
    error('Smooth SCAD requires every a candidate to exceed 1.');
end
 
lambdaU = sigma*sqrt(2*log(max(Nglobal,2)));
 
lamGrid = [0,linspace(lambdaU/nLambda,lambdaU,nLambda)];
 
r = abs(d);
 
bestSURE   = Inf;
aHat       = aGrid(1);
lambdaSURE = 0;
 
for ia = 1:numel(aGrid)
 
    a = aGrid(ia);
 
    % lambda=0 => identity estimator; SURE = n sigma^2
    if n*sigma^2 < bestSURE
        bestSURE   = n*sigma^2;
        aHat       = a;
        lambdaSURE = 0;
    end
 
    for il = 2:length(lamGrid)
 
        lam = lamGrid(il);
 
        low = (r <= lam);
        mid = (r > lam & r < a*lam);
        hi  = (r >= a*lam);
 
        F = zeros(size(d));
 
        % rho=0
        F(low) = d(low).^2 - sigma^2;
 
        % rho=d-h
        if any(mid)
 
            s = (r(mid)-lam)/((a-1)*lam);
 
            hmag = lam*cos(0.5*pi*s).^2;
 
            % derivative dh/dd is even and equals:
            hp = -pi/(2*(a-1))*sin(pi*s);
 
            F(mid) = hmag.^2 + sigma^2 - 2*sigma^2*hp;
 
        end
 
        % rho=d
        F(hi) = sigma^2;
 
        sureVal = sum(F);
 
        if sureVal < bestSURE
            bestSURE   = sureVal;
            aHat       = a;
            lambdaSURE = lam;
        end
 
    end
end
 
lambdaHat = min(lambdaU,max(lambdaSURE,cFloor*sigma));
 
dhat = smooth_scad_apply(d,lambdaHat,aHat);
dhat = reshape(dhat,sz);
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function out = smooth_scad_apply(d,lambda,a)
 
if lambda <= 0
    out = d;
    return;
end
 
r = abs(d);
 
out = zeros(size(d));
 
mid = (r > lambda & r < a*lambda);
hi  = (r >= a*lambda);
 
if any(mid)
 
    s = (r(mid)-lambda)/((a-1)*lambda);
 
    hmag = lambda*cos(0.5*pi*s).^2;
 
    out(mid) = d(mid) - sign(d(mid)).*hmag;
 
end
 
out(hi) = d(hi);
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NCPmean-2: WANG-WOOD EMPIRICAL-BAYES BLOCK SHRINKAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [Dhat,pHat,lambdaHat,sigmaHat] = ncpmean2_shrink( ...
                     D,L,J,sigmaKnown,useKnownSigma,xGL,wGL)
% Wang-Wood NCPmean with block size m=2.
%
% Reference construction:
%   z_B = ||d_B||^2
% is modeled through the noncentral chi-square representation.  The prior
% for the true block energy has a point mass at zero with probability p
% plus the Wang-Wood power prior.  The posterior mean Pmean(z) is used.
%
% Hyperparameters:
%   * sigma^2 is estimated from finest-level block energies unless
%     useKnownSigma=true;
%   * one global power-prior shape lambda is obtained by the
%     Wang-Wood quick moment/median construction;
%   * p_j is estimated separately at each level by marginal MLE.
%
% The estimated block vector is
%
%       d_B_hat = sqrt(Pmean(z_B)/z_B) d_B.
%
% No cap at one is imposed: the published posterior-mean rule is not a
% strict shrinkage rule, although it has bounded-shrinkage behavior.
 
m = 2;
 
% Gauss-Legendre rule on (0,1) for the C_eta integrals.
u  = 0.5*(xGL(:) + 1);
wu = 0.5*wGL(:);
 
% -------------------------------------------------------------------------
% Block energies at every processed level.
% -------------------------------------------------------------------------
Z = cell(1,J);
allz = [];
 
for j = L:J-1
 
    d = D{j+1}(:).';
 
    if mod(numel(d),m) ~= 0
        error('NCPmean-2 requires an even number of coefficients at each level.');
    end
 
    db = reshape(d,m,[]);
    z  = sum(db.^2,1);
 
    Z{j+1} = z;
    allz = [allz,z]; %#ok<AGROW>
 
end
 
% -------------------------------------------------------------------------
% Noise variance.  For m=2, median(chi-square_2) = 2 log 2 exactly.
% -------------------------------------------------------------------------
if useKnownSigma
 
    sigma2Hat = sigmaKnown^2;
 
else
 
    dFine = D{J}(:).';
    zFine = sum(reshape(dFine,m,[]).^2,1);
 
    chi2med_m2 = 2*log(2);
    sigma2Hat = median(zFine)/chi2med_m2;
 
    sigma2Hat = max(sigma2Hat,1e-12);
 
end
 
sigmaHat = sqrt(sigma2Hat);
 
% -------------------------------------------------------------------------
% Wang-Wood quick-and-dirty global power-prior shape.
%
% Preliminary zero-mass estimate:
%
%   p0 = #{z_i < 2 m sigma^2 log(Nblocks)} / Nblocks.
%
% Moment matching gives the central beta:
%
%   sum z_i = sigma^2 sum m_i + (1-p0) sum m_i / betaBar.
%
% Transform uBar = betaBar*sigma^2/(1+betaBar*sigma^2).  Under the
% power prior U~Beta(lambda+1,1), so choosing betaBar as the prior median
% gives
%
%   lambda = -log(2)/log(uBar) - 1.
% -------------------------------------------------------------------------
Nblocks = numel(allz);
 
if Nblocks == 0
    error('NCPmean-2 found no detail blocks.');
end
 
threshold0 = 2*m*sigma2Hat*log(max(Nblocks,2));
p0 = mean(allz < threshold0);
 
sumM = m*Nblocks;
excess = sum(allz) - sigma2Hat*sumM;
 
if excess <= 1e-12*sigma2Hat*sumM
    % Essentially no positive empirical excess energy.  The limiting
    % moment estimate betaBar is very large; use a finite numerical proxy.
    betaBar = 1e8/max(sigma2Hat,1e-12);
else
    betaBar = (1-p0)*sumM/excess;
    betaBar = max(betaBar,1e-12/max(sigma2Hat,1e-12));
end
 
uBar = betaBar*sigma2Hat/(1 + betaBar*sigma2Hat);
uBar = min(max(uBar,1e-10),1-1e-10);
 
lambdaHat = -log(2)/log(uBar) - 1;
lambdaHat = max(lambdaHat,-1+1e-8);
 
% -------------------------------------------------------------------------
% Levelwise marginal-MLE p_j and posterior-mean block shrinkage.
% -------------------------------------------------------------------------
Dhat = D;
pHat = NaN(1,J);
 
for j = L:J-1
 
    d = D{j+1}(:).';
    db = reshape(d,m,[]);
    z  = Z{j+1};
 
    [logR,A1,A2] = ncp_power_quantities( ...
                      z,sigma2Hat,lambdaHat,m,u,wu);
 
    pJ = ncp_fit_zero_mass_mle(logR);
    pHat(j+1) = pJ;
 
    % Posterior probability that the true block energy is zero.
    logitP = log(pJ) - log(1-pJ);
    pi0 = logistic_stable(logitP - logR);
 
    % Wang-Wood posterior mean of true block energy, Eq. (3.43):
    %
    % (1-pi0) [m*sigma^2 + z
    %          -(m*sigma^2+2z) A_{eta,1} + z A_{eta,2}].
    Pmean = (1-pi0) .* ...
            (m*sigma2Hat + z ...
             - (m*sigma2Hat + 2*z).*A1 ...
             + z.*A2);
 
    % Numerical guard only; theoretical posterior energy is nonnegative.
    Pmean = max(Pmean,0);
 
    fac = zeros(size(z));
    nz = (z > 1e-14);
    fac(nz) = sqrt(Pmean(nz)./z(nz));
 
    dbhat = db .* fac;
 
    Dhat{j+1} = reshape(dbhat,1,[]);
    Dhat{j+1} = reshape(Dhat{j+1},size(D{j+1}));
 
end
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [logR,A1,A2] = ncp_power_quantities(z,sigma2,lambda,m,u,wu)
% Quantities for the Wang-Wood power prior.
%
% Let r=z/(2 sigma^2), eta=m/2+lambda, and
%
%   C_eta(r) = integral_0^1 u^eta exp(-r u) du.
%
% Then
%
%   R(z) = f_1(z)/f_0(z)
%        = (lambda+1) exp(r) C_eta(r),
%
% and
%
%   A_{eta,1}=C_{eta+1}/C_eta,   A_{eta,2}=C_{eta+2}/C_eta.
%
% We evaluate these by stable vectorized Gauss-Legendre quadrature.
 
z = z(:).';
r = z/(2*sigma2);
eta = m/2 + lambda;
 
if eta <= -1
    error('NCP power-prior integral requires eta > -1.');
end
 
logu = log(u);
logwu = log(wu);
 
% nq x nb matrix
logBase = logwu + eta*logu - u*r;
 
mx = max(logBase,[],1);
W = exp(logBase - mx);
 
den = sum(W,1);
den = max(den,realmin);
 
A1 = sum(W.*u,1)./den;
A2 = sum(W.*(u.^2),1)./den;
 
logC = mx + log(den);
 
logR = log(lambda+1) + r + logC;
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function pHat = ncp_fit_zero_mass_mle(logR)
% Marginal MLE of p, the prior probability of zero block energy.
%
% Apart from the central-chi-square factor, the level likelihood is
%
%   prod_i [ p + (1-p) R_i ].
%
% The likelihood is concave in p.  fminbnd is used on the log scale
% expression, and the two boundary limits are checked explicitly.
 
logR = logR(:).';
 
obj = @(p) ncp_neg_profile_loglik(p,logR);
 
opts = optimset('Display','off','TolX',1e-9);
pInt = fminbnd(obj,1e-9,1-1e-9,opts);
 
cand = [1e-10,pInt,1-1e-10];
val  = [obj(cand(1)),obj(cand(2)),obj(cand(3))];
 
[~,ii] = min(val);
pHat = cand(ii);
 
pHat = min(max(pHat,1e-10),1-1e-10);
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function nll = ncp_neg_profile_loglik(p,logR)
 
a = log(p);
b = log(1-p) + logR;
 
mx = max(a,b);
ll = sum(mx + log(exp(a-mx) + exp(b-mx)));
 
nll = -ll;
 
if ~isfinite(nll)
    nll = realmax/1e100;
end
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function y = logistic_stable(x)
 
y = zeros(size(x));
 
ip = (x >= 0);
im = ~ip;
 
y(ip) = 1./(1+exp(-x(ip)));
 
ex = exp(x(im));
y(im) = ex./(1+ex);
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAMBDA-NEIGHBORHOOD WAVELET SHRINKAGE (LNWS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [Dhat,epsHat,bHat] = lnws_posterior_mean( ...
              D,L,J,sigma,bMin,bMax,nBgrid)
% Remenyi-Vidakovic Lambda-neighborhood empirical-Bayes shrinkage.
%
% The likelihood is exactly the marginal likelihood in Eq. (15) of the
% LNWS paper.  Instead of requiring Optimization Toolbox/fmincon, we profile
% out epsilon_j at each candidate b.  This maximizes the same likelihood:
%
%   for fixed b, each epsilon_j occurs only in level j and is optimized
%   independently; b is then optimized globally over the profile likelihood.
%
% Posterior-mean energy (Eq. 9) is used, followed by the coefficient
% shrinkage rule in Eq. (11).
 
s2 = sigma^2;
 
% -------------------------------------------------------------------------
% Form all Lambda-neighborhood energy statistics x_{j,l}.
% Child pairs are nonoverlapping: (1,2), (3,4), ...
% Their common parents are positions 1,2,...
% -------------------------------------------------------------------------
 
X = cell(1,J);
 
allx = [];
 
for j = L:J-1
 
    child  = D{j+1}(:).';
    parent = D{j}(:).';      % detail level j-1
 
    nchild = length(child);
    npair  = nchild/2;
 
    c1 = child(1:2:end);
    c2 = child(2:2:end);
 
    if length(parent) ~= npair
        error('LNWS parent/child indexing mismatch at level j=%d.',j);
    end
 
    x = (c1.^2 + c2.^2 + parent.^2)/s2;
 
    x = max(x,1e-14);
 
    X{j+1} = x;
 
    allx = [allx,x]; %#ok<AGROW>
 
end
 
% -------------------------------------------------------------------------
% Global b: paper suggests starting b0 = 1/mean(x).
% We use a broad log-grid and then refine around the best profile point.
% This avoids requiring Optimization Toolbox and is substantially more
% stable than a single unconstrained multidimensional optimization.
% -------------------------------------------------------------------------
 
b0 = 1/max(mean(allx),1e-12);
 
bLo = min(bMin,b0/100);
bHi = max(bMax,b0*100);
 
bLo = max(bLo,1e-10);
bHi = min(bHi,1e4);
 
lg = linspace(log(bLo),log(bHi),nBgrid);
 
prof = zeros(size(lg));
 
for ib = 1:length(lg)
    prof(ib) = lnws_profile_nll(exp(lg(ib)),X,L,J);
end
 
[~,imin] = min(prof);
 
ilo = max(1,imin-1);
ihi = min(length(lg),imin+1);
 
if ilo == ihi
 
    logbHat = lg(imin);
 
else
 
    obj = @(z) lnws_profile_nll(exp(z),X,L,J);
 
    opts = optimset('Display','off','TolX',1e-7);
 
    logbHat = fminbnd(obj,lg(ilo),lg(ihi),opts);
 
end
 
bHat = exp(logbHat);
 
% Final levelwise epsilon estimates at bHat
epsHat = NaN(1,J);
 
for j = L:J-1
 
    x = X{j+1};
 
    epsHat(j+1) = lnws_fit_epsilon(x,bHat);
 
end
 
% -------------------------------------------------------------------------
% Apply posterior-mean energy shrinkage.
%
% Only child coefficients are shrunk at level j.  Their parent is used to
% construct x but is left unchanged until its own level is processed.
% Detail level L-1 is never shrunk.
% -------------------------------------------------------------------------
 
Dhat = D;
 
for j = L:J-1
 
    child = D{j+1}(:).';
 
    x   = X{j+1};
    eps = epsHat(j+1);
 
    delta = lnws_delta_mean(x,eps,bHat);
 
    fac = sqrt(max(delta,0)./max(x,1e-14));
 
    % The published rule is exactly sqrt(delta/x); no artificial cap at 1.
    out = child;
 
    out(1:2:end) = fac .* child(1:2:end);
    out(2:2:end) = fac .* child(2:2:end);
 
    Dhat{j+1} = reshape(out,size(D{j+1}));
 
end
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function nll = lnws_profile_nll(b,X,L,J)
% Profile negative marginal log-likelihood for global rate b.
 
if ~(isfinite(b) && b>0)
 
    nll = realmax/1e100;
    return;
 
end
 
nll = 0;
 
for j = L:J-1
 
    x = X{j+1};
 
    epsj = lnws_fit_epsilon(x,b);
 
    nll = nll + lnws_level_nll(epsj,b,x);
 
end
 
if ~isfinite(nll)
    nll = realmax/1e100;
end
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function epsHat = lnws_fit_epsilon(x,b)
% Marginal-MLE epsilon_j for fixed b.
%
% epsilon is the PRIOR MASS AT ZERO, matching the LNWS notation.
 
obj = @(e) lnws_level_nll(e,b,x);
 
opts = optimset('Display','off','TolX',1e-8);
 
% Permit solutions arbitrarily close to either boundary without evaluating
% exactly at 0 or 1.
epsHat = fminbnd(obj,1e-8,1-1e-8,opts);
 
% Explicitly check the two limiting boundary cases.
vals = [obj(1e-10),obj(epsHat),obj(1-1e-10)];
[~,ii] = min(vals);
 
if ii==1
    epsHat = 1e-10;
elseif ii==3
    epsHat = 1-1e-10;
end
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function nll = lnws_level_nll(epsj,b,x)
% Negative log of
%
% m(x) = eps_j m0(x) + (1-eps_j) m1(x)
%
% with
%
% m0 = central chi-square_3 density,
%
% m1 = marginal under lambda ~ Exp(rate=b).
 
[logm0,logm1] = lnws_log_marginals(x,b);
 
A = log(max(epsj,realmin))       + logm0;
B = log(max(1-epsj,realmin))     + logm1;
 
mx = max(A,B);
 
logm = mx + log(exp(A-mx)+exp(B-mx));
 
nll = -sum(logm);
 
if ~isfinite(nll)
    nll = realmax/1e100;
end
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [logm0,logm1] = lnws_log_marginals(x,b)
% Stable log marginals from Eq. (8).
%
% m0(x) = sqrt[x/(2*pi)] exp(-x/2)
%
% m1(x) = b/sqrt(1+2b)
%         exp[-b*x/(1+2b)]
%         erf( sqrt[x/(2+4b)] ).
 
x = max(x,1e-300);
 
q = 1+2*b;
 
logm0 = 0.5*log(x) - 0.5*log(2*pi) - 0.5*x;
 
z = sqrt(x/(2*q));
 
E = erf(z);
E = max(E,realmin);
 
logm1 = log(b) - 0.5*log(q) - (b/q)*x + log(E);
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function delta = lnws_delta_mean(x,epsj,b)
% Posterior mean of noncentrality lambda, Eq. (9).
%
% delta(x) = [1-p_j(x)] delta_E(x),
%
% where p_j is posterior probability of lambda=0 and
%
% delta_E(x) =
%
% [1+2b+x
%   + sqrt{2x(1+2b)/pi}
%     exp{-x/(2+4b)} / erf{sqrt[x/(2+4b)]}]
% / (1+2b)^2.
 
x = max(x,1e-14);
 
q = 1+2*b;
 
[logm0,logm1] = lnws_log_marginals(x,b);
 
A = log(max(epsj,realmin))   + logm0;
B = log(max(1-epsj,realmin)) + logm1;
 
mx = max(A,B);
 
p0 = exp(A-mx)./(exp(A-mx)+exp(B-mx));
 
z = sqrt(x/(2*q));
 
E = erf(z);
E = max(E,realmin);
 
ratioTerm = sqrt(2*x*q/pi).*exp(-x/(2*q))./E;
 
deltaE = (q + x + ratioTerm)/(q^2);
 
delta = (1-p0).*deltaE;
 
end
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SURESHRINK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [dhat,lambdaHat] = sureshrink_level(d,sigma)
% Donoho--Johnstone HYBRID levelwise soft SureShrink.
% Standardize d by estimated sigma, use the sparse-level test
%   mean((d/sigma).^2)-1 <= (log2(n))^(3/2)/sqrt(n),
% and revert to the levelwise universal threshold if sparse.
% Otherwise minimize exact soft-threshold SURE for t in
% [0,sqrt(2*log(n))].  This checks 0, each |d|/sigma<=cap,
% and the universal cap, including the cap when it is not a knot.
sz=size(d); d=d(:).'; n=numel(d);
if n==0, dhat=reshape(d,sz); lambdaHat=0; return; end
assert(isfinite(sigma) && sigma>0);
z=abs(d)/sigma;
cap=sqrt(2*log(max(n,2)));
gamma=(log2(max(n,2)))^(3/2)/sqrt(n);
if mean(z.^2)-1 <= gamma
    threshold=cap;
else
    zs=sort(z(:),'ascend');
    z2=zs.^2; cum2=cumsum(z2);
    valid=find(zs<=cap); k=valid(:);
    % SURE/sigma^2 = n-2k+sum_{i<=k}z_i^2+(n-k)t^2.
    risk=[n; n-2*k+cum2(k)+(n-k).*z2(k)];
    cand=[0;zs(k)];
    kCap=numel(k);
    if kCap==0, sumBelow=0; else, sumBelow=cum2(kCap); end
    capRisk=n-2*kCap+sumBelow+(n-kCap)*cap^2;
    risk=[risk;capRisk];cand=[cand;cap];
    [~,ix]=min(risk);
    threshold=cand(ix);
end
lambdaHat=sigma*threshold;
dhat=reshape(soft_threshold(d,lambdaHat),sz);
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAUSS-LEGENDRE RULE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [x,w] = gauss_legendre_rule(n)
% n-point Gauss-Legendre rule on [-1,1], computed by Golub-Welsch.
 
k = 1:(n-1);
 
beta = k./sqrt(4*k.^2-1);
 
T = diag(beta,1) + diag(beta,-1);
 
[V,D] = eig(T);
 
x = diag(D);
 
[x,idx] = sort(x);
 
V = V(:,idx);
 
w = 2*(V(1,:).^2).';
 
x = x(:);
w = w(:);
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERIODIC ORTHONORMAL DWT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [aL,D] = my_periodic_dwt(x,h,g,L)
 
x = x(:).';
 
N = length(x);
 
J = round(log2(N));
 
D = cell(1,J);
 
a = x;
 
for j = J-1:-1:L
 
    [a,d] = analysis_step_periodic(a,h,g);
 
    D{j+1} = d;
 
end
 
aL = a;
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function x = my_periodic_idwt(aL,D,h,g,L,J)
 
a = aL(:).';
 
for j = L:J-1
 
    d = D{j+1};
 
    a = synthesis_step_periodic(a,d,h,g);
 
end
 
x = a;
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [a,d] = analysis_step_periodic(x,h,g)
 
x = x(:).';
 
n = length(x);
 
if mod(n,2) ~= 0
    error('Signal length must be even.');
end
 
nh = length(h);
 
n2 = n/2;
 
a = zeros(1,n2);
d = zeros(1,n2);
 
offset = 0:(nh-1);
 
for k = 0:(n2-1)
 
    ind = mod(2*k+offset,n)+1;
 
    xx = x(ind);
 
    a(k+1) = sum(h.*xx);
    d(k+1) = sum(g.*xx);
 
end
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function x = synthesis_step_periodic(a,d,h,g)
 
a = a(:).';
d = d(:).';
 
n2 = length(a);
 
n = 2*n2;
 
nh = length(h);
 
x = zeros(1,n);
 
offset = 0:(nh-1);
 
for k = 0:(n2-1)
 
    ind = mod(2*k+offset,n)+1;
 
    x(ind) = x(ind) ...
           + a(k+1)*h ...
           + d(k+1)*g;
 
end
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cau/Cau SECOND POSTERIOR MODE (2mod / SPoM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [me,bre] = bams_bayesrule(d,mu,tau,epsi)
% BAMS marginal distribution and posterior-mean Bayes rule.
 
de = 0.5*sqrt(2*mu)*exp(-abs(d)*sqrt(2*mu));
 
m = (tau*exp(-abs(d)/tau) - 1/sqrt(2*mu)* ...
     exp(-sqrt(2*mu)*abs(d))) / (2*tau^2 - 1/mu);
 
me = epsi*de + (1-epsi)*m;
 
if d >= 0
 
    br = (tau*(tau^2 - 1/(2*mu))*d*exp(-d/tau) + ...
          tau^2/mu*(exp(-sqrt(2*mu)*d)-exp(-d/tau))) / ...
         ((tau^2 - 1/(2*mu))* ...
          (tau*exp(-d/tau)-1/sqrt(2*mu)*exp(-sqrt(2*mu)*d)));
 
else
 
    br = (tau*(tau^2 - 1/(2*mu))*d*exp(d/tau) - ...
          tau^2/mu*(exp(sqrt(2*mu)*d)-exp(d/tau))) / ...
         ((tau^2 - 1/(2*mu))* ...
          (tau*exp(d/tau)-1/sqrt(2*mu)*exp(sqrt(2*mu)*d)));
 
end
 
bre = (1-epsi)*m*br / ((1-epsi)*m + epsi*de);
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function q = caucau_qprior(tgrid)
 
q = tgrid./(pi^2*sinh(tgrid));
iz = abs(tgrid) < 1e-10;
q(iz) = 1/pi^2;
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOUBLE-EXPONENTIAL POSTERIOR MEAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
function muhat = caucau_mean_vectorized(d,tau,sigma,tgrid)
 
original_size = size(d);
d = d(:).';
sgn  = sign(d);
dabs = abs(d);
tgrid = tgrid(:);
 
avec = tau*exp(tgrid);
q = caucau_qprior(tgrid);
 
Zplus  = (avec-dabs)/sigma;
Zminus = (avec+dabs)/sigma;
 
Wplus  = q.*exp(-0.5*Zplus.^2);
Wminus = q.*exp(-0.5*Zminus.^2);
 
den = trapz(tgrid,Wplus+Wminus,1);
num = trapz(tgrid,avec.*(Wplus-Wminus),1);
 
muabs = zeros(size(dabs));
good = den > realmin & isfinite(den) & isfinite(num);
muabs(good) = num(good)./den(good);
muabs = max(muabs,0);
 
muhat = sgn.*muabs;
muhat = reshape(muhat,original_size);
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAUSSIAN Cau/Cau POSTERIOR MEDIAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
function estimates=E3_extended_methods(d,sigma)
% E3_EXTENDED_METHODS   Twelve additional classical / adapted benchmarks.
%
% Column index in this function -> Experiment III table index:
% [7 FDR, 8 EBBlock, 10 NeighCoeff, 11 NeighShrinkSURE, 12 SURELET,
%  14 NeighBlock, 15 SureBlock, 16 BlockJS, 17 Jeffreys,
%  18 LocalWiener, 19 Context, 20 HMT].
%
% IMPORTANT IMPLEMENTATION AUDIT:
% FDR, NeighCoeff, NeighBlock and BlockJS implement established formulas.
% Others explicitly use 1D reductions or PROVISIONAL implementations.
% In particular EBBlock and HMT are research-level variants; the
% code DOES NOT assert that they reproduce the original published software.
% The resulting comparisons cannot be submitted as exact-published-method
% comparisons without a separate implementation audit.
N=numel(d); J=round(log2(N)); L=3;
estimates=repmat(d(:),1,12);
sigma2=max(sigma^2,realmin);
Tref=sigma*sqrt(2*log(N));
% FDR: Abramovich et al., BH hard threshold, two-sided Gaussian p values.
absz=abs(d(2^L+1:end))/max(sigma,eps);
p=erfc(absz/sqrt(2));
ps=sort(p); K=numel(ps); q=.05;
i=find(ps<=q*(1:K)'/K,1,'last');
if isempty(i), estimates(2^L+1:end,1)=0;
else
  cut=ps(i); xx=d(2^L+1:end);
  estimates(2^L+1:end,1)=xx.*(p<=cut);
end
 
for j=L:J-1
 ind=(2^j+1):2^(j+1); x=d(ind).'; n=numel(x);
 [letz,blockz,abez]=E3_repaired_three(x,sigma);
 % 8: empirical-Bayes GAUSSIAN block-mixture prototype (not exact ABS02)
 estimates(ind,2)=ebblock_provisional(x,sigma2).';
 % 10: Cai--Silverman NeighCoeff; lambda*3=2 log N
 local=x.^2+circshift(x.^2,[0,1])+circshift(x.^2,[0,-1]);
 estimates(ind,3)=(x.*max(1-2*log(N)*sigma2./max(local,realmin),0)).';
 % 11: NeighShrink-SURE 1D periodic adaptation with levelwise SURE.
 estimates(ind,4)=neighshrink_sure_1d(x,sigma,Tref).';
 % 12: 1D orthonormal SURE-LET with the published 12*sigma^2 Gaussian
 % kernel; this is NOT the original 2D undecimated/interscale algorithm.
 estimates(ind,5)=letz.';
 % 14: Cai--Silverman NeighBlock, overlapping windows.
 estimates(ind,6)=neighblock_cs(x,sigma2,N).';
 % 15: Cai--Zhou 2009 hybrid SureBlock with constrained SURE search.
 estimates(ind,7)=blockz.';
 % 16: Cai 1999 BlockJS
 estimates(ind,8)=block_js(x,sigma2,N).';
 % 17: exact Figueiredo--Nowak 2001 ABE shrinkage (equation 17):
 % (x^2 - 3*sigma^2)_+/x, with value zero at x=0.
 estimates(ind,9)=abez.';
 % 18: Mihcak-style locally adaptive variance Wiener, 1D reduction.
 estimates(ind,10)=local_wiener(x,sigma2).';
 % 19: CYV context is computed using the full wavelet vector below,
 % so that each level has access to its parent coefficients.
end
% 19: Chang--Yu--Vetterli (2000) context grouping, 1D orthonormal variant.
estimates(:,11)=cyv_context_coefficients(d,sigma,J,L);
% 20: a two-state Gaussian hidden-Markov TREE model trained by EM in 1D.
% This is an experimental implementation of the HMT class, NOT an assertion
% of identity to a published Crouse/Nowak/Baraniuk codebase.
estimates(:,12)=hmt_1d_EM(d,sigma2,L,J);
assert(all(isfinite(estimates(:))),'Nonfinite extended benchmark output');
end
 
function z=ebblock_provisional(x,sigma2)
n=numel(x); blockSize=4; z=x;
for i=1:blockSize:n
 idx=i:min(i+blockSize-1,n); a=x(idx); m=numel(a);
 % Empirical Bayes spike and Gaussian slab in each block.
 tau2=max(mean(x.^2)-sigma2, 1e-5*sigma2);
 b0=exp(-sum(a.^2)/(2*sigma2)) / (sigma2^(m/2));
 b1=exp(-sum(a.^2)/(2*(sigma2+tau2))) / ((sigma2+tau2)^(m/2));
 pi1=min(.95,max(.05,mean(x.^2>2*sigma2)));
 post=(pi1*b1)/(pi1*b1+(1-pi1)*b0+realmin);
 z(idx)=a.*post.*tau2/(sigma2+tau2);
end
end
 
function z=neighshrink_sure_1d(x,sigma,Tref)
sigma2=sigma^2; best=Inf; z=zeros(size(x));
for r=1:3
  m=x.^2;
  for lag=1:r, m=m+circshift(x.^2,[0 lag])+circshift(x.^2,[0 -lag]); end
  for c=[.5 .75 1 1.25]
    T=(c*Tref)^2;
    w=max(1-T./max(m,realmin),0);
    t=x.*w;
    div=zeros(size(x)); active=m>T;
    div(active)=w(active)+2*T*x(active).^2./m(active).^2;
    sure=sum((t-x).^2+2*sigma2*div-sigma2);
    if sure<best, best=sure; z=t; end
  end
end
end
 
 
function z=neighblock_cs(x,sigma2,N)
n=numel(x); L0=max(1,floor(log(N)/2));
L1=max(1,floor(L0/2)); L=L0+2*L1; z=zeros(size(x));
for i=1:L0:n
  e=min(n,i+L0-1);
  ii=mod(((i-L1):(i+L0+L1-1))-1,n)+1; % fixed-size periodic window
  SS=sum(x(ii).^2);
  w=max(1-4.50524*L*sigma2/max(SS,realmin),0);
  z(i:e)=w*x(i:e);
end
end
 
function z=block_js(x,sigma2,N)
% Cai 1999 uses L=[log(n)]=floor(natural log n); at n=1024, L=6.
n=numel(x); L=max(1,floor(log(N))); z=zeros(size(x));
for i=1:L:n
  idx=i:min(n,i+L-1); m=numel(idx); v=sum(x(idx).^2);
  z(idx)=x(idx)*max(1-4.50524*m*sigma2/max(v,realmin),0);
end
end
 
 
function z=local_wiener(x,sigma2)
v=max(local_average(x.^2,7)-sigma2,0);
z=x.*(v./max(v+sigma2,realmin));
end
 
function a=local_average(x,n)
a=zeros(size(x));h=(n-1)/2;
for lag=-h:h, a=a+circshift(x,[0 lag]); end
a=a/n;
end
 
function z=hmt_1d_EM(d,sigma2,L,J)
% Two-state zero-mean Gaussian HMT, level-dependent observations.
% Tree has 2^L roots (eight), each node has 2 descendants in next level.
% An EM-like update uses normalized upward/downward messages and global
% shared variances/transitions. Denominator and weight floors safeguard EM.
% Train from noisy data ONLY. This is a 1D adaptation requiring audit.
x=cell(J-L,1);
for j=L:J-1, x{j-L+1}=d(2^j+1:2^(j+1)).'; end
c=J-L; tau2=[0.05, max(mean(d(2^L+1:end).^2)-sigma2,0.2*sigma2)*4];
pi0=[.85;.15]; A=[.94 .06; .14 .86];
for iter=1:15
  Up=cell(c,1); Em=cell(c,1); ratios=cell(c,1);
  for jj=c:-1:1
    obs=x{jj}; v=sigma2+tau2;
    logB=[-obs.^2/(2*v(1))-.5*log(v(1)); ...
          -obs.^2/(2*v(2))-.5*log(v(2))];
    logB=bsxfun(@minus,logB,max(logB,[],1));
    B=exp(logB);
    if jj<c
      kid=Up{jj+1};
      msg=A*kid;
      B=B.*msg(:,1:2:end).*msg(:,2:2:end);
    end
    B=B./max(sum(B,1),realmin);
    Up{jj}=max(B,realmin); Em{jj}=obs;
  end
  Gam=cell(c,1);
  G=pi0.*Up{1}; G=G./max(sum(G,1),realmin);
  Gam{1}=G; counts=zeros(2,2);
  for jj=2:c
    B=Up{jj}; parents=ceil((1:size(B,2))/2);
    oldG=Gam{jj-1}(:,parents);
    denom=A*B;
    xi=zeros(2,size(B,2));
    for child=1:2
      v=oldG(child,:).*A(child,:)';
      v=v.*B./max(denom(child,:),realmin);
      counts(child,:)=counts(child,:)+sum(v,2)';
      xi=xi+v;
    end
    xi=xi./max(sum(xi,1),realmin);
    Gam{jj}=xi;
  end
  pi0=max(mean(Gam{1},2),1e-5);pi0=pi0/sum(pi0);
  A=max(counts,1e-6);A=A./sum(A,2);
  for s=1:2
    weighted=0; mass=0;
    for jj=1:c
      weighted=weighted+sum(Gam{jj}(s,:).*(x{jj}.^2));
      mass=mass+sum(Gam{jj}(s,:));
    end
    tau2(s)=max(weighted/max(mass,1e-6)-sigma2,.01*sigma2);
  end
  tau2(2)=max(tau2(2),tau2(1)+.05*sigma2);
end
z=d;
for jj=1:c
 j=L+jj-1; ind=(2^j+1):2^(j+1);
 shrink=Gam{jj}(1,:)*tau2(1)/(sigma2+tau2(1))+ ...
        Gam{jj}(2,:)*tau2(2)/(sigma2+tau2(2));
 z(ind)=x{jj}.*shrink;
end
end
 
 
function [letOut,blockOut,abeOut]=E3_repaired_three(x,sigma)
% Repaired Experiment III literature competitors, using the identical
% coefficient input and plug-in noise estimate as the main benchmark.
% 1. 1D TWO-basis orthonormal SURE-LET, exponent 12 sigma^2;
% 2. Cai & Zhou (2009) complete hybrid SureBlock with exact SURE search;
% 3. Figueiredo & Nowak (2001) ABE with factor 3 sigma^2.
% The SURE-LET here is explicitly a one-dimensional reduction, not a
% reproduction of the published 2D/interscale image algorithm.
sz=size(x);x=x(:).';sigma2=max(sigma^2,realmin);
letOut=reshape(sure_let_1d(x,sigma),sz);
blockOut=reshape(sureblock_1d(x,sigma2),sz);
abeOut=reshape(x.*max(1-3*sigma2./max(x.^2,realmin),0),sz);
end
 
function z=sure_let_1d(x,sigma)
% 1D orthonormal LET, a TWO-basis reduction of Luisier--Blu--Unser.
% h(x) = a1*x + a2*x*exp(-x^2/(12*sigma^2)).
% Its optimal unconstrained (a1,a2) minimize plug-in SURE.
% This is NOT the original 2D / undecimated / interscale SURE-LET.
sigma2=max(sigma^2,realmin);
t=x(:)/sqrt(sigma2); e=exp(-t.^2/12);
A=[x(:),x(:).*e];
V=[ones(numel(x),1),e.*(1-t.^2/6)];
G=A'*A; rhs=A'*x(:)-sigma2*sum(V,1)';
% Pseudoinverse only for nearly dependent columns; no arbitrary ridge.
if rcond(G)<1e-11
  a=pinv(G)*rhs;
else
  a=G\rhs;
end
z=reshape(A*a,size(x));
end
 
function z=sureblock_1d(x,sigma2)
% Cai--Zhou (2009) hybrid SureBlock, complete search over integer L.
% T_d=mean((x/sigma)^2-1), gamma_d=log2(d)^(3/2)/sqrt(d).
% For T_d<=gamma_d use garrote lambda=2*log(d); otherwise minimize
% exact blockwise SURE over integer 1<=L<=floor(sqrt(d)),
% max(L-2,0)<=lambda<=2*L*log(d) (lambda in sigma^2 units).
% The last block is shortened if L does not divide d. Its actual
% dimension is used in the block SURE and divergence formula.
% For fixed L, SURE is piecewise quadratic; check each breakpoint
% and the stationary point in each interval. No coarse lambda grid.
sz=size(x);x=x(:).';d=numel(x);
sigma=max(sqrt(sigma2),sqrt(realmin));u=x/sigma;
if d<2,z=reshape(x,sz);return;end
Td=mean(u.^2)-1; gamma=(log2(d)^1.5)/sqrt(d);
if Td<=gamma
 z=reshape(x.*max(1-2*log(d)./max(u.^2,realmin),0),sz);
 return;
end
best=Inf;bestL=1;bestLambda=0;
for L=1:floor(sqrt(d))
 nBlocks=ceil(d/L);energy=zeros(nBlocks,1);m=zeros(nBlocks,1);
 for ib=1:nBlocks
   idx=(ib-1)*L+1:min(d,ib*L);
   energy(ib)=sum(u(idx).^2);m(ib)=numel(idx);
 end
 [E,order]=sort(energy); m=m(order);
 lo=max(L-2,0);hi=2*L*log(d);
 preE=[0;cumsum(E)];preM=[0;cumsum(m)];
 invE=1./max(E,realmin);
 tailInv=[flipud(cumsum(flipud(invE)));0];
 tailBeta=[flipud(cumsum(flipud((m-2).*invE)));0];
 [eu,kidx]=unique(E,'last');good=(eu>=lo & eu<=hi);
 cand=[lo;eu(good);hi];
 inactive=[sum(E<=lo);kidx(good);sum(E<=hi)];
 % stationary points in intervals with exactly k inactive blocks;
 % the convex quadratic is minimized at tailBeta/tailInv.
 intervalLow=max(lo,[0;E]);
 intervalHigh=min(hi,[E;Inf]);
 k=(0:nBlocks)';
 valid=(intervalHigh>intervalLow) & (tailInv>0);
 stationary=zeros(nBlocks+1,1);
 stationary(valid)=min(max(tailBeta(valid)./tailInv(valid), ...
          intervalLow(valid)),intervalHigh(valid));
 % Do not use an upper endpoint at a breakpoint with the wrong active
 % set; that breakpoint is included separately with all ties inactive.
 valid=valid & (stationary<intervalHigh | intervalHigh>=hi);
 cand=[cand;stationary(valid)];inactive=[inactive;k(valid)];
 risk=d+preE(inactive+1)-2*preM(inactive+1) + ...
      cand.^2.*tailInv(inactive+1)-2*cand.*tailBeta(inactive+1);
 [value,ind]=min(risk);
 if value<best,best=value;bestL=L;bestLambda=cand(ind);end
end
z=zeros(1,d);
for ib=1:ceil(d/bestL)
 idx=(ib-1)*bestL+1:min(ib*bestL,d);
 E=sum(u(idx).^2);
 z(idx)=x(idx)*max(1-bestLambda/max(E,realmin),0);
end
z=reshape(z,sz);
end