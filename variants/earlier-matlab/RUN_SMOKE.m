% Run from the folder containing SCOPE_Exp123_FINAL.m.
SCOPE_Exp123_FINAL('smoke','all');
base=fullfile(fileparts(which('SCOPE_Exp123_FINAL')),'SCOPE_Exp123_results');
assert(exist(fullfile(base,'Exp1_N512_SMOKE','plots','E1_AMSE_Doppler_SNR5.png'),'file')==2);
assert(exist(fullfile(base,'Exp2_N512_SMOKE','plots','E2_AMSE_Doppler_SNR5.png'),'file')==2);
assert(exist(fullfile(base,'Exp3_N512_SMOKE','plots','E3_AMSE_Doppler_SNR5.png'),'file')==2);
fprintf('All three smoke experiments completed and figures created.\n');
