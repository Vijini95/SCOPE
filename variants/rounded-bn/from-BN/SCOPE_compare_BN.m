function SCOPE_compare_BN(originalResultsDir,roundedResultsDir)
% SCOPE_compare_BN  Compare original versus rounded-B_N SCOPE AMSE results.
%
% USAGE:
%   SCOPE_compare_BN
%   SCOPE_compare_BN('C:\Brani2\WorkTAMU\SCOPE\SCOPE_Exp123_results')
%   SCOPE_compare_BN('C:\path\to\old_results','C:\path\to\new_results')
%
% With no arguments, the rounded results are expected beside this .m file.
% The original results are sought in this folder and up to three parent
% folders (use explicit directory arguments if they are elsewhere).
%
% This compares existing summary CSVs only; it does not run simulations.
% Required source files: Experiment1/2/3_AMSE_all.csv in corresponding
% Exp<e>_N<N>_M<M>[_SNR...] run subdirectories. MATLAB R2018b or later.

here = fileparts(mfilename('fullpath'));

if nargin < 2 || isempty(roundedResultsDir)
    roundedResultsDir = fullfile(here,'SCOPE_Exp123_results_BNrounded');
end
roundedResultsDir = char(roundedResultsDir);
assert(isfolder(roundedResultsDir), ...
    ['Cannot find rounded-B_N results: %s\n' ...
     'Pass the original and rounded results folder paths as arguments.'], ...
     roundedResultsDir);

if nargin < 1 || isempty(originalResultsDir)
    originalResultsDir = '';
    parent = here;
    for level = 0:3
        tryDir = fullfile(parent,'SCOPE_Exp123_results');
        if isfolder(tryDir)
            originalResultsDir = tryDir;
            break;
        end
        nextParent = fileparts(parent);
        if strcmp(nextParent,parent), break; end
        parent = nextParent;
    end
end
originalResultsDir = char(originalResultsDir);
assert(isfolder(originalResultsDir), ...
    ['Cannot find original results automatically.\n' ...
     'Use: SCOPE_compare_BN(''C:\path\to\SCOPE_Exp123_results'')\n' ...
     'Or supply both folder paths as two arguments.']);

fprintf('Original results: %s\n',originalResultsDir);
fprintf('Rounded-B_N results: %s\n',roundedResultsDir);

runDirs = dir(fullfile(roundedResultsDir,'Exp*_N*_M*'));
combined = table();  % IMPORTANT: do not name this variable "all"!
matchedFolders = 0;

for ii = 1:numel(runDirs)
    if ~runDirs(ii).isdir, continue; end
    folderName = runDirs(ii).name;
    parts = regexp(folderName,'^Exp([123])_N(\d+)_M(\d+)(.*)$','tokens','once');
    if isempty(parts), continue; end

    expNo = str2double(parts{1});
    originalFile = fullfile(originalResultsDir,folderName, ...
                        sprintf('Experiment%d_AMSE_all.csv',expNo));
    roundedFile = fullfile(roundedResultsDir,folderName, ...
                       sprintf('Experiment%d_AMSE_all.csv',expNo));
    if ~isfile(originalFile) || ~isfile(roundedFile)
        fprintf('Skipping %s: matching original/new AMSE CSV missing.\n',folderName);
        continue;
    end

    original = readtable(originalFile);
    rounded = readtable(roundedFile);

    % Experiment III identifies a procedure in "Method"; Experiments I and
    % II identify it in "Procedure". Normalize to one join key.
    if expNo == 3
        nameField = 'Method';
    else
        nameField = 'Procedure';
    end
    mustHave = {'N','Signal','SNR',nameField,'AMSE'};
    assert(all(ismember(mustHave,original.Properties.VariableNames)), ...
        'Missing AMSE columns in original CSV: %s',originalFile);
    assert(all(ismember(mustHave,rounded.Properties.VariableNames)), ...
        'Missing AMSE columns in rounded CSV: %s',roundedFile);

    original.Procedure = string(original.(nameField));
    rounded.Procedure = string(rounded.(nameField));
    original.Signal = string(original.Signal);
    rounded.Signal = string(rounded.Signal);
    keys = {'N','Signal','SNR','Procedure'};

    original = original(:,[keys,{'AMSE'}]);
    rounded = rounded(:,[keys,{'AMSE'}]);
    original.Properties.VariableNames{end} = 'AMSE_original';
    rounded.Properties.VariableNames{end} = 'AMSE_roundedBN';

    % These keys are unique in the Experiment I/II/III summary tables.
    % Defensive check prevents accidental Cartesian-product joins.
    assert(height(unique(original(:,keys),'rows')) == height(original), ...
        'Duplicate original summary keys in %s',originalFile);
    assert(height(unique(rounded(:,keys),'rows')) == height(rounded), ...
        'Duplicate rounded summary keys in %s',roundedFile);

    comparison = innerjoin(original,rounded,'Keys',keys);
    if height(comparison) ~= height(original) || ...
       height(comparison) ~= height(rounded)
        warning('Different method/condition sets in %s; only matching rows compared.', ...
                folderName);
    end
    if isempty(comparison)
        fprintf('Skipping %s: no matching rows.\n',folderName);
        continue;
    end

    comparison.Delta_AMSE = comparison.AMSE_roundedBN - comparison.AMSE_original;
    comparison.PctChange_AMSE = 100*comparison.Delta_AMSE ./ ...
                                max(comparison.AMSE_original,eps);
    comparison.Experiment = repmat(expNo,height(comparison),1);
    comparison.RunFolder = repmat({folderName},height(comparison),1);
    comparison = comparison(:,{'Experiment','N','Signal','SNR','Procedure', ...
        'AMSE_original','AMSE_roundedBN','Delta_AMSE', ...
        'PctChange_AMSE','RunFolder'});

    if isempty(combined)
        combined = comparison;
    else
        combined = [combined;comparison]; %#ok<AGROW>
    end
    matchedFolders = matchedFolders + 1;
    fprintf('%s: %d matched rows; median |AMSE change| %.4f%%; max %.4f%%.\n', ...
       folderName,height(comparison),median(abs(comparison.PctChange_AMSE)), ...
       max(abs(comparison.PctChange_AMSE)));
end

if isempty(combined)
    fprintf(['No matching completed runs found. Check that both results folders ' ...
             'contain identically named Exp*_N*_M* run subfolders.\n']);
    return;
end

combined = sortrows(combined,{'Experiment','N','Signal','SNR','Procedure'});
outFile = fullfile(roundedResultsDir,'BNrounded_versus_original_AMSE.csv');
writetable(combined,outFile);
fprintf('\nCompared %d rows in %d run folders.\n',height(combined),matchedFolders);
fprintf('Comparison saved to:\n%s\n',outFile);
fprintf(['PctChange_AMSE = 100*(rounded - original)/original; ' ...
         'negative values favor rounded B_N.\n']);
fprintf('These are differences between means, not paired Monte Carlo confidence intervals.\n');
end
