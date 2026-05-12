function [reports, halted, haltPass] = kiaSort_drift_merge_posthoc_iterative(outputPath, varargin)
%KIASORT_DRIFT_MERGE_POSTHOC_ITERATIVE Progressive drift merge in up to 4 passes.
%
% Each pass relaxes the waveform-correlation threshold while tightening
% the other gates (density anti-corr, amp ratio, ACG similarity,
% coverage gain). During iteration the canonical unifiedLabels.h5 is
% overwritten with the cumulative merged state so pass N operates on
% the output of passes 1..N-1.
%
% A pass is COMMITTED only if it produces strictly more than
% `minMergesToContinue` unit merges (default 10). The first pass that
% does not clear that threshold halts the iteration and its proposed
% merges are discarded.
%
% Output policy (single final output, no per-pass files):
%   'overwrite' = false (default)
%       The canonical unifiedLabels.h5 is left untouched at the end.
%       The FINAL cumulative merged labels are written to
%       unifiedLabels_merged.h5.
%   'overwrite' = true
%       The final cumulative merged labels REPLACE unifiedLabels.h5.
%       The pre-iteration state is backed up to
%       unifiedLabels_predrift.h5 (created only if no such backup
%       already exists).
%
% The final driftReport is saved to drift_merge_report.mat regardless
% (it reflects only the LAST executed pass, whose diag/evidence arrays
% are the most informative for debugging; the full per-pass reports
% are returned in the `reports` output for programmatic inspection).
%
% NAME/VALUE PAIRS
%   'minMergesToContinue'  (default 10)    strict-greater-than gate
%   'isiThresholdPct'      (default 1.0)   constant across passes
%   'overwrite'            (default false) commit final merge into
%                                          unifiedLabels.h5
%   'passOpts'             (default [])    override the 4-pass schedule
%                                          with a cell array of structs
%                                          containing fields:
%                                            corrThreshold
%                                            densityAntiCorr
%                                            ampRatioMax
%                                            acgCorrThreshold
%                                            coverageGainMin
%   'mainArgs'             (default {})    extra name/value args
%                                          forwarded to
%                                          kiaSort_drift_merge_posthoc
%                                          on every pass (e.g. cfg,
%                                          maxDriftUm, debugFigs)
%   'verbose'              (default true)
%
% OUTPUTS
%   reports    cell array of driftReport structs for each COMMITTED
%              pass (pass that exceeded minMergesToContinue).
%   halted     true if iteration stopped early due to a low merge count.
%   haltPass   index of the pass that halted (NaN if all passes ran).

p = inputParser;
p.addRequired('outputPath', @(x) ischar(x) || isstring(x));
p.addParameter('minMergesToContinue', 10, @isscalar);
p.addParameter('isiThresholdPct',     1.0, @isscalar);
p.addParameter('overwrite',           false, @islogical);
p.addParameter('passOpts',            [], @(x) isempty(x) || iscell(x));
p.addParameter('mainArgs',            {}, @iscell);
p.addParameter('verbose',             true, @islogical);
p.parse(outputPath, varargin{:});
opt = p.Results;

if isempty(opt.passOpts)
    opt.passOpts = {
        struct('corrThreshold', 0.90, 'densityAntiCorr', -0.20, ...
               'ampRatioMax',   2.00, 'acgCorrThreshold', 0.20, ...
               'coverageGainMin', 0.05);
        struct('corrThreshold', 0.85, 'densityAntiCorr', -0.50, ...
               'ampRatioMax',   1.75, 'acgCorrThreshold', 0.25, ...
               'coverageGainMin', 0.10);
        struct('corrThreshold', 0.80, 'densityAntiCorr', -0.65, ...
               'ampRatioMax',   1.50, 'acgCorrThreshold', 0.30, ...
               'coverageGainMin', 0.15);
        struct('corrThreshold', 0.75, 'densityAntiCorr', -0.75, ...
               'ampRatioMax',   1.50, 'acgCorrThreshold', 0.40, ...
               'coverageGainMin', 0.25);
    };
end

outPath      = char(outputPath);
resSortedDir = fullfile(outPath, 'RES_Sorted');
uLabelsPath  = fullfile(resSortedDir, 'unifiedLabels.h5');
uMergedPath  = fullfile(resSortedDir, 'unifiedLabels_merged.h5');

if ~exist(uLabelsPath, 'file')
    error('kiaSort_drift_merge_posthoc_iterative:missingLabels', ...
        'unifiedLabels.h5 not found at %s', uLabelsPath);
end

% In-memory-ish pristine backup: copy the canonical labels aside so the
% wrapper can restore them at the end when 'overwrite' is false (or use
% them as the content of unifiedLabels_predrift.h5 when 'overwrite' is
% true). The file is temporary and deleted before returning.
pristineTmp = fullfile(resSortedDir, 'unifiedLabels_iter_pristine_TMP.h5');
if exist(pristineTmp, 'file'), delete(pristineTmp); end
copyfile(uLabelsPath, pristineTmp);

reports  = {};
halted   = false;
haltPass = NaN;
anyCommit = false;

try
    for ip = 1:numel(opt.passOpts)
        po = opt.passOpts{ip};

        if opt.verbose
            fprintf('\n======================================================\n');
            fprintf(' Iterative drift-merge: pass %d of %d\n', ip, numel(opt.passOpts));
            fprintf('   corrThreshold    = %.2f\n', po.corrThreshold);
            fprintf('   densityAntiCorr  = %.2f\n', po.densityAntiCorr);
            fprintf('   ampRatioMax      = %.2f\n', po.ampRatioMax);
            fprintf('   acgCorrThreshold = %.2f\n', po.acgCorrThreshold);
            fprintf('   coverageGainMin  = %.2f\n', po.coverageGainMin);
            fprintf('   isiThresholdPct  = %.2f  (constant)\n', opt.isiThresholdPct);
            fprintf('======================================================\n');
        end

        % Run the pass without touching unifiedLabels.h5 from inside the
        % main function. The wrapper will decide whether to commit by
        % promoting unifiedLabels_merged.h5 afterwards.
        report = kiaSort_drift_merge_posthoc(outPath, ...
            'corrThreshold',     po.corrThreshold, ...
            'densityAntiCorr',   po.densityAntiCorr, ...
            'ampRatioMax',       po.ampRatioMax, ...
            'acgCorrThreshold',  po.acgCorrThreshold, ...
            'coverageGainMin',   po.coverageGainMin, ...
            'isiThresholdPct',   opt.isiThresholdPct, ...
            'overwriteOriginal', false, ...
            'dryRun',            false, ...
            'verbose',           opt.verbose, ...
            opt.mainArgs{:});

        nMerges = report.nUnitsBefore - report.nUnitsAfter;
        if opt.verbose
            fprintf('\nPass %d: %d unit merges, %d groups.\n', ...
                    ip, nMerges, numel(report.groups));
        end

        if nMerges <= opt.minMergesToContinue
            if opt.verbose
                fprintf(['Pass %d merge count (%d) <= minMergesToContinue (%d); ' ...
                         'NOT committing, halting iteration.\n'], ...
                        ip, nMerges, opt.minMergesToContinue);
            end
            halted   = true;
            haltPass = ip;
            break;
        end

        % Commit: promote this pass's merged output to the canonical
        % labels so the NEXT pass reads the cumulative state. This is
        % in-place; the wrapper undoes it at the end when overwrite is
        % false by restoring from pristineTmp.
        if exist(uLabelsPath, 'file'), delete(uLabelsPath); end
        copyfile(uMergedPath, uLabelsPath);
        anyCommit = true;
        reports{end+1} = report; %#ok<AGROW>

        if opt.verbose
            fprintf('Pass %d committed to unifiedLabels.h5 (in-place, may be restored at end).\n', ip);
        end
    end

    % ---- Final handling -------------------------------------------------
    % Output policy:
    %   * Final merged state lives in unifiedLabels.h5 (same file name
    %     as the input).
    %   * The pre-iteration state is preserved in
    %     unifiedLabels_predrift.h5 (created only when at least one
    %     pass committed and overwrite is true).
    %   * No auxiliary unifiedLabels_merged.h5 and no .mat report are
    %     left behind -- they are deleted at the end if present.
    if ~anyCommit
        % No pass ever committed. Nothing to do for unifiedLabels.h5
        % (it was never modified). Leave things untouched.
        if opt.verbose
            fprintf('\nNo pass committed; unifiedLabels.h5 is unchanged.\n');
        end
    else
        if opt.overwrite
            % Final cumulative state is already in unifiedLabels.h5
            % from the in-place commits during iteration. Just record
            % the pristine pre-iteration backup if it doesn't exist.
            preDriftBackup = fullfile(resSortedDir, 'unifiedLabels_predrift.h5');
            if ~exist(preDriftBackup, 'file')
                copyfile(pristineTmp, preDriftBackup);
                if opt.verbose
                    fprintf('Backed up pre-iteration labels -> %s\n', preDriftBackup);
                end
            end
            if opt.verbose
                fprintf('Overwrote unifiedLabels.h5 with final merged labels.\n');
            end
        else
            % Restore canonical labels to pristine -- nothing was
            % committed permanently.
            delete(uLabelsPath);
            copyfile(pristineTmp, uLabelsPath);
            if opt.verbose
                fprintf('Restored unifiedLabels.h5 to pre-iteration state (no commit).\n');
            end
        end
    end

    % ---- Cleanup -------------------------------------------------------
    % Remove the auxiliary files the iterative wrapper does not want to
    % leave behind: the per-pass merged H5 and the per-pass driftReport
    % .mat (both written unconditionally by the main function on every
    % call).
    if exist(uMergedPath, 'file'), delete(uMergedPath); end
    matReportPath = fullfile(resSortedDir, 'drift_merge_report.mat');
    if exist(matReportPath, 'file'), delete(matReportPath); end

catch ME
    % On any failure, restore canonical labels so we don't leave the
    % output directory in a half-merged state.
    if exist(pristineTmp, 'file')
        if exist(uLabelsPath, 'file'), delete(uLabelsPath); end
        copyfile(pristineTmp, uLabelsPath);
        delete(pristineTmp);
    end
    rethrow(ME);
end

% Clean up temp pristine backup.
if exist(pristineTmp, 'file'), delete(pristineTmp); end

if opt.verbose
    if halted
        fprintf('\nIteration halted at pass %d.\n', haltPass);
    elseif ~isempty(reports)
        fprintf('\nAll %d passes committed successfully.\n', numel(reports));
    end
    totalAbsorbed = 0;
    for i = 1:numel(reports)
        totalAbsorbed = totalAbsorbed + ...
            (reports{i}.nUnitsBefore - reports{i}.nUnitsAfter);
    end
    fprintf('Total units absorbed across committed passes: %d\n', totalAbsorbed);
end
end
