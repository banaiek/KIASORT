function bk = kiaSort_backup_results(outputPath, tag, files)
%KIASORT_BACKUP_RESULTS  Copy result files before a pass rewrites them.
%
%   bk = kiaSort_backup_results(outputPath, tag, files)
%
%   Copies each existing entry of `files` (full paths) into
%   <outputPath>/Backup/<tag>/ under its own name. Returns
%       bk.ok    true when every existing file was copied
%       bk.dir   the backup directory
%       bk.pairs n-by-2 cellstr, {original, copy}
%   so a caller can put the originals back with kiaSort_restore_results(bk).
%
%   The post-hoc passes rewrite unifiedLabels.h5 and sorted_samples.mat as a
%   pair. They are separate files, so a failure between the two writes leaves
%   spike labels that disagree with the unit table. Backing up first makes
%   that recoverable, and rolling back in the caller's catch makes the pair
%   effectively atomic.

bk = struct('ok', false, 'dir', '', 'pairs', {cell(0,2)});
if nargin < 3 || isempty(files), return; end
if ischar(files) || isstring(files), files = {char(files)}; end

bk.dir = fullfile(char(outputPath), 'Backup', char(tag));
if ~exist(bk.dir, 'dir')
    [made, msg] = mkdir(bk.dir);
    if ~made
        warning('kiaSort:backupFailed', 'Could not create %s (%s)', bk.dir, msg);
        return;
    end
end

pairs = cell(0,2);
for i = 1:numel(files)
    src = char(files{i});
    if ~exist(src, 'file'), continue; end
    [~, nm, ext] = fileparts(src);
    dst = fullfile(bk.dir, [nm ext]);
    [okCopy, msg] = copyfile(src, dst, 'f');
    if ~okCopy
        warning('kiaSort:backupFailed', 'Could not back up %s (%s)', src, msg);
        return;
    end
    pairs(end+1, :) = {src, dst}; %#ok<AGROW>
end

bk.pairs = pairs;
bk.ok    = true;
end
