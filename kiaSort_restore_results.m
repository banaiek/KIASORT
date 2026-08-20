function ok = kiaSort_restore_results(bk)
%KIASORT_RESTORE_RESULTS  Put back the files kiaSort_backup_results copied.
%
%   ok = kiaSort_restore_results(bk)
%
%   Used in the catch of a two-file write so a half-finished pass leaves the
%   tree as it was rather than with spike labels that disagree with the unit
%   table. Returns false if any file could not be restored, in which case the
%   copies are still in bk.dir.

ok = false;
if ~isstruct(bk) || ~isfield(bk, 'pairs') || isempty(bk.pairs), return; end

ok = true;
for i = 1:size(bk.pairs, 1)
    src = bk.pairs{i, 1};
    dst = bk.pairs{i, 2};
    if ~exist(dst, 'file'), ok = false; continue; end
    [okCopy, msg] = copyfile(dst, src, 'f');
    if ~okCopy
        ok = false;
        warning('kiaSort:restoreFailed', 'Could not restore %s (%s)', src, msg);
    end
end
end
