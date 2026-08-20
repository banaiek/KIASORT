function [W, usedRows] = kiaSort_read_waveforms(src, rows, spk, chn, capN, drawMode)
%KIASORT_READ_WAVEFORMS  Main-channel waveforms for a set of spike rows.
%
%   [W, usedRows] = kiaSort_read_waveforms(src, rows, spk, chn, capN, drawMode)
%
%   src       from kiaSort_waveform_source
%   rows      global row indices into the RES_Sorted H5 outputs
%   spk, chn  full-length spike_idx / channelNum vectors
%   capN      max waveforms to return ([] or Inf for all; the raw path is
%             always capped because it costs one filtered read per spike)
%   drawMode  'random' (default) or 'spread' -- how rows are subsampled
%
%   W is (n x T) on each spike's own channel; usedRows are the rows that
%   produced them, in the same order. Rows whose window falls off the end of
%   the recording, or whose channel is unmapped, are dropped.
%
%   'random' uses a private RandStream so the draw is reproducible and the
%   global RNG is left alone.

W = []; usedRows = [];
if ~isstruct(src) || ~src.ok || isempty(rows), return; end
if nargin < 5, capN = []; end
if nargin < 6 || isempty(drawMode), drawMode = 'random'; end

rows = rows(:);
if ~isempty(capN) && isfinite(capN) && numel(rows) > capN
    if strcmpi(drawMode, 'spread')
        pick = round(linspace(1, numel(rows), capN))';
    else
        rs   = RandStream('threefry', 'Seed', 20240811);
        pick = sort(randperm(rs, numel(rows), capN))';
    end
    rows = rows(pick);
end

switch src.mode
    case 'h5'
        rows = sort(rows);
        W = zeros(numel(rows), src.nT);
        runs = localRuns(rows);
        at = 0;
        for r = 1:size(runs,1)
            s = runs(r,1); c = runs(r,2);
            if src.is3d
                blk = h5read(src.file, src.dset, [s 1 1], [c src.sz(2) src.sz(3)]);
                mid = ceil(size(blk,2)/2);
                blk = reshape(blk(:,mid,:), c, []);
            else
                blk = h5read(src.file, src.dset, [s 1], [c src.sz(2)]);
            end
            W(at+1:at+c, :) = double(blk);
            at = at + c;
        end
        keep = any(W, 2);
        W = W(keep, :);
        usedRows = rows(keep);

    case 'raw'
        half = src.half;
        pad  = max(4*half, 256);
        W    = nan(numel(rows), 2*half+1);
        for i = 1:numel(rows)
            ri = rows(i);
            c  = chn(ri);
            if ~isfinite(c) || c < 1 || c > numel(src.chanMap), continue; end
            row = src.chanMap(c);
            if ~isfinite(row) || row < 1, continue; end
            s = spk(ri);
            a = s - half - pad;
            b = s + half + pad;
            if a < 1 || b > src.nSamp, continue; end
            seg = localBandpass(double(src.map.Data.data(row, a:b)), src.cfg);
            c0  = half + pad + 1;
            W(i,:) = seg(c0-half : c0+half);
        end
        keep = all(isfinite(W), 2);
        W = W(keep, :);
        usedRows = rows(keep);
end
end


function y = localBandpass(x, cfg)
% Same design as kiaSort_filter_signal, with the coefficients cached and
% neither the pool nor the GPU touched -- this runs once per spike.
persistent b a key
k = sprintf('%g_%g_%g', cfg.bandpass(1), cfg.bandpass(2), cfg.samplingFrequency);
if isempty(key) || ~strcmp(key, k)
    fsl = cfg.samplingFrequency;
    up  = min(1.05 * cfg.bandpass(2), 0.95 * (fsl/2));
    [b, a] = butter(4, [cfg.bandpass(1), up] / (fsl/2), 'bandpass');
    key = k;
end
y = filtfilt(b, a, x(:))';
end


function runs = localRuns(idx)
% Contiguous [start count] runs so scattered rows are read in few calls.
idx = idx(:);
brk = [1; find(diff(idx) ~= 1) + 1; numel(idx)+1];
runs = zeros(numel(brk)-1, 2);
for i = 1:numel(brk)-1
    runs(i,1) = idx(brk(i));
    runs(i,2) = brk(i+1) - brk(i);
end
end
