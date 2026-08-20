function [corrThr, ampThr] = kiaSort_channel_scaled_thresholds(nCh, corrThr0, ampThr0, cfg)
%KIASORT_CHANNEL_SCALED_THRESHOLDS  Tighten merge gates on narrow footprints.
%
%   [corrThr, ampThr] = kiaSort_channel_scaled_thresholds(nCh, corrThr0, ampThr0, cfg)
%
%   Two different neurons correlate highly on a single channel far more
%   often than across a wide footprint, so the same correlation is much
%   weaker evidence when nCh is small. Below corrScaleRefChannels the
%   correlation gate is raised and the amplitude-variance gate lowered;
%   at or above it both are returned unchanged.
%
%   cfg fields (all optional):
%       corrScaleRefChannels  reference footprint width (default 9)
%       mergeScaleStrength    0..1, how hard to tighten both gates (0.25; 0 disables)

if nargin < 4, cfg = struct(); end

nRef  = local_get(cfg, 'corrScaleRefChannels', 9);
kScale = local_get(cfg, 'mergeScaleStrength', 0.25);

corrThr = corrThr0;
ampThr  = ampThr0;
if ~isfinite(nCh) || nCh <= 0 || nRef <= 0
    return;
end

f = min(1, nCh / nRef);
if f >= 1
    return;
end

if ~isempty(corrThr0) && isfinite(corrThr0)
    corrThr = corrThr0 + (1 - corrThr0) * kScale * (1 - f);
end
if ~isempty(ampThr0) && isfinite(ampThr0)
    ampThr = ampThr0 * (1 - kScale * (1 - f));
end

end


function v = local_get(cfg, name, dflt)
if isstruct(cfg) && isfield(cfg, name) && ~isempty(cfg.(name))
    v = cfg.(name);
else
    v = dflt;
end
end
