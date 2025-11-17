function out = kiaSort_detect_spike(inputSignal, cfg, sample)

% Main spike detection function
    if nargin < 3
        sample = 0;
    end

    bandpass_signal    = inputSignal.bandpass_signal;
    scale_factor       = inputSignal.scale_factor;
    rms_bands          = inputSignal.rms_bands;

    trial_length = length(bandpass_signal);
    fs = cfg.samplingFrequency;
    spk_Distance = round(cfg.spikeDistance * fs / 1000);
    original_spk_Distance = spk_Distance;
    border_Margin = round(cfg.borderMargin * fs / 1000);

    % Parameters
    min_thr = cfg.min_threshold;
    badChannel_factor = cfg.badChannel_factor;
    min_rate = cfg.minRate;

    % adjust threshold by scaling factor
    madBP = median(abs(bandpass_signal-median(bandpass_signal)));
    mad_Thresh_init = min_thr * madBP;

    % Initial peak detection
    bp_Loc   = detectPeaks(bandpass_signal, -mad_Thresh_init, 0, -1);
    bp_Loc_r = detectPeaks(bandpass_signal,  mad_Thresh_init, 0,  1);
    initial_rate = max([length(bp_Loc), length(bp_Loc_r)]) / (trial_length / fs);

    % Early termination conditions
    if sample
        if min(rms_bands(1) ./ rms_bands(2:end)) > badChannel_factor
            out = set_null_output(scale_factor, "bad", 0, 0);
            return;
        elseif scale_factor > badChannel_factor
            out = set_null_output(scale_factor, "low", 0, 0);
            return;
        elseif initial_rate < min_rate
            out = set_null_output(scale_factor, "low", 0, 1);
            return;
        elseif scale_factor > 1 && scale_factor < badChannel_factor
            snr_status = "fair"; inclusion_flag = 1;
        else
            snr_status = "good"; inclusion_flag = 1;
        end

    bp_rms_adj = std_exclude(bandpass_signal, [bp_Loc; bp_Loc_r], spk_Distance);
    mad_Thresh = mad_Thresh_init;%min_thr * bp_rms_adj;
    else
        mad_Thresh = .9*inputSignal.mad_Thresh;
        snr_status = "N/A"; inclusion_flag = 1;
    end
    
    Loc   = detectPeaks(bandpass_signal, -mad_Thresh, 0, -1);
    Loc_r = detectPeaks(bandpass_signal, mad_Thresh, 0, 1);

    spk_Loc   = Loc(Loc >= border_Margin & Loc <= trial_length - border_Margin);
    spk_Loc_r = Loc_r(Loc_r >= border_Margin & Loc_r <= trial_length - border_Margin);

    if sample

        spk_Loc = spk_Loc(~ismember(spk_Loc,cfg.edge_exclusion));
        spk_Loc_r = spk_Loc_r(~ismember(spk_Loc_r,cfg.edge_exclusion));

        spike_Vals = abs(bandpass_signal(spk_Loc));
        spike_Vals_r = abs(bandpass_signal(spk_Loc_r));

        extremeVal = prctile(spike_Vals, 90);
        extremeVal_r = prctile(spike_Vals_r, 90);

        spk_dist   = spike_Vals > max(extremeVal,1.3 * mad_Thresh);
        spk_dist_r = spike_Vals_r > max(extremeVal_r,1.3 * mad_Thresh);

        [ext_d, ~] = nearest_index(spk_Loc(spk_dist), spk_Loc_r);
        [ext_d_r, ~] = nearest_index(spk_Loc_r(spk_dist_r), spk_Loc);

        [ext_mod_dist, ext_freq_mod] = mode(ext_d(ext_d > round(.5*spk_Distance)));
        [ext_mod_dist_r, ext_freq_mod_r] = mode(ext_d_r(ext_d_r>round(.5*spk_Distance)));

        mod_dist_ext = 0;                              

        candidates     = [ext_mod_dist   ext_mod_dist_r];
        candidateFreqs = [ext_freq_mod   ext_freq_mod_r];
        distArrays     = {ext_d,         ext_d_r};

        isValid = false(1,2);                     

        for k = 1:2
            d = distArrays{k};

            [u,~,idx] = unique(d);
            cnts      = accumarray(idx,1);        

            otherMask   = (u>round(.5*spk_Distance) & u<=100 & cnts > 0 & u~=candidates(k));
            meanOther   = mean(cnts(otherMask));
            if isnan(meanOther)                   
                meanOther = inf;
            end

            isValid(k) = (candidates(k) <= 2*spk_Distance+1) && ...
                (candidateFreqs(k) > 3*meanOther);
        end

        if any(isValid)
            mod_dist_ext = max(candidates(isValid));  
        end
    end

    if sample
        [d, ~] = nearest_index(spk_Loc, spk_Loc_r);
        if ~isempty(d)
            if sum(d< 2 * spk_Distance) > length(spk_Loc)/20
                d = d(d<100 & d>round(.5*spk_Distance)+1);
                if ~isempty(d)
                    [mod_dist, freq_mod] = mode(d);
                    if freq_mod > (5 * length(d)/(range(d))) && mod_dist < 2 * spk_Distance
                        spk_Distance =  max(spk_Distance, mod_dist + 3);
                    end
                end
            end
        end
        spk_Distance = max([spk_Distance,mod_dist_ext+3]);
    else
        spk_Distance = inputSignal.adj_distant+1;
    end

    spk_ID = [-1*ones(size(spk_Loc)); ones(size(spk_Loc_r))];
    
    if sample
        [spk_All_adj, spk_Val_adj, spk_ID_adj] = keep_highest_value([spk_Loc;spk_Loc_r], abs(bandpass_signal([spk_Loc; spk_Loc_r])), spk_Distance, spk_ID);
    else
        [spk_All_adj, spk_Val_adj, spk_ID_adj] = keep_highest_value([spk_Loc;spk_Loc_r], abs(bandpass_signal([spk_Loc; spk_Loc_r])), spk_Distance, spk_ID, original_spk_Distance);
    end

    rms_values = struct(...
        "bandpass_min_threshold", mad_Thresh, ...
        "bandpass_init_threshold", mad_Thresh_init ...
    );

    out.spk_idx      = spk_All_adj;
    out.spk_Val      = spk_Val_adj;
    out.spk_ID       = spk_ID_adj;
    out.scale_factor = scale_factor;
    out.rms_values   = rms_values;
    out.SNR          = snr_status;
    out.inclusion    = inclusion_flag;
    out.waveInclusion= inclusion_flag;
    out.adj_distant  = spk_Distance;

    if sample
        sample_spike_distance = cfg.spikeSampleDistance * fs /1000;
        left_diff       = [Inf; diff(spk_All_adj)];
        right_diff      = [diff(spk_All_adj); Inf];
        valid_idx_dist  = (left_diff > sample_spike_distance) & (right_diff > sample_spike_distance);
        numValidSpikes  = sum(valid_idx_dist);
        numSpikes       = length(spk_All_adj);
        if (numSpikes/ (trial_length / fs)) < min_rate || numValidSpikes <= 2 * cfg.nPCAcomp
            out = set_null_output(scale_factor, "low", 0, 1);
        end
    end
end

function s = std_exclude(signal, Idx, N)
% std excluding samples around initial spikes
    idxExclude = unique(Idx(:) + (-N:N));
    idxExclude = idxExclude(idxExclude >= 1 & idxExclude <= length(signal));
    include = true(size(signal));
    include(idxExclude) = false;
    s = median(abs(signal(include)-median(signal(include))));
end

function out = set_null_output(scale_factor, SNR_status, inclusion_flag, waveInc_flag)
% null spike data
    out.spk_idx       = [];
    out.spk_Val       = [];
    out.spk_ID        = [];
    out.scale_factor  = scale_factor;
    out.rms_values    = [];
    out.SNR           = SNR_status;
    out.inclusion     = inclusion_flag;
    out.waveInclusion = waveInc_flag;
    out.adj_distant   = 0;
end