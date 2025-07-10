function bp_data = bandpass_filter_GUI(signal, passband, fs)

% Bandpass filter the signal using zero-phase butter 
[b, a] = butter(4, passband / (fs / 2), 'bandpass');
bp_data = filtfilt(b, a, signal);

end