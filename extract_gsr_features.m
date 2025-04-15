function gsr_features = extract_gsr_features(gsr_signal, fs)
% Inputs:
%   gsr_signal : 1D array of GSR data
%   fs         : Sampling rate
%
% Output:
%   gsr_features : Struct of features

    gsr_signal = detrend(gsr_signal);  % Remove linear drift

    % Basic stats
    gsr_features.mean     = mean(gsr_signal);
    gsr_features.std      = std(gsr_signal);
    gsr_features.range    = range(gsr_signal);
    gsr_features.min      = min(gsr_signal);
    gsr_features.max      = max(gsr_signal);

    % Area under the curve
    gsr_features.auc      = trapz(gsr_signal);

    % Number of peaks (simple SCR detection)
    [pks, locs] = findpeaks(gsr_signal, 'MinPeakHeight', std(gsr_signal), ...
                            'MinPeakDistance', fs);  % SCRs >1 sec apart
    gsr_features.num_scrs = length(pks);
    gsr_features.mean_scr_amplitude = mean(pks);

    % Rise time estimation (approx.)
    rise_times = [];
    for i = 1:length(pks)
        start_idx = max(locs(i) - round(1*fs), 1);
        baseline = gsr_signal(start_idx);
        rise_times(i) = (pks(i) - baseline) / fs;  % crude rise time
    end
    gsr_features.mean_rise_time = mean(rise_times);

    % Frequency-domain: low frequency power (<0.5Hz)
    [pxx, f] = pwelch(gsr_signal, hamming(fs*2), [], [], fs);
    gsr_features.low_freq_power = bandpower(pxx, f, [0.01 0.5], 'psd');
end