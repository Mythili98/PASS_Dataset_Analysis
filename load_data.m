clear all;
close all;
clc;
%% add eeglab to the path
addpath('/home/desktop/Desktop/22104412_Docs/MATLAB/EEGLAB/eeglab2023.1');
eeglab;
%% Processing
file_path = '/home/desktop/Desktop/22104412_Docs/Missing_Modalities/PASS_stress/';
output_file_path = '/home/desktop/Desktop/22104412_Docs/Missing_Modalities/PASS_stress/feat/';
files = dir(fullfile(file_path, '*.mat'));
for fi=1:length(files)
    filename = files(fi);
    f_name = [filename.folder,filesep,filename.name];
    name = filename.name;
    %% Load data
    data = load(f_name);
    eeg_data = data.eeg_data.eeg;
    ecg_data = data.ecg_data.ecg;
    gsr_data = data.gsr_data.gsr;
    resp_data = data.ecg_data.resp;
    %% Clean EEG and extract features
    fs = eeg_data.fs; % Sampling rate = 220Hz
    eegData = eeg_data.X';
    EEG = pop_importdata('dataformat','array','data','eegData','srate',fs); % Modify 'srate' accordingly
    % Manually assign channel labels
    custom_labels = {'TP9', 'AF7', 'AF8', 'TP10'};
    for i = 1:length(custom_labels)
        EEG.chanlocs(i).labels = custom_labels{i};
    end
    EEG = pop_chanedit(EEG, 'lookup','standard-10-5-cap385.elp');
    % Band pass
    EEG = pop_eegfiltnew(EEG, 1, 40);
    eegData = extract_trials(eeg_data);
    
    for i=1:length(eegData)
        data = eegData(i).data;
        EEG = pop_importdata('dataformat','array','data','data','srate',fs);
        EEG = pop_resample(EEG, 128);
        data = EEG.data;
        % Wavelet ICA
        [icaEEG, A, W] = fastica(data,'stabilization','off','verbose','off'); 
        nICs = 1:size(icaEEG,1); 
        Kthr = 1;         
        ArtefThreshold = 7; 
        verbose = 'off';         
        
        icaEEG2 = RemoveStrongArtifacts(icaEEG, nICs, Kthr, ArtefThreshold, fs, verbose); 
        Data_wICA = A*icaEEG2;
        eegData(i).data = Data_wICA;
    end
    %segment into 1 minute data and remove the extra if needed from the start
    %of the trial
    eegData = segment_epochs(eegData,128, 60, 0, []);
    %% Clean ECG and extract features
    % Pan-Tompkins Algorithm
    ecgData = extract_trials(ecg_data);
    ecg_fs = ecg_data.fs;
    %segment data - 1 minute epochs - non overlapping
    ecgData = segment_epochs(ecgData,ecg_fs, 60, 0, []);
    ecg = [];
    for i = 1:length(ecgData)
        [ecg(i).qrs_amp_raw,ecg(i).qrs_i_raw,ecg(i).Q_peaks,ecg(i).Q_peaks_val,ecg(i).S_peaks,ecg(i).S_peaks_val,ecg(i).T_peaks,ecg(i).T_peaks_val,ecg(i).delay] = pan_tompkin(ecgData(i).data, 250, 0);
        ecg(i).label = ecgData(i).labels;
        ecg(i).y = ecgData(i).y;
    end
    
    %% Clean GSR and extract features
    % Parameters
    Fs = gsr_data.fs;               % Sampling frequency (change to your actual Fs)
    Fc = 1;                % Cutoff frequency (Hz)
    Rp = 0.5;              % Passband ripple (dB)
    N = 8;                 % Filter order
    Wn = Fc / (Fs/2);
    [b, a] = cheby1(N, Rp, Wn, 'low');
    gsrfeat = {};
    % Apply zero-phase filtering to avoid phase distortion
    gsr_filtered = filtfilt(b, a, gsr_data.X');
    % extract segments
    gsrData = extract_trials(gsr_data);
    gsrData = segment_epochs(gsrData, Fs, 60, 0, []);
    for i=1:length(gsrData)
        data = gsrData(i).data;
        label = gsrData(i).labels;
        fs = gsrData(i).fs;
        gsr = data(1,:);
        gsrfeat(i).feature =  extract_gsr_features(gsr, fs);
        gsrfeat(i).label = label;
        gsrfeat(i).y = gsrData(i).y;
    end
    %% Extract EEG Features
    fs = 128;  % sampling frequency
    eegfeat = [];
    for i = 1:length(eegData)
        eeg = eegData(i).data;  % should be 4 x 30000
        
        % Define bands
        bands = struct( ...
            'delta', [1 4], ...
            'theta', [4 8], ...
            'alpha', [8 13], ...
            'beta',  [13 30], ...
            'gamma', [30 45]);
        alpha_asymmetry_pairs = { {'AF7', 'AF8'}, {'TP9', 'TP10'} };
        
        n_channels = size(eeg, 1);
        band_names = fieldnames(bands);
        
        for ch = 1:n_channels
            signal = eeg(ch, :);
            ch_name = custom_labels{ch};
            % Welch PSD estimation
            win_len = 1 * fs;  % 2-second window
            [pxx, f] = pwelch(signal, hamming(win_len), win_len/2, [], fs);
        
            total_pow = trapz(f, pxx);  % total power for relative calc
        
            for b = 1:length(band_names)
                band = bands.(band_names{b});
                idx = f >= band(1) & f <= band(2);
                eegfeat(i).abs_power.(ch_name).(band_names{b}) = trapz(f(idx), pxx(idx));
                eegfeat(i).abs_power_dB.(ch_name).(band_names{b}) = 10* log10(trapz(f(idx), pxx(idx)));
                eegfeat(i).rel_power.(ch_name).(band_names{b}) = eegfeat(i).abs_power.(ch_name).(band_names{b})/ total_pow;
            end
        end
    
        for pair_idx = 1:length(alpha_asymmetry_pairs)
        left_ch = alpha_asymmetry_pairs{pair_idx}{1};
        right_ch = alpha_asymmetry_pairs{pair_idx}{2};
        
        % Get the absolute power for each channel in the alpha band
        left_alpha_power = eegfeat(i).abs_power.(left_ch).alpha;
        right_alpha_power = eegfeat(i).abs_power.(right_ch).alpha;
        
        % Calculate the asymmetry
        eegfeat(i).alpha_asymmetry.(sprintf('%s_%s', left_ch, right_ch)) = ...
            (left_alpha_power - right_alpha_power) / (left_alpha_power + right_alpha_power);
    end
    % Calculate Coherence for different bands (e.g., alpha, beta, gamma, theta)
    for band_idx = 1:length(band_names)
        band_range = bands.(band_names{band_idx});
        
        % Coherence between TP9 and TP10
        signal_1 = eeg(strcmp(custom_labels, 'TP9'), :);  % TP9
        signal_2 = eeg(strcmp(custom_labels, 'TP10'), :);  % TP10
        
        % Calculate coherence using mscohere
        [Cxy, f_cohere] = mscohere(signal_1, signal_2, hamming(2*fs), fs, [], fs);
        
        % Store coherence values for the given frequency band
        idx_band = f_cohere >= band_range(1) & f_cohere <= band_range(2);
        eegfeat(i).coherence_results.(band_names{band_idx}) = mean(Cxy(idx_band));  % Average coherence in the band
    end
        eegfeat(i).label = eegData(i).labels;
        eegfeat(i).y = eegData(i).y;
    end
    
    %% Saving data
    % ecg, gsr and eeg feature structs
    file_name_feat = [output_file_path,name(1:4),'_feat.mat'];
    file_name_segments = [output_file_path,name(1:4),'_segment.mat'];
    save(file_name_feat, 'eegfeat', 'ecg', 'gsrfeat');
    save(file_name_segments,'eegData', 'ecgData', 'gsrData');
end