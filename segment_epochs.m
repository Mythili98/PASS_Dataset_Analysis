function epoched_data = segment_epochs(signal_data, fs, epoch_seconds, overlapping, overlaping_ratio)
%SEGMENT_EPOCHS Segments data into fixed-length epochs with optional overlap.
%
% Inputs:
%   signal_data : struct array with fields 'data', 'labels', and optionally 'fs'
%   epoch_seconds : number of seconds in one epoch
%   overlapping : 1 (enable overlap) or 0 (no overlap)
%   overlaping_ratio : [optional] ratio of overlap (default = 0.5)
%
% Output:
%   epoched_data : struct array with segmented data and labels

    if nargin < 4
        overlaping_ratio = 0.5;  % default overlap ratio
    end

    k = 1;  % counter for epoched_data
    epoched_data = struct();

    for i = 1:length(signal_data)
        data = signal_data(i).data;
        labels = signal_data(i).label;
        
        epoch_len = epoch_seconds * fs;
        
        if overlapping
            step = round(epoch_len * (1 - overlaping_ratio));
        else
            step = epoch_len;
        end

        num_samples = size(data, 2);
        
        for start_idx = 1:step:(num_samples - epoch_len + 1)
            end_idx = start_idx + epoch_len - 1;

            epoched_data(k).data = data(:, start_idx:end_idx);
            epoched_data(k).labels = labels;
            epoched_data(k).fs = fs;
            k = k + 1;
        end
    end
end
