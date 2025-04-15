function trial_structs = extract_trials(signal_data)
    markers = signal_data.markers;
    data = signal_data.X';  % channels x samples

    [marker_indices, ~, marker_values] = find(markers);

    trials = {};
    i = 1;
    latest_practice = [];
    latest_practice_idx = [];

    while i <= length(marker_values) - 1
        val_i = marker_values(i);
        val_next = marker_values(i+1);

        if val_i == 10 && val_next == 20
            latest_practice = [val_i, val_next];
            latest_practice_idx = [marker_indices(i), marker_indices(i+1)];
            i = i + 2;
            continue;
        end

        if val_i >= 11 && val_next == val_i + 10
            repeated_trials = {};
            code = val_i;

            while i <= length(marker_values) - 1 && ...
                  marker_values(i) == code && ...
                  marker_values(i+1) == code + 10

                idx_start = marker_indices(i);
                idx_end = marker_indices(i+1);

                repeated_trials{end+1} = struct( ...
                    'practice', latest_practice, ...
                    'practice_idx', latest_practice_idx, ...
                    'trial', [val_i, val_next], ...
                    'trial_idx', [idx_start, idx_end], ...
                    'trial_code', code, ...
                    'span', idx_end - idx_start ...
                );

                i = i + 2;
                if i <= length(marker_values) - 1
                    val_i = marker_values(i);
                    val_next = marker_values(i+1);
                end
            end

            if ~isempty(repeated_trials)
                [~, max_idx] = max(cellfun(@(t) t.span, repeated_trials));
                trials{end+1} = repeated_trials{max_idx};
                latest_practice = [];
                latest_practice_idx = [];
            end
        else
            i = i + 1;
        end
    end

trial_structs = struct('data', {}, 'label', {});  % Proper initialization

for k = 1:length(trials)
    % Practice (tutorial) data
    if ~isempty(trials{k}.practice_idx)
        idx_range = trials{k}.practice_idx(1):trials{k}.practice_idx(2);
        entry.data = data(:, idx_range);
        entry.label = trials{k}.trial_code * 10;  % e.g., 110
        entry.y = [];
        trial_structs(end+1) = entry;
    end

    % Trial data
    idx_range = trials{k}.trial_idx(1):trials{k}.trial_idx(2);
    entry.data = data(:, idx_range);
    entry.label = trials{k}.trial_code * 10 + 1;  % e.g., 111
    entry.y = [];
    trial_structs(end+1) = entry;
end

fprintf('Extracted %d trials (returned %d struct entries).\n', ...
    length(trials), length(trial_structs));
end
