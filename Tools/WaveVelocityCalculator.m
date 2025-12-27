% This program is used to compare multiple CFPP samples and attempt to find
% out what the velocity of the material should be.
% Kostiantyn Savuk | 06/13/2025

%% Transducer Properties
OutputFreq_Mhz = 10;
SamplingRate_Mhz = 100; % 10x the output frequency

%% CFPP Material Properties (Those important for Wave Velocity)
Velocity_Range = [2000, 3500];

%% Front Wall Detection Parameters
% Estimated front wall time ranges (in sample points)
FrontWall_Range_Samples = [4, 15]; % Typical range for front wall detection

%% Data Quality Control
% Samples 3, 4, 5 have purposeful defects that may affect layer detection
% Enable data filtering to test with/without problematic samples
Enable_Data_Filtering = true;
Problematic_Samples = [3, 4, 5]; % Sample indices with known defects



%% Sample 1: L0P5S25_10Mhz
% Number of Plys: 6
% Overall Thickness (mm): 1.146
% Average Layer Thickness (mm): 0.191

% Aprox. Wall time seen from B scans (10e-6seconds)
% Visual Layers Seen in Heatmaps @ times
Layer1_1        = [17,21];
Layer1_1_Value  = [1];
Layer1_2        = [22,28];
Layer1_2_Value  = [-1];
Layer1_3        = [29,36];
Layer1_3_Value  = [1];
Layer1_4        = [37,44];
Layer1_4_Value  = [-1];
Layer1_5        = [45,53];
Layer1_5_Value  = [1];
Layer1_6        = [54,60];
Layer1_6_Value  = [-1];
Layer1_7        = [61,67];
Layer1_7_Value  = [1];
Layer1_8        = [68,74];
Layer1_8_Value  = [-1];
Layer1_9        = [75,81];
Layer1_9_Value  = [1];
Layer1_10       = [82,87];
Layer1_10_Value = [-1];
Layer1_11       = [88,92];
Layer1_11_Value = [1];
Layer1_12       = [93,98];
Layer1_12_Value = [-1];
Layer1_13       = [99,105];
Layer1_13_Value = [1];
Layer1_14       = [106,113];
Layer1_14_Value = [-1];
Layer1_15       = [114,122];
Layer1_15_Value = [1];
Layer1_16       = [123,133];
Layer1_16_Value = [-1];
Layer1_17       = [134,141];
Layer1_17_Value = [1];
Layer1_18       = [142,147];
Layer1_18_Value = [-1];
Layer1_19       = [148,152];
Layer1_19_Value = [1];
Layer1_20       = [153,159];
Layer1_20_Value = [-1];
Layer1_21       = [160,166];
Layer1_21_Value = [1];
Layer1_22       = [167,173];
Layer1_22_Value = [-1];
Layer1_23       = [174,183];
Layer1_23_Value = [1];
Layer1_24       = [184,192];
Layer1_24_Value = [-1];
Layer1_25       = [193,200];
Layer1_25_Value = [1];
Layer1_26       = [201,210];
Layer1_26_Value = [-1];

%% Sample 2: L8P12S5_10Mhz
% Number of Plys: N/A
% Overall Thickness: N/A
% Average Layer Thickness: N/A

% Aprox. Wall time seen from B scans (10e-6seconds)
FrontWallTime2 = 0;
BackWallTime2  = 0;

%% Sample 3: L16P1S9_10MHz
% Number of Plys: 10
% Overall Thickness (mm): 1.805
% Average Layer Thickness (mm): 0.168

% Aprox. Wall time seen from B scans (10e-6seconds)
% Visual Layers Seen in Heatmaps @ times
Layer3_1        = [2,8];
Layer3_1_Value  = [-1];
Layer3_2        = [9,16];
Layer3_2_Value  = [1];
Layer3_3        = [17,25];
Layer3_3_Value  = [-1];
Layer3_4        = [26,33];
Layer3_4_Value  = [1];
Layer3_5        = [33,40];
Layer3_5_Value  = [-1];
Layer3_6        = [41,46];
Layer3_6_Value  = [1];
Layer3_7        = [47,57];
Layer3_7_Value  = [-1];
Layer3_8        = [58,66];
Layer3_8_Value  = [1];
Layer3_9        = [67,72];
Layer3_9_Value  = [-1];
Layer3_10       = [73,79];
Layer3_10_Value = [1];
Layer3_11       = [80,86];
Layer3_11_Value = [-1];
Layer3_12       = [87,92];
Layer3_12_Value = [1];
Layer3_13       = [93,99];
Layer3_13_Value = [-1];
Layer3_14       = [100,106];
Layer3_14_Value = [1];
Layer3_15       = [107,112];
Layer3_15_Value = [-1];
Layer3_16       = [113,119];
Layer3_16_Value = [1];
Layer3_17       = [120,125];
Layer3_17_Value = [-1];
Layer3_18       = [126,132];
Layer3_18_Value = [1];
Layer3_19       = [133,139];
Layer3_19_Value = [-1];
Layer3_20       = [140,149];
Layer3_20_Value = [1];
Layer3_21       = [150,157];
Layer3_21_Value = [-1];
Layer3_22       = [158,166];
Layer3_22_Value = [1];
Layer3_23       = [167,174];
Layer3_23_Value = [-1];
Layer3_24       = [175,182];
Layer3_24_Value = [1];
Layer3_25       = [183,189];
Layer3_25_Value = [-1];


%% Sample 4: L16P1S10_10MHz
% Number of Plys: 10
% Overall Thickness (mm): 1.738
% Average Layer Thickness (mm): 0.162

% Aprox. Wall time seen from B scans (10e-6seconds)
% Visual Layers Seen in Heatmaps @ times
Layer4_1        = [4,9];
Layer4_1_Value  = [-1];
Layer4_2        = [10,17];
Layer4_2_Value  = [1];
Layer4_3        = [18,26];
Layer4_3_Value  = [-1];
Layer4_4        = [27,34];
Layer4_4_Value  = [1];
Layer4_5        = [35,41];
Layer4_5_Value  = [-1];
Layer4_6        = [42,47];
Layer4_6_Value  = [1];
Layer4_7        = [48,60];
Layer4_7_Value  = [-1];
Layer4_8        = [61,67];
Layer4_8_Value  = [1];
Layer4_9        = [68,72];
Layer4_9_Value  = [-1];
Layer4_10       = [73,80];
Layer4_10_Value = [1];
Layer4_11       = [81,86];
Layer4_11_Value = [-1];
Layer4_12       = [87,91];
Layer4_12_Value = [1];
Layer4_13       = [92,98];
Layer4_13_Value = [-1];
Layer4_14       = [99,102];
Layer4_14_Value = [1];
Layer4_15       = [103,109];
Layer4_15_Value = [-1];
Layer4_16       = [110,119];
Layer4_16_Value = [1];
Layer4_17       = [120,125];
Layer4_17_Value = [-1];
Layer4_18       = [126,131];
Layer4_18_Value = [1];
Layer4_19       = [132,139];
Layer4_19_Value = [-1];
Layer4_20       = [140,148];
Layer4_20_Value = [1];
Layer4_21       = [149,157];
Layer4_21_Value = [-1];
Layer4_22       = [158,165];
Layer4_22_Value = [1];
Layer4_23       = [166,173];
Layer4_23_Value = [-1];
Layer4_24       = [174,181];
Layer4_24_Value = [1];
Layer4_25       = [182,188];
Layer4_25_Value = [-1];
Layer4_26       = [189,195];
Layer4_26_Value = [1];

%% Sample 5: L16P34S3_10MHz
% Number of Plys: N/A
% Overall Thickness (mm): N/A
% Average Layer Thickness (mm): N/A

% Aprox. Wall time seen from B scans (10e-6seconds)
% Visual Layers Seen in Heatmaps @ times
Layer5_1        = [4,9];
Layer5_1_Value  = [-1];
Layer5_2        = [10,17];
Layer5_2_Value  = [1];
Layer5_3        = [18,26];
Layer5_3_Value  = [-1];
Layer5_4        = [27,34];
Layer5_4_Value  = [1];
Layer5_5        = [35,41];
Layer5_5_Value  = [-1];
Layer5_6        = [42,47];
Layer5_6_Value  = [1];
Layer5_7        = [48,54];
Layer5_7_Value  = [-1];
Layer5_8        = [55,60];
Layer5_8_Value  = [1];
Layer5_9        = [61,64];
Layer5_9_Value  = [-1];
Layer5_10       = [65,70];
Layer5_10_Value = [1];
Layer5_11       = [71,76];
Layer5_11_Value = [-1];
Layer5_12       = [77,83];
Layer5_12_Value = [1];
Layer5_13       = [84,90];
Layer5_13_Value = [-1];
Layer5_14       = [91,96];
Layer5_14_Value = [1];
Layer5_15       = [97,102];
Layer5_15_Value = [-1];
Layer5_16       = [103,107];
Layer5_16_Value = [1];
Layer5_17       = [108,114];
Layer5_17_Value = [-1];
Layer5_18       = [115,119];
Layer5_18_Value = [1];
Layer5_19       = [120,125];
Layer5_19_Value = [-1];
Layer5_20       = [126,136];
Layer5_20_Value = [1];
Layer5_21       = [137,145];
Layer5_21_Value = [-1];
Layer5_22       = [146,155];
Layer5_22_Value = [1];
Layer5_23       = [156,164];
Layer5_23_Value = [-1];
Layer5_24       = [165,175];
Layer5_24_Value = [1];
Layer5_25       = [176,184];
Layer5_25_Value = [-1];
Layer5_26       = [185,191];
Layer5_26_Value = [1];

%% Data Organization and Processing
% Convert manual data into structured format for easier processing

% Initialize sample data structure
samples = struct();

% Sample 1: L0P5S25_10Mhz
samples(1).name = 'L0P5S25_10Mhz';
samples(1).num_plys = 6;
samples(1).thickness_mm = 1.146;
samples(1).avg_layer_thickness_mm = 0.191;
samples(1).layers = [];
samples(1).values = [];

% Collect Sample 1 data
layer_times_1 = [17,21; 22,28; 29,36; 37,44; 45,53; 54,60; 61,67; 68,74; 75,81; 82,87; 88,92; 93,98; 99,105; 106,113; 114,122; 123,133; 134,141; 142,147; 148,152; 153,159; 160,166; 167,173; 174,183; 184,192; 193,200; 201,210];
layer_values_1 = [1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1];
samples(1).layers = layer_times_1;
samples(1).values = layer_values_1;

% Sample 3: L16P1S9_10MHz
samples(3).name = 'L16P1S9_10MHz';
samples(3).num_plys = 10;
samples(3).thickness_mm = 1.805;
samples(3).avg_layer_thickness_mm = 0.168;
samples(3).layers = [];
samples(3).values = [];

% Collect Sample 3 data
layer_times_3 = [2,8; 9,16; 17,25; 26,33; 33,40; 41,46; 47,57; 58,66; 67,72; 73,79; 80,86; 87,92; 93,99; 100,106; 107,112; 113,119; 120,125; 126,132; 133,139; 140,149; 150,157; 158,166; 167,174; 175,182; 183,189];
layer_values_3 = [-1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1];
samples(3).layers = layer_times_3;
samples(3).values = layer_values_3;

% Sample 4: L16P1S10_10MHz
samples(4).name = 'L16P1S10_10MHz';
samples(4).num_plys = 10;
samples(4).thickness_mm = 1.738;
samples(4).avg_layer_thickness_mm = 0.162;
samples(4).layers = [];
samples(4).values = [];

% Collect Sample 4 data
layer_times_4 = [4,9; 10,17; 18,26; 27,34; 35,41; 42,47; 48,60; 61,67; 68,72; 73,80; 81,86; 87,91; 92,98; 99,102; 103,109; 110,119; 120,125; 126,131; 132,139; 140,148; 149,157; 158,165; 166,173; 174,181; 182,188; 189,195];
layer_values_4 = [-1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1];
samples(4).layers = layer_times_4;
samples(4).values = layer_values_4;

% Sample 5: L16P34S3_10MHz
samples(5).name = 'L16P34S3_10MHz';
samples(5).num_plys = NaN; % Not available
samples(5).thickness_mm = NaN; % Not available
samples(5).avg_layer_thickness_mm = NaN; % Not available
samples(5).layers = [];
samples(5).values = [];

% Collect Sample 5 data
layer_times_5 = [4,9; 10,17; 18,26; 27,34; 35,41; 42,47; 48,54; 55,60; 61,64; 65,70; 71,76; 77,83; 84,90; 91,96; 97,102; 103,107; 108,114; 115,119; 120,125; 126,136; 137,145; 146,155; 156,164; 165,175; 176,184; 185,191];
layer_values_5 = [-1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1];
samples(5).layers = layer_times_5;
samples(5).values = layer_values_5;

% Note: Sample 2 has no data, so we skip it

%% Manual TOF and Layer Timing Calculations

% Calculate expected TOF and layer timing based on known material properties
fprintf('\n=== MANUAL TOF CALCULATIONS ===\n');
fprintf('Using estimated velocity range: %.0f - %.0f m/s\n', Velocity_Range(1), Velocity_Range(2));

% For each sample with known thickness, calculate expected timing
manual_results = struct();
for i = 1:length(samples)
    if ~isempty(samples(i).name) && ~isnan(samples(i).thickness_mm)
        thickness_m = samples(i).thickness_mm * 1e-3;

        % Calculate expected TOF range
        tof_min_s = (2 * thickness_m) / Velocity_Range(2); % Faster velocity = shorter TOF
        tof_max_s = (2 * thickness_m) / Velocity_Range(1); % Slower velocity = longer TOF

        tof_min_samples = tof_min_s * (SamplingRate_Mhz * 1e6);
        tof_max_samples = tof_max_s * (SamplingRate_Mhz * 1e6);

        % Calculate expected layer spacing (if layers are uniform)
        if ~isnan(samples(i).num_plys)
            layer_spacing_m = thickness_m / samples(i).num_plys;
            layer_spacing_tof_s = (2 * layer_spacing_m) / mean(Velocity_Range);
            layer_spacing_samples = layer_spacing_tof_s * (SamplingRate_Mhz * 1e6);
        else
            layer_spacing_samples = NaN;
        end

        manual_results(i).name = samples(i).name;
        manual_results(i).thickness_mm = samples(i).thickness_mm;
        manual_results(i).expected_tof_range_samples = [tof_min_samples, tof_max_samples];
        manual_results(i).expected_tof_range_us = [tof_min_samples, tof_max_samples] / SamplingRate_Mhz;
        manual_results(i).expected_layer_spacing_samples = layer_spacing_samples;

        fprintf('\nSample %s (%.3f mm, %d plys):\n', samples(i).name, samples(i).thickness_mm, samples(i).num_plys);
        fprintf('  Expected TOF: %.1f - %.1f samples (%.2f - %.2f μs)\n', ...
            tof_min_samples, tof_max_samples, tof_min_samples/SamplingRate_Mhz, tof_max_samples/SamplingRate_Mhz);
        if ~isnan(layer_spacing_samples)
            fprintf('  Expected layer spacing: %.1f samples (%.2f μs)\n', ...
                layer_spacing_samples, layer_spacing_samples/SamplingRate_Mhz);
        end
    end
end

%% Data Quality Assessment and Filtering

% Function to assess data quality and optionally filter problematic samples
function [filtered_samples, quality_report] = assess_and_filter_data(samples, enable_filtering, problematic_indices)
    quality_report = struct();
    filtered_samples = samples;

    fprintf('\n=== DATA QUALITY ASSESSMENT ===\n');

    for i = 1:length(samples)
        if ~isempty(samples(i).name)
            layer_values = samples(i).values;
            layer_times = samples(i).layers;

            % Calculate quality metrics
            num_transitions = sum(abs(diff(layer_values)) > 0);
            total_layers = length(layer_values);
            transition_rate = num_transitions / total_layers;

            % Check for excessive transitions (sign of defects)
            % For alternating layers, expect ~50% transitions, but allow more flexibility
            expected_transitions = total_layers * 0.5;
            excess_transitions = num_transitions > expected_transitions * 2.0; % More lenient threshold

            % Check for irregular layer spacing
            layer_centers = mean(layer_times, 2);
            layer_spacings = diff(layer_centers);
            spacing_cv = std(layer_spacings) / mean(layer_spacings); % Coefficient of variation
            irregular_spacing = spacing_cv > 0.8; % More lenient threshold for spacing variation

            quality_report(i).name = samples(i).name;
            quality_report(i).num_transitions = num_transitions;
            quality_report(i).transition_rate = transition_rate;
            quality_report(i).spacing_cv = spacing_cv;
            quality_report(i).has_excess_transitions = excess_transitions;
            quality_report(i).has_irregular_spacing = irregular_spacing;
            quality_report(i).is_problematic = ismember(i, problematic_indices);
            % Only filter if sample is in problematic list AND has quality issues
            quality_report(i).recommended_filter = ismember(i, problematic_indices) && (excess_transitions || irregular_spacing);

            fprintf('Sample %s:\n', samples(i).name);
            fprintf('  Transitions: %d/%d (%.1f%%)\n', num_transitions, total_layers, transition_rate*100);
            fprintf('  Spacing CV: %.2f (%.1f%% variation)\n', spacing_cv, spacing_cv*100);
            % Build quality issues list
            issues = {};
            if excess_transitions
                issues{end+1} = 'Excess transitions';
            end
            if irregular_spacing
                issues{end+1} = 'Irregular spacing';
            end
            if ismember(i, problematic_indices)
                issues{end+1} = 'Known defects';
            end

            if isempty(issues)
                issues_str = 'None';
            else
                issues_str = strjoin(issues, ', ');
            end

            fprintf('  Quality issues: %s\n', issues_str);

            if quality_report(i).recommended_filter && enable_filtering
                fprintf('  -> FILTERED OUT due to quality issues\n');
                filtered_samples(i).name = ''; % Mark for exclusion
            end
        end
    end
end

% Apply data filtering
[filtered_samples, quality_report] = assess_and_filter_data(samples, Enable_Data_Filtering, Problematic_Samples);

%% Velocity Calculations

% Function to calculate velocity based on time-of-flight and thickness
% Formula: v = (2 * thickness) / time_of_flight
% Convert: thickness_mm to meters (*1e-3), time_samples to seconds (samples/sampling_rate)
calculate_velocity = @(time_samples, thickness_mm) (2 * thickness_mm * 1e-3) / (time_samples / (SamplingRate_Mhz * 1e6)); % m/s

% Initialize results structure
results = struct();
velocity_range = Velocity_Range; % Use the defined range

% Process each sample with known thickness (using filtered data)
valid_samples = [];
for i = 1:length(filtered_samples)
    if ~isempty(filtered_samples(i).name) && ~isnan(filtered_samples(i).thickness_mm)
        valid_samples = [valid_samples, i];

        % Calculate center times for each layer
        layer_centers = mean(filtered_samples(i).layers, 2);

        % Find the approximate back wall time (last significant layer)
        back_wall_time = max(layer_centers);

        % Calculate velocity based on total thickness and back wall time
        velocity = calculate_velocity(back_wall_time, filtered_samples(i).thickness_mm);

        % Store results
        results(i).name = filtered_samples(i).name;
        results(i).thickness_mm = filtered_samples(i).thickness_mm;
        results(i).back_wall_time_us = back_wall_time;
        results(i).calculated_velocity_ms = velocity;
        results(i).layer_centers = layer_centers;
        results(i).layer_values = filtered_samples(i).values;

        fprintf('Sample %s:\n', filtered_samples(i).name);
        fprintf('  Thickness: %.3f mm\n', filtered_samples(i).thickness_mm);
        fprintf('  Back wall time: %.1f samples (%.2f μs)\n', back_wall_time, back_wall_time/(SamplingRate_Mhz));
        fprintf('  Calculated velocity: %.0f m/s\n\n', velocity);
    end
end

% Calculate statistics for velocity estimates
velocities = [results(valid_samples).calculated_velocity_ms];
mean_velocity = mean(velocities);
std_velocity = std(velocities);
median_velocity = median(velocities);

fprintf('=== VELOCITY ANALYSIS SUMMARY ===\n');
fprintf('Mean velocity: %.0f ± %.0f m/s\n', mean_velocity, std_velocity);
fprintf('Median velocity: %.0f m/s\n', median_velocity);
fprintf('Range: %.0f - %.0f m/s\n', min(velocities), max(velocities));
fprintf('Coefficient of variation: %.1f%%\n', (std_velocity/mean_velocity)*100);

% Test velocity range to find optimal value
test_velocities = linspace(velocity_range(1), velocity_range(2), 100);
error_scores = zeros(size(test_velocities));

for v_idx = 1:length(test_velocities)
    test_vel = test_velocities(v_idx);
    total_error = 0;

    for s_idx = valid_samples
        % Calculate expected back wall time for this velocity (in sample points)
        expected_time = (2 * results(s_idx).thickness_mm * 1e-3) / test_vel * (SamplingRate_Mhz * 1e6);
        actual_time = results(s_idx).back_wall_time_us;

        % Calculate relative error
        rel_error = abs(expected_time - actual_time) / actual_time;
        total_error = total_error + rel_error^2;
    end

    error_scores(v_idx) = sqrt(total_error / length(valid_samples));
end

% Find optimal velocity
[min_error, opt_idx] = min(error_scores);
optimal_velocity = test_velocities(opt_idx);

fprintf('\n=== OPTIMAL VELOCITY ANALYSIS ===\n');
fprintf('Optimal velocity: %.0f m/s\n', optimal_velocity);
fprintf('RMS error: %.3f\n', min_error);

%% Enhanced TOF Analysis (inspired by the reference script)

% Calculate Time-of-Flight (TOF) for each sample
fprintf('\n=== TIME-OF-FLIGHT ANALYSIS ===\n');
for i = 1:length(valid_samples)
    s_idx = valid_samples(i);

    % Estimate front wall time (first significant layer or minimum range)
    layer_centers = results(s_idx).layer_centers;
    front_wall_estimate = min(layer_centers); % First layer center
    front_wall_estimate = max(front_wall_estimate, FrontWall_Range_Samples(1)); % Apply minimum constraint

    % Calculate TOF
    back_wall_time = results(s_idx).back_wall_time_us;
    tof_samples = back_wall_time - front_wall_estimate;
    tof_us = tof_samples / SamplingRate_Mhz;

    % Store enhanced results
    results(s_idx).front_wall_time_samples = front_wall_estimate;
    results(s_idx).tof_samples = tof_samples;
    results(s_idx).tof_us = tof_us;

    % Recalculate velocity using TOF
    thickness_m = results(s_idx).thickness_mm * 1e-3;
    tof_s = tof_samples / (SamplingRate_Mhz * 1e6);
    velocity_tof = (2 * thickness_m) / tof_s;
    results(s_idx).velocity_tof_ms = velocity_tof;

    fprintf('Sample %s:\n', results(s_idx).name);
    fprintf('  Front wall: %.1f samples (%.2f μs)\n', front_wall_estimate, front_wall_estimate/SamplingRate_Mhz);
    fprintf('  Back wall: %.1f samples (%.2f μs)\n', back_wall_time, back_wall_time/SamplingRate_Mhz);
    fprintf('  TOF: %.1f samples (%.2f μs)\n', tof_samples, tof_us);
    fprintf('  Velocity (TOF-based): %.0f m/s\n\n', velocity_tof);
end

%% Velocity Range Intersection Analysis

% Extract TOF-based velocities
tof_velocities = [results(valid_samples).velocity_tof_ms];
tof_mean = mean(tof_velocities);
tof_std = std(tof_velocities);

fprintf('=== TOF-BASED VELOCITY ANALYSIS ===\n');
fprintf('TOF-based velocities: ');
for i = 1:length(tof_velocities)
    fprintf('%.0f ', tof_velocities(i));
end
fprintf('m/s\n');
fprintf('Mean TOF velocity: %.0f ± %.0f m/s\n', tof_mean, tof_std);

% Calculate velocity ranges for each sample (with uncertainty)
velocity_ranges = zeros(length(valid_samples), 2);
for i = 1:length(valid_samples)
    s_idx = valid_samples(i);
    vel = results(s_idx).velocity_tof_ms;
    uncertainty = 0.1 * vel; % 10% uncertainty estimate
    velocity_ranges(i, :) = [vel - uncertainty, vel + uncertainty];

    fprintf('Sample %s velocity range: %.0f - %.0f m/s\n', ...
        results(s_idx).name, velocity_ranges(i, 1), velocity_ranges(i, 2));
end

% Find overlapping velocity range
v_lower_bound = max(velocity_ranges(:, 1));
v_upper_bound = min(velocity_ranges(:, 2));

% Apply material property constraints
v_lower_bound = max(v_lower_bound, Velocity_Range(1));
v_upper_bound = min(v_upper_bound, Velocity_Range(2));

if v_lower_bound <= v_upper_bound
    v_best_intersection = (v_lower_bound + v_upper_bound) / 2;
    fprintf('\n=== VELOCITY INTERSECTION ANALYSIS ===\n');
    fprintf('Overlapping velocity range: %.0f - %.0f m/s\n', v_lower_bound, v_upper_bound);
    fprintf('Best intersection velocity: %.0f m/s\n', v_best_intersection);
    fprintf('Range width: %.0f m/s (%.1f%% of mean)\n', ...
        v_upper_bound - v_lower_bound, (v_upper_bound - v_lower_bound)/v_best_intersection*100);
else
    fprintf('\n=== VELOCITY INTERSECTION ANALYSIS ===\n');
    fprintf('WARNING: No overlapping velocity range found!\n');
    fprintf('Samples may have inconsistent properties or measurement errors.\n');
    v_best_intersection = tof_mean; % Fallback to mean
end

%% Enhanced Visualization

% Create velocity vs TOF analysis plot (inspired by reference script)
figure('Position', [50, 50, 1200, 800], 'Name', 'Velocity vs TOF Analysis');

% Create TOF range for plotting
all_tof_samples = [results(valid_samples).tof_samples];
tof_range_samples = linspace(min(all_tof_samples)*0.8, max(all_tof_samples)*1.2, 100);
tof_range_us = tof_range_samples / SamplingRate_Mhz;

% Plot velocity vs TOF curves for each sample
subplot(2, 2, 1);
% Define safe colors to avoid RGB errors
safe_colors = [0.0, 0.4, 0.8; 0.8, 0.2, 0.2; 0.2, 0.7, 0.2; 0.8, 0.6, 0.0; 0.6, 0.2, 0.8];
colors = safe_colors(1:min(length(valid_samples), size(safe_colors, 1)), :);
for i = 1:length(valid_samples)
    s_idx = valid_samples(i);
    thickness_m = results(s_idx).thickness_mm * 1e-3;

    % Calculate velocity curve: v = 2*h/TOF
    tof_s = tof_range_samples / (SamplingRate_Mhz * 1e6);
    velocity_curve = (2 * thickness_m) ./ tof_s;

    plot(tof_range_us, velocity_curve, 'Color', colors(i,:), 'LineWidth', 2, ...
         'DisplayName', results(s_idx).name);
    hold on;

    % Mark actual measurement point
    actual_tof_us = results(s_idx).tof_us;
    actual_velocity = results(s_idx).velocity_tof_ms;
    plot(actual_tof_us, actual_velocity, 'o', 'Color', colors(i,:), ...
         'MarkerSize', 10, 'MarkerFaceColor', colors(i,:), 'HandleVisibility', 'off');
end

% Add velocity bounds
yline(Velocity_Range(1), 'k--', 'LineWidth', 1, 'DisplayName', 'Material Bounds');
yline(Velocity_Range(2), 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

if exist('v_best_intersection', 'var')
    yline(v_best_intersection, 'g-', 'LineWidth', 2, 'DisplayName', 'Best Intersection');
end

xlabel('Time-of-Flight (μs)');
ylabel('Velocity (m/s)');
title('Velocity vs TOF Analysis');
legend('Location', 'best');
grid on;
xlim([min(tof_range_us), max(tof_range_us)]);
ylim([Velocity_Range(1)*0.8, Velocity_Range(2)*1.2]);

% TOF comparison
subplot(2, 2, 2);
tof_values = [results(valid_samples).tof_us];
sample_names = {results(valid_samples).name};
bar(1:length(valid_samples), tof_values, 'FaceColor', [0.3, 0.7, 0.9]);
xlabel('Sample');
ylabel('Time-of-Flight (μs)');
title('TOF Comparison');
set(gca, 'XTickLabel', sample_names, 'XTickLabelRotation', 45);
grid on;

% Velocity comparison (both methods)
subplot(2, 2, 3);
x_pos = 1:length(valid_samples);
width = 0.35;

original_velocities = [results(valid_samples).calculated_velocity_ms];
tof_velocities = [results(valid_samples).velocity_tof_ms];

bar(x_pos - width/2, original_velocities, width, 'FaceColor', [0.8, 0.4, 0.4], 'DisplayName', 'Back Wall Only');
hold on;
bar(x_pos + width/2, tof_velocities, width, 'FaceColor', [0.4, 0.8, 0.4], 'DisplayName', 'TOF-based');

xlabel('Sample');
ylabel('Velocity (m/s)');
title('Velocity Calculation Comparison');
set(gca, 'XTickLabel', sample_names, 'XTickLabelRotation', 45);
legend('Location', 'best');
grid on;

% Velocity range intersection plot
subplot(2, 2, 4);
for i = 1:length(valid_samples)
    y_pos = i;
    vel_range = velocity_ranges(i, :);

    % Plot range as horizontal line
    plot(vel_range, [y_pos, y_pos], 'o-', 'Color', colors(i,:), 'LineWidth', 3, ...
         'MarkerSize', 8, 'DisplayName', sample_names{i});
    hold on;

    % Plot center point
    center_vel = mean(vel_range);
    plot(center_vel, y_pos, 's', 'Color', colors(i,:), 'MarkerSize', 10, ...
         'MarkerFaceColor', colors(i,:), 'HandleVisibility', 'off');
end

% Highlight intersection range
if exist('v_best_intersection', 'var') && v_lower_bound <= v_upper_bound
    fill([v_lower_bound, v_upper_bound, v_upper_bound, v_lower_bound], ...
         [0.5, 0.5, length(valid_samples)+0.5, length(valid_samples)+0.5], ...
         'g', 'FaceAlpha', 0.3, 'EdgeColor', 'g', 'LineWidth', 2, ...
         'DisplayName', 'Intersection Range');
end

xlabel('Velocity (m/s)');
ylabel('Sample');
title('Velocity Range Intersection');
set(gca, 'YTick', 1:length(valid_samples), 'YTickLabel', sample_names);
legend('Location', 'best');
grid on;
xlim([Velocity_Range(1)*0.9, Velocity_Range(2)*1.1]);
ylim([0.5, length(valid_samples)+0.5]);

sgtitle('Enhanced CFPP Velocity Analysis', 'FontSize', 16, 'FontWeight', 'bold');

% Create comprehensive comparison figure
figure('Position', [100, 100, 1400, 1000], 'Name', 'CFPP Wave Velocity Analysis');

% Define safe colors to avoid RGB errors
safe_colors = [0.0, 0.4, 0.8; 0.8, 0.2, 0.2; 0.2, 0.7, 0.2; 0.8, 0.6, 0.0; 0.6, 0.2, 0.8];

% Subplot 1: Individual sample velocity calculations
subplot(2, 3, 1);
bar(valid_samples, velocities, 'FaceColor', [0.2, 0.6, 0.8]);
hold on;
yline(mean_velocity, 'r--', 'LineWidth', 2, 'Label', sprintf('Mean: %.0f m/s', mean_velocity));
yline(mean_velocity + std_velocity, 'r:', 'LineWidth', 1, 'Label', sprintf('+1σ: %.0f m/s', mean_velocity + std_velocity));
yline(mean_velocity - std_velocity, 'r:', 'LineWidth', 1, 'Label', sprintf('-1σ: %.0f m/s', mean_velocity - std_velocity));
xlabel('Sample Number');
ylabel('Calculated Velocity (m/s)');
title('Individual Sample Velocities');
grid on;
legend('Location', 'best');

% Subplot 2: Velocity optimization curve
subplot(2, 3, 2);
plot(test_velocities, error_scores, 'b-', 'LineWidth', 2);
hold on;
plot(optimal_velocity, min_error, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('Test Velocity (m/s)');
ylabel('RMS Error');
title('Velocity Optimization');
grid on;
legend(sprintf('Optimal: %.0f m/s', optimal_velocity), 'Location', 'best');

% Subplot 3: Layer pattern comparison
subplot(2, 3, 3);
colors = safe_colors(1:min(length(valid_samples), size(safe_colors, 1)), :);
for i = 1:length(valid_samples)
    s_idx = valid_samples(i);
    layer_centers = results(s_idx).layer_centers;
    layer_values = results(s_idx).layer_values;

    % Plot as step function
    for j = 1:length(layer_centers)
        if j == 1
            plot([0, layer_centers(j)], [layer_values(j), layer_values(j)], ...
                 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', results(s_idx).name);
        else
            plot([layer_centers(j-1), layer_centers(j)], [layer_values(j), layer_values(j)], ...
                 'Color', colors(i,:), 'LineWidth', 2, 'HandleVisibility', 'off');
        end
        hold on;
    end
end
xlabel('Time (samples)');
ylabel('Layer Value');
title('Layer Patterns Comparison');
legend('Location', 'best');
grid on;
ylim([-1.5, 1.5]);

% Subplot 4: Thickness vs Back Wall Time
subplot(2, 3, 4);
thicknesses = [results(valid_samples).thickness_mm];
back_wall_times = [results(valid_samples).back_wall_time_us];
scatter(back_wall_times, thicknesses, 100, velocities, 'filled');
colorbar;
colormap(jet);
xlabel('Back Wall Time (samples)');
ylabel('Thickness (mm)');
title('Thickness vs Time (colored by velocity)');
grid on;

% Add theoretical line for optimal velocity
time_theory = linspace(min(back_wall_times)*0.8, max(back_wall_times)*1.2, 100);
thickness_theory = (optimal_velocity * time_theory / (SamplingRate_Mhz * 1e6)) / 2 * 1e3; % Convert to mm
hold on;
plot(time_theory, thickness_theory, 'k--', 'LineWidth', 2, 'DisplayName', sprintf('Theory (%.0f m/s)', optimal_velocity));
legend('Location', 'best');

% Subplot 5: Velocity histogram
subplot(2, 3, 5);
histogram(velocities, 'FaceColor', [0.6, 0.8, 0.6], 'EdgeColor', 'black');
hold on;
xline(mean_velocity, 'r--', 'LineWidth', 2, 'Label', sprintf('Mean: %.0f m/s', mean_velocity));
xline(optimal_velocity, 'g--', 'LineWidth', 2, 'Label', sprintf('Optimal: %.0f m/s', optimal_velocity));
xlabel('Velocity (m/s)');
ylabel('Count');
title('Velocity Distribution');
legend('Location', 'best');
grid on;

% Subplot 6: Error analysis
subplot(2, 3, 6);
sample_names = {results(valid_samples).name};
actual_times = [results(valid_samples).back_wall_time_us];
predicted_times = (2 * [results(valid_samples).thickness_mm] * 1e-3) ./ optimal_velocity * (SamplingRate_Mhz * 1e6); % Convert to samples
errors = ((predicted_times - actual_times) ./ actual_times) * 100; % Percentage error

bar(1:length(valid_samples), errors, 'FaceColor', [0.8, 0.4, 0.4]);
xlabel('Sample');
ylabel('Prediction Error (%)');
title('Time Prediction Errors (Optimal Velocity)');
set(gca, 'XTickLabel', sample_names, 'XTickLabelRotation', 45);
grid on;
yline(0, 'k--', 'LineWidth', 1);

% Adjust layout
sgtitle(sprintf('CFPP Wave Velocity Analysis - Optimal Velocity: %.0f m/s', optimal_velocity), 'FontSize', 16, 'FontWeight', 'bold');

%% Additional Analysis Functions

% Function to analyze layer periodicity
function analyze_layer_periodicity(samples, results, valid_samples)
    figure('Position', [200, 200, 1200, 800], 'Name', 'Layer Periodicity Analysis');

    for i = 1:length(valid_samples)
        s_idx = valid_samples(i);
        subplot(2, 2, i);

        layer_centers = results(s_idx).layer_centers;
        layer_values = results(s_idx).layer_values;

        % Calculate layer transitions (from -1 to 1 or 1 to -1)
        transitions = find(diff(layer_values) ~= 0);
        if ~isempty(transitions)
            transition_times = layer_centers(transitions);
            layer_periods = diff(transition_times);

            % Create a clearer plot showing layer spacing consistency
            transition_numbers = 1:length(layer_periods);

            plot(transition_numbers, layer_periods, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.2, 0.6, 0.8]);
            hold on;

            % Add mean line for reference
            mean_period = mean(layer_periods);
            yline(mean_period, '--', 'Color', [0.8, 0.2, 0.2], 'LineWidth', 2, ...
                  'Label', sprintf('Mean: %.1f samples', mean_period));

            % Add ±1 std deviation bands
            std_period = std(layer_periods);
            fill([transition_numbers, fliplr(transition_numbers)], ...
                 [ones(size(transition_numbers))*(mean_period + std_period), ...
                  fliplr(ones(size(transition_numbers))*(mean_period - std_period))], ...
                 [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
                 'DisplayName', '±1σ band');

            xlabel('Layer Transition Number');
            ylabel('Distance Between Transitions (samples)');
            title(sprintf('%s - Layer Spacing Consistency', strrep(results(s_idx).name, '_', '\_')));
            grid on;

            % Add detailed statistics
            cv = std_period / mean_period;
            text(0.05, 0.95, sprintf('Mean: %.1f samples\nStd: %.1f samples\nCV: %.1f%%', ...
                 mean_period, std_period, cv*100), ...
                 'Units', 'normalized', 'VerticalAlignment', 'top', 'BackgroundColor', 'white', ...
                 'FontSize', 10);

            % Add interpretation
            if cv < 0.2
                quality_text = 'Good uniformity';
                quality_color = [0, 0.7, 0];
            elseif cv < 0.4
                quality_text = 'Moderate uniformity';
                quality_color = [0.8, 0.6, 0];
            else
                quality_text = 'Poor uniformity';
                quality_color = [0.8, 0, 0];
            end

            text(0.05, 0.75, quality_text, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
                 'BackgroundColor', 'white', 'Color', quality_color, 'FontWeight', 'bold');
        else
            text(0.5, 0.5, 'No layer transitions detected', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 12, 'Color', [0.6, 0.6, 0.6]);
            title(sprintf('%s - No Transitions', strrep(results(s_idx).name, '_', '\_')));
        end
    end

    sgtitle('Layer Spacing Uniformity Analysis', 'FontSize', 14, 'FontWeight', 'bold');
end

% Function to create detailed layer visualization
function create_layer_visualization(samples, results, valid_samples, optimal_velocity, sampling_rate_mhz)
    figure('Position', [300, 300, 1400, 600], 'Name', 'Detailed Layer Visualization');

    for i = 1:length(valid_samples)
        s_idx = valid_samples(i);
        subplot(1, length(valid_samples), i);

        layers = samples(s_idx).layers;
        values = samples(s_idx).values;

        % Create a detailed visualization of each layer
        y_pos = 1;
        colors = [0.8, 0.2, 0.2; 0.2, 0.2, 0.8]; % Red for -1, Blue for 1

        for j = 1:size(layers, 1)
            time_start = layers(j, 1);
            time_end = layers(j, 2);
            value = values(j);

            color_idx = (value == 1) + 1; % 1 for -1 values, 2 for 1 values

            rectangle('Position', [time_start, y_pos-0.4, time_end-time_start, 0.8], ...
                     'FaceColor', colors(color_idx, :), 'EdgeColor', 'black');
            hold on;

            % Add value text
            text(mean([time_start, time_end]), y_pos, num2str(value), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontWeight', 'bold', 'Color', 'white');
        end

        % Calculate and show theoretical layer positions based on optimal velocity
        thickness = results(s_idx).thickness_mm;
        expected_back_wall = (2 * thickness * 1e-3) / optimal_velocity * (sampling_rate_mhz * 1e6); % Convert to samples

        xline(expected_back_wall, 'g--', 'LineWidth', 3, 'Label', 'Theoretical Back Wall');

        xlabel('Time (samples)');
        ylabel('Layer');
        title(sprintf('%s\n%.3f mm thick', strrep(results(s_idx).name, '_', '\_'), thickness));
        ylim([0.5, 1.5]);
        xlim([0, max(layers(:)) * 1.1]);
        grid on;
        legend('Location', 'best');
    end

    sgtitle(sprintf('Layer Structure Visualization (Optimal Velocity: %.0f m/s)', optimal_velocity), ...
            'FontSize', 14, 'FontWeight', 'bold');
end

% Call additional analysis functions
analyze_layer_periodicity(samples, results, valid_samples);
create_layer_visualization(samples, results, valid_samples, optimal_velocity, SamplingRate_Mhz);

%% Data Export and Summary Report

% Create summary table
summary_table = table();
for i = 1:length(valid_samples)
    s_idx = valid_samples(i);
    summary_table.Sample{i} = results(s_idx).name;
    summary_table.Thickness_mm(i) = results(s_idx).thickness_mm;
    summary_table.BackWallTime_us(i) = results(s_idx).back_wall_time_us;
    summary_table.CalculatedVelocity_ms(i) = results(s_idx).calculated_velocity_ms;
    summary_table.NumLayers(i) = length(results(s_idx).layer_centers);

    % Calculate prediction error using optimal velocity
    predicted_time = (2 * results(s_idx).thickness_mm * 1e-3) / optimal_velocity * (SamplingRate_Mhz * 1e6); % Convert to samples
    summary_table.PredictionError_percent(i) = ((predicted_time - results(s_idx).back_wall_time_us) / results(s_idx).back_wall_time_us) * 100;
end

% Display summary table
fprintf('\n=== DETAILED SAMPLE ANALYSIS ===\n');
disp(summary_table);

% Save results to file
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
filename = sprintf('CFPP_Velocity_Analysis_%s.mat', timestamp);
save(filename, 'samples', 'results', 'summary_table', 'optimal_velocity', 'mean_velocity', 'std_velocity');
fprintf('\nResults saved to: %s\n', filename);

% Export summary to CSV
csv_filename = sprintf('CFPP_Velocity_Summary_%s.csv', timestamp);
writetable(summary_table, csv_filename);
fprintf('Summary table exported to: %s\n', csv_filename);

% Generate final report
fprintf('\n=== FINAL RECOMMENDATIONS ===\n');
fprintf('Based on the analysis of %d CFPP samples:\n', length(valid_samples));

fprintf('\n--- BACK WALL ONLY METHOD ---\n');
fprintf('• Recommended wave velocity: %.0f m/s\n', optimal_velocity);
fprintf('• Velocity range observed: %.0f - %.0f m/s\n', min(velocities), max(velocities));
fprintf('• Standard deviation: %.0f m/s (%.1f%% variation)\n', std_velocity, (std_velocity/mean_velocity)*100);

fprintf('\n--- TIME-OF-FLIGHT METHOD ---\n');
fprintf('• TOF-based mean velocity: %.0f m/s\n', tof_mean);
fprintf('• TOF velocity range: %.0f - %.0f m/s\n', min(tof_velocities), max(tof_velocities));
fprintf('• TOF standard deviation: %.0f m/s (%.1f%% variation)\n', tof_std, (tof_std/tof_mean)*100);

if exist('v_best_intersection', 'var') && v_lower_bound <= v_upper_bound
    fprintf('\n--- INTERSECTION ANALYSIS ---\n');
    fprintf('• Best intersection velocity: %.0f m/s\n', v_best_intersection);
    fprintf('• Intersection range: %.0f - %.0f m/s\n', v_lower_bound, v_upper_bound);
    fprintf('• Range confidence: %.1f%% (narrower is better)\n', (v_upper_bound - v_lower_bound)/v_best_intersection*100);

    recommended_velocity = v_best_intersection;
    fprintf('\n*** RECOMMENDED VELOCITY: %.0f m/s (intersection method) ***\n', recommended_velocity);
else
    recommended_velocity = tof_mean;
    fprintf('\n*** RECOMMENDED VELOCITY: %.0f m/s (TOF mean - no intersection) ***\n', recommended_velocity);
end

fprintf('\n--- PREDICTION ACCURACY ---\n');
fprintf('• Maximum prediction error: %.1f%%\n', max(abs(summary_table.PredictionError_percent)));
fprintf('• Average prediction error: %.1f%%\n', mean(abs(summary_table.PredictionError_percent)));

% Assessment based on TOF method (more accurate)
if tof_std/tof_mean < 0.1
    fprintf('• Assessment: GOOD - Low velocity variation between samples\n');
elseif tof_std/tof_mean < 0.2
    fprintf('• Assessment: MODERATE - Some velocity variation between samples\n');
else
    fprintf('• Assessment: HIGH - Significant velocity variation between samples\n');
end

fprintf('\n--- METHOD COMPARISON ---\n');
fprintf('• Back wall method tends to be less accurate (ignores front wall)\n');
fprintf('• TOF method is more physically accurate (accounts for full wave travel)\n');
fprintf('• Intersection method provides best estimate when samples agree\n');

if Enable_Data_Filtering
    fprintf('\n--- DATA FILTERING IMPACT ---\n');
    fprintf('• Data filtering was ENABLED\n');
    fprintf('• Filtered out samples with known defects or quality issues\n');
    fprintf('• Remaining samples: %d out of %d total\n', length(valid_samples), sum(~cellfun(@isempty, {samples.name})));

    % Show which samples were filtered
    filtered_names = {};
    for i = 1:length(quality_report)
        if ~isempty(quality_report(i).name) && quality_report(i).recommended_filter
            filtered_names{end+1} = quality_report(i).name;
        end
    end
    if ~isempty(filtered_names)
        fprintf('• Filtered samples: %s\n', strjoin(filtered_names, ', '));
    end

    fprintf('• This may improve velocity consistency by removing defective data\n');
    fprintf('• To include all samples, set Enable_Data_Filtering = false\n');
else
    fprintf('\n--- DATA FILTERING IMPACT ---\n');
    fprintf('• Data filtering was DISABLED - using all available samples\n');
    fprintf('• Results may include effects from samples with known defects\n');
    fprintf('• To filter problematic samples, set Enable_Data_Filtering = true\n');
end

fprintf('\nNote: All low values have been set to -1 and high values to 1 as requested.\n');
fprintf('Time values are interpreted as sample points at %.0f MHz sampling rate.\n', SamplingRate_Mhz);

%% Additional Utility Functions

% Function to predict layer structure for new samples
function predict_layer_structure(thickness_mm, velocity_ms, output_freq_mhz, sampling_rate_mhz)
    fprintf('\n=== LAYER STRUCTURE PREDICTION ===\n');
    fprintf('For a sample with thickness %.3f mm:\n', thickness_mm);

    % Calculate expected back wall time
    back_wall_time_samples = (2 * thickness_mm * 1e-3) / velocity_ms * (sampling_rate_mhz * 1e6); % Convert to samples
    back_wall_time_us = back_wall_time_samples / sampling_rate_mhz; % Convert to μs for display
    fprintf('Expected back wall time: %.1f samples (%.2f μs)\n', back_wall_time_samples, back_wall_time_us);

    % Estimate number of observable layers based on frequency resolution
    wavelength_mm = velocity_ms / (output_freq_mhz * 1e6) * 1e3;
    min_resolvable_thickness = wavelength_mm / 4; % Quarter wavelength resolution
    max_layers = floor(thickness_mm / min_resolvable_thickness);

    fprintf('Wavelength in material: %.3f mm\n', wavelength_mm);
    fprintf('Minimum resolvable thickness: %.3f mm\n', min_resolvable_thickness);
    fprintf('Maximum resolvable layers: %d\n', max_layers);
end

% Example prediction for a new sample
fprintf('\n=== EXAMPLE: New Sample Prediction ===\n');
predict_layer_structure(2.0, optimal_velocity, OutputFreq_Mhz, SamplingRate_Mhz);

fprintf('\n=== ANALYSIS COMPLETE ===\n');
fprintf('Check the generated figures for visual comparisons.\n');