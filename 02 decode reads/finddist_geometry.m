% Usage:
% [pairwise_distances, pairwise_peak_heights] = finddist_geometry(pairwise_records_list,calibration_fun)
% where:
% pairwise_records_list is the output of extract_pairwise_record_lengths,
% calibration_fun is a cfit object that holds the calibration function of the ruler, which maps bases to nanometers,
% pairwise_distances is one output, an array of measured distances for each pair of points, and
% pairwise_peak_heights is the peak height (in bases) corresponding to the distances measured. It is a measure of the confidence in the measurement.

function [pairwise_distances, pairwise_peak_heights] = finddist_geometry(pairwise_records_list,calibration_fun)

%useful variables
offset = 66; % includes primer regions and TYE665 dye
n = size(pairwise_records_list,1);

[pairwise_peaks, pairwise_peak_heights] = findpeaks_geometry(pairwise_records_list); % obtain peak locations in bp
pairwise_distances = tril(vec2mat(calibration_fun(pairwise_peaks+offset),n)'); % convert peaks from bases to nanometers
pairwise_distances(isnan(pairwise_distances)) = 0; % replace NaNs by 0s
end