% Called by finddist_geometry
% What it does: finds the peak of a distribution of record lengths and the corresponding
% peak height
% Input: pairwise record list is a 3D matrix containing lists of pairwise records
function [pairwise_record_peaks, pairwise_peak_heights] = findpeaks_geometry(pairwise_record_list)

% useful variables
n = size(pairwise_record_list,1);
m = size(pairwise_record_list,3); % max(domain) of record lengths
x = 1:m;
min_guess_threshold = 0;

pairwise_record_peaks = zeros(n); % preallocate
pairwise_peak_heights = zeros(n);

% trim to delete short record lengths which are typically background
pairwise_record_list(:,:,x <= 12) = 0;

for i=1:n
    for j=1:i
        y = pairwise_record_list(i,j,:);
        y = double(y(:));
        
        % smooth data
         span = 8;
         y = smooth(y,span);
        
        % extract location and height of major peak as follows
        guess_prominence = 2^7; % guess a high value for the peak prominence
        a_guess = []; % initialization
        while isempty(a_guess) % while no peak found
            guess_prominence = max(min_guess_threshold, floor(guess_prominence/2)); % reduce the guess value, since no peak has been found so far
            [a_guess, b_guess, c_guess] = findpeaks(y,x,'MinPeakProminence',guess_prominence); % find peaks and store height, location and width

            if size(a_guess,1) > 1 % if more than one peak found
                [a_guess, ix] = max(a_guess); % pick only the highest peak. This will exit while loop
                b_guess = b_guess(ix); % and the corresponding position
                c_guess = c_guess(ix); % and the corresponding peak width
            end
            
            if guess_prominence == min_guess_threshold && isempty(a_guess) % threshold is reached and no peaks found
                % stop searching
                a_guess = 0; % preset value: no peak found. This will exit the while loop
                b_guess = NaN; % preset value: no peak found
                c_guess = 0; % preset value: no peak found
            end
        end
        % store value
        pairwise_record_peaks(i,j) = b_guess;
        pairwise_peak_heights(i,j) = a_guess;
        
    end
end

end