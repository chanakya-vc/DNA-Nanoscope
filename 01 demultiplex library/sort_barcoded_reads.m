% Abstract
% 1. Reads are in the form of fastq files and are situated in the
% fastq_directory. 
% 2. fastq files have filenames of the form fastq_prefix_0, fastq_prefix_1
% .... fastq_prefix_numfiles
% 3. Output structure: Directories BC1, BC2, ... , BC48 in
% fastq_directory/demux. BC1 contains files fastq_prefix_0, fastq_prefix_1
% .... fastq_prefix_numfiles with BC1 reads only. Ditto for the other
% barcodes
%

function sort_barcoded_reads(fastq_directory, barcodes)
total_time = tic; % start timing the total time to process all files

barcodes = char(barcodes);
num_barcodes = size(barcodes,1);
barcodes = nt2int(barcodes); % convert barcodes from char to uint
file_list = dir(fastq_directory); % Get a list of fastq filenames
numfiles = size(file_list,1)-2; % number of fastq files in the path

% creat directories to store the demultiplexed reads
for i=1:(num_barcodes/2)
%mkdir(strcat(fastq_directory,'demux/'),strcat('BC',int2str(i)));% Mac
mkdir(strcat(fastq_directory,'demux\'),strcat('BC',int2str(i))); % PC
end

% for each read file, identify the correct barcode of the read and write it to
% the appropriate location
parfor k=1:numfiles
    for_time = tic; % start timing current loop
    fastq_filename = file_list(k+2).name;
    split_barcodes(fastq_directory,fastq_filename,barcodes); % append barcode id to
    % print progress
    elapsed_time = toc(for_time);   
    %file_id = k;
    fprintf('Processed file %s in %.0f seconds\n',fastq_filename,elapsed_time);
end

% stop timing and print total elapsed time
total_elapsed_time = toc(total_time);
fprintf('Processed all files in %.0f seconds\n',total_elapsed_time);

end