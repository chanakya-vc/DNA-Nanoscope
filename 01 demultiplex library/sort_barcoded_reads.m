% Usage: 
% sort_barcoded_reads(fastq_dir,ONT_barcodes)
% where:
% fastq_dir is the path of the fastq sequence files and
% barcodes (see file ONT_library_barcodes.mat) specifies the sequencing barcodes (from Oxford Nanotech i.e. ONT) and their reverse complements.
% The output is written into fastq files with reads that are sorted into sub-directories corresponding to the identity of the sequencing barcodes.

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
    fprintf('Processed file %s in %.0f seconds\n',fastq_filename,elapsed_time);
end

% stop timing and print total elapsed time
total_elapsed_time = toc(total_time);
fprintf('Processed all files in %.0f seconds\n',total_elapsed_time);

end