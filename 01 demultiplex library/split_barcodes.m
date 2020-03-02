% Called by sort_barcoded_reads

function [fastqstruct, max_score, best_alignment] = split_barcodes(fastq_directory, fastq_file_name, barcodes)

% read a fastq file
fastqfile_fullpath = strcat(fastq_directory,fastq_file_name); % location of the read file
fastqstruct = fastqread(fastqfile_fullpath); % import all reads

% useful short hand variables
n = size(barcodes,1); % number of barcodes
m = length(fastqstruct); % number of reads

% convert nucleotides to int representation to save time later
fastqstruct_nt = fastqstruct; % preallocate memory
for i=1:m
fastqstruct_nt(i).Sequence = nt2int(fastqstruct(i).Sequence); % convert reads to int
end

% max_score is the best alignment score for each sequence
max_score = ones(m,1) - 1; % preallocate
best_alignment = cell(m,1); % will store the alignment found
barcode_id = zeros(m,1); % will store the id of the barcode for each sequence

% align barcodes to sequences, retain best score and write to fastq
for i=1:m % for each read in file
    
    for j=1:n % for each barcode
        [max_score_temp, alignment_temp] = nwalign(fastqstruct_nt(i),barcodes(j,:),'Alphabet','NT','Glocal','true'); % align current barcode
        
        if max_score_temp >= max_score(i) % if this is the best alignment found so far
            max_score(i) = max_score_temp; % update max score
            best_alignment{i} = alignment_temp; % update corresponding alignment
            barcode_id(i) = j; % update barcode
        end
    end
    
    fastqstruct(i).Header = strcat(fastqstruct(i).Header,' barcode=',int2str(barcode_id(i)), ' alignment_score=',num2str(max_score(i))); % save barcode id to header

    if barcode_id(i) > n/2 % bc_(n/2 + 1) to bc_n are reverse complements of bc1 to bc_n/2
        barcode_id(i) = barcode_id(i) - n/2; % store the read with bc_(i-n/2)
    end
    
end

% split fastqstruct into different fastqstruct based on barcodes and write
% it to the appropriate directory.
for j=1:(n/2) % for each barcode
    fastqstruct_filtered = fastqstruct(barcode_id == j); % filter only barcode_j reads
    
    % write to appropriate directory
    %warnState = warning; %Save the current warning state
    %warning('off','Bioinfo:fastqwrite:AppendToFile');
    %full_filename = strcat(fastq_directory,'demux/BC',int2str(barcode_id),'/',fastq_prefix,int2str(k),'.fastq'); % mac
    full_filename = strcat(fastq_directory,'demux\BC',int2str(j),'\',fastq_file_name); % PC
    fastqwrite(full_filename,fastqstruct_filtered);
    %warning(warnState) %Reset warning state to previous settings
    
end

end
