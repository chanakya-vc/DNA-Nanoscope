% Usage:
% pairwise_record_list = extract_pairwise_record_lengths(target_barcodes,path,library_size,color_length)
% where:
% target_barcodes specifies the staple barcode sequences,
% path specifies the path of the fastq reads,
% library_size specifies the number of sequencing libraries that were
% combined for a single run
% color_length specifies the length of the auxiliary tag sequence
% pairwise_record_list is the output, a matrix of size (n, n, 2001) where n is the number of target points and cell (i, j, k) holds the number of distance records of length k bases (only counting the repeat region) between points i and j. All distance records of length > 2000 are stored in the slice (:, :, 2001).
% pairwise_record_list variables from different sequencing runs of the same experiment were combined by simply adding them.

function [pairwise_record_list] = extract_pairwise_record_lengths(staple_barcodes,path, library_size, color_length)

% staple_barcodes is a string column vectors containing the staple barcode sequences
% We require that each barcode be the same length (staple_barcode_length)

% preset
% LongPrimers were used in calibration experiments
% LongPrimer1 = fastnt2int('TCGTGCGAGTATAGAAAGTGAGGGATTAATGG'); %
% LongPrimer2 = fastnt2int('TCTACCCCATGAAGAGTAAATAGGTTGTGGGA');
% RCLongPrimer1 = fastnt2int('CCATTAATCCCTCACTTTCTATACTCGCACGA');
% RCLongPrimer2 = fastnt2int('TCCCACAACCTATTTACTCTTCATGGGGTAGA');

% These primers are used in pattern reconstruction experiments
Primer1 = 'GTCTTGAGCAAATAGCAGGTGACA';
Primer2 = 'TCCATCTTGTCTGTTAGCAAGCTG';
RCPrimer1 = seqrcomplement(Primer1);
RCPrimer2 = seqrcomplement(Primer2);
primer_length = length(Primer1);
staple_barcode_length = strlength(staple_barcodes(1));
threshold_primer = floor(0.5*5*primer_length); % best possible score = 5 * primer_length
threshold_staple = floor(0.5*5*staple_barcode_length); % best possible score = 5 * staple_length
max_reclen = 2001;
% library_size = 1; % set to number of samples in the library
% color_length = 0;

total_time = tic; % start timing the total time to process all files

% convert staple sequences to int
num_of_points = size(staple_barcodes,1);
staple_barcodes_int{num_of_points,2} = []; % preallocate cell array to hold int-converted sequences
for n=1:num_of_points
    staple_barcodes_int{n,1} = fastnt2int(char(staple_barcodes(n)));
    staple_barcodes_int{n,2} = fastnt2int(seqrcomplement(char(staple_barcodes(n)))); % reverse complement
end

file_list = dir(path); % list of all file names, sizes etc
num_files = size(file_list,1) - 2; % number of fastq files in the path
% preallocate
%pairwise_record_list = zeros(num_of_points,num_of_points,max_reclen,num_files,'uint16'); %% SLICING VARIABLE USES TOO MUCH MEMORY
pairwise_record_list = zeros(num_of_points,num_of_points,max_reclen,'uint32');

parfor i=1:num_files % for each file
    for_time = tic;
    % read sequencing file into memory
    full_filename = strcat(path,file_list(i+2).name);
    [header, sequence] = fastqread(full_filename); % returns cell row array
    
    local_record_list = zeros(num_of_points,num_of_points,max_reclen,'uint32'); % instantiate a local copy
    num_reads = size(sequence,2);
    
    for j=1:num_reads
        % extract barcode_id
      
        index_start = strfind(header{j},'barcode='); % !!!returns cell array!!!
        index = index_start + length('barcode=');
        barcode_id = str2double(header{j}(index:index+1));
        seq = fastnt2int(sequence{j});
        
        % find locations of the staple ssequences: step 1
        % !!! Need to fix barcode ids so that they are consistent across experiments !!!
        if barcode_id >= 1 && barcode_id <= library_size % barcode_id is in (1,k), i.e. top_strand was sequenced
            strand_sequenced = 1; % i.e. top strand 
            % align with LongPrimer1 --- RevLongPrimer2
            [zalign_score1, alignment1, start] = fastswalign(seq,Primer1,'Alphabet','NT');
            [zalign_score2, alignment2, finish] = fastswalign(seq,RCPrimer2,'Alphabet','NT');

        elseif barcode_id >= library_size+1 && barcode_id <= 3*library_size % bottom_strand was sequenced
            strand_sequenced = 2; % i.e. bottom strand was sequenced
            % align with LongPrimer2 --- RevLongPrimer1
            [zalign_score1, alignment1, start] = fastswalign(seq,Primer2,'Alphabet','NT');
            [zalign_score2, alignment2, finish] = fastswalign(seq,RCPrimer1,'Alphabet','NT');
        
        elseif barcode_id >= 3*library_size + 1 && barcode_id <= 4*library_size % top_strand was sequenced
            strand_sequenced = 1; % i.e. top strand was sequenced
            % align with LongPrimer1 --- RevLongPrimer2
            [zalign_score1, alignment1, start] = fastswalign(seq,Primer1,'Alphabet','NT');
            [zalign_score2, alignment2, finish] = fastswalign(seq,RCPrimer2,'Alphabet','NT');
        
        end
        % showalignment(alignment1,'StartPointers',start); showalignment(alignment2,'StartPointers',finish);
        %seq_nt = int2nt(seq);
        %full_alignment = [seq_nt(1:start(1)-1) alignment1(1,:) seq_nt(start(1)+size(alignment1,2):
        
%         showalignment(alignment1); showalignment(alignment2);
        
        % find locations of the staple sequences: step 2
        if zalign_score1 > threshold_primer && zalign_score2 > threshold_primer && finish(1) > (start(1) + primer_length - start(2) + 2*color_length + 2*staple_barcode_length) % if good alignments found
            staple1_start_pos = start(1) + primer_length - (start(2) - 1) + color_length;
            staple2_start_pos = finish(1) - color_length - staple_barcode_length;
            
            epsilon = 5; % wiggle room to allow for a better alignment
            staple1_supersequence = seq(staple1_start_pos-epsilon:staple1_start_pos+staple_barcode_length+epsilon-1); % extract sequence for alignment
            staple2_supersequence = seq(staple2_start_pos-epsilon:staple2_start_pos+staple_barcode_length+epsilon-1);
            
            % find the best 'staple' match to the sequence read
            zstaple1_align_score = -1; % store alignment scores here
            zstaple2_align_score = -1;
            staple1_offset = [1; 1]; % store offsets here
            staple2_offset = [1; 1];
            staple1_id = 0; % store ID (1 to num_points) of best aligned staples here
            staple2_id = 0;
            for k=1:num_of_points % for each staple
                [temp_align_score1, temp_alignment1, temp_start1] = fastswalign(staple1_supersequence,staple_barcodes_int{k,1},'Alphabet','NT'); % align against staple
                if temp_align_score1 > zstaple1_align_score % if better alignment found
                    zstaple1_align_score = temp_align_score1; % update score
                    staple1_offset = temp_start1; % update location
                    staple1_id = k; % update staple id
                    staple_alignment1 = temp_alignment1; % update alignment
                    if zstaple1_align_score >= 0.9*5*staple_barcode_length % if a really good alignment is found
                        break; % stop searching and exit the for loop
                    end
                end
            end
            
            for k=1:num_of_points
                [temp_align_score2, temp_alignment2, temp_start2] = fastswalign(staple2_supersequence,staple_barcodes_int{k,2},'Alphabet','NT');
                if temp_align_score2 > zstaple2_align_score % if better alignment found
                    zstaple2_align_score = temp_align_score2; % update score
                    staple2_offset = temp_start2; % update location
                    staple2_id = k; % update staple id
                    staple_alignment2 = temp_alignment2; % update alignment
                    if zstaple2_align_score >= 0.95*5*staple_barcode_length % if a really good alignment is found
                        break; % stop searching and exit the for loop
                    end
                end
            end
            
            %showalignment(staple_alignment1); showalignment(staple_alignment2);
            repeat_start = staple1_start_pos - epsilon + staple1_offset(1) + staple_barcode_length - staple1_offset(2);
            repeat_end = staple2_start_pos - epsilon + staple2_offset(1) - staple2_offset(2) - 1;
            repeat_length = repeat_end - repeat_start +  1; % extract length
                
            % if the alignment scores are good, extract the length of the record
            if zstaple1_align_score > threshold_staple && zstaple2_align_score > threshold_staple && repeat_length >= 4

                % check if we have the right repeat sequence
                record_repeat_sequence = seq(repeat_start:repeat_end);
                if strand_sequenced == 1 % top strand
                        ideal_repeat_seq = repmat('AAAT',1,floor(repeat_length/4));
                elseif strand_sequenced == 2 % bottom strand
                        ideal_repeat_seq = repmat('ATTT',1,floor(repeat_length/4));
                end
                zzrepeat_pdist = seqpdist({fastint2nt(record_repeat_sequence);ideal_repeat_seq},'Alphabet','NT','Method','p-distance');
                
                if zzrepeat_pdist < 0.25 % If sequences are less than 25% different, then proceed to record length
                    if staple2_id > staple1_id % first id should be bigger
                        temp_swap_variable = staple1_id;
                        staple1_id = staple2_id;
                        staple2_id = temp_swap_variable;
                    end
                
                if repeat_length > 2000
                    repeat_length = 2001; % catch all case for extra long records
                end
                
                local_record_list(staple1_id,staple2_id,repeat_length) = local_record_list(staple1_id,staple2_id,repeat_length) + 1;
                                
                end
                
            end
            
        end
        
        
    end
    %pairwise_record_list(:,:,:,i) = local_record_list ; % DEPRECATED. NOT USING SLICED VARS ANYMORE % update the sliced copy
    pairwise_record_list = pairwise_record_list + local_record_list; % add local computation to global reduction variable
    
    % print progress
    elapsed_time = toc(for_time);   
    fprintf('Processed file %s in %.0f seconds\n',file_list(i+2).name,elapsed_time);
    
end

% print progress
total_elapsed_time = toc(total_time);
fprintf('Processed all files in %.0f seconds\n',total_elapsed_time);
end
