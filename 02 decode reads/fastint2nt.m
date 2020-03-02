% Forked from MATLAB function int2nt
% Edited by nikhil314 on 2020.03.02
% Certain section deleted since they slow down the code and are not used in
% our specific application

function nt = fastint2nt(sequence,varargin)
%INT2NT converts a sequence of integers to a character array of nucleotides.
%
%   NT = INT2NT(SEQ) converts the sequence SEQ from integers into a sequence
%   of characters representing the nucleotides A,C,T,G. Unknown nucleotides
%   are represented by the symbol '*'.
%
%   N = INT2NT(SEQ,...,'ALPHABET',A) defines which nucleotide alphabet to
%   use. The default value is 'DNA' which uses the symbols A,C,T,G. If
%   ALPHABET is set to 'RNA', then A,C,U,G are used instead.
%
%   N = INT2NT(SEQ,...,'UNKNOWN',C) defines the symbol used to represent an
%   unknown nucleotide. The default value is '*'.
%
%   N = INT2NT(SEQ,...,'CASE',case) sets the output case of the nucleotide
%   string. Default is uppercase.
%
%   Example:
%
%       s = int2nt([1 2 4 3 2 4 1 3 2])
%
%   See also AA2INT, INT2AA, NT2INT.

%   Copyright 2002-2012 The MathWorks, Inc.


if isempty(sequence)
    nt = '';
    return
end

%%% DELETED SECTION FOR SPEED
% % If the input is a structure then extract the Sequence data.
% if isstruct(sequence)
%     sequence = bioinfoprivate.seqfromstruct(sequence);
% end
%%% DELETED SECTION FOR SPEED

%%% DELETED SECTION FOR SPEED
% origsize = size(sequence);
% sequence = sequence(:);
% sequence = sequence';
%%% DELETED SECTION FOR SPEED

% % % A C G T U R Y K M S  W  B  D  H  V N  -(gap)
% % % 1 2 3 4 4 5 6 7 8 9 10 11 12 13 14 15 16
% % map = '*ACGTRYKMSWBDHVN-';

%%% new code %%%%%
% A C G T
% 1 2 3 4
map = '*ACGT';
%%% new code %%%%%

if nargin > 1
    unknown = '*';
    tORu = 'T';
    lowercase = false;
    if rem(nargin,2)== 0
        error(message('bioinfo:int2nt:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'alphabet','unknownsymbol','case'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:int2nt:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:int2nt:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % alphabet
                    [isAminoAcid,isRNA] = bioinfoprivate.optAlphabet(pval,okargs{k}, mfilename);
                    if isRNA
                        tORu = 'U';
                    elseif isAminoAcid
                        error(message('bioinfo:int2nt:InvalidAlphabet'));
                    end
                case 2  % unknown
                    unknown = pval;
                    if ~ischar(pval)
                        error(message('bioinfo:int2nt:UnknownSymbolMustBeChar'));
                    end
                    if any(pval == 'acgtrykmswbdhvn')
                        error(message('bioinfo:int2nt:UnknownSymbolNotACGT'));
                    end
                case 3  % case
                    if ~isempty(strmatch(lower(pval),'lower'))
                        lowercase = true;
                    end
            end
        end
    end
    map = strrep(map,'*',unknown);
    map = strrep(map,'T',tORu);
    if lowercase
        map = lower(map);
    end
end

%sequence(sequence>16) = 0; %%% DELETED FOR SPEED
%sequence(sequence<0) = 0; %%% DELETED FOR SPEED
seqLength = length(sequence);


% pre-allocate the char

nt = char(0);nt(1,seqLength) = char(0);

for count = 1:seqLength
    nt(count) = map(double(sequence(count))+1);
end
% nt = reshape(nt,origsize); %%% DELETED BECAUSE DEPENDENT ON EARLIER
                             %%% DELETED CODE
