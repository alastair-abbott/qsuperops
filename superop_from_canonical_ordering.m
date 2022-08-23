function Wr = superop_from_canonical_ordering(Wr_canonical,dims,parties)
%superop_from_canonical_ordering Reorders spaces so that they're in the desired order
%   Wr = superop_from_canonical_ordering(Wr_canonical,dims,parties)
%
%   This is effectively the inverse of superop_from_canonical_ordering
%   The canonical groups the input and output spaces for each party into a single
%   space, and them orders them as P,AI,AO,BI,BO,...,F
%
%   Careful: the input is the desired dim and parties, the canonical ones
%   being canonically known!
%
% Requires QETLAB for PermuteSystems

% Written by Alastair Abbott, last modified 29 April 2021

    % Everything here works equally well for witnesses as for process matrices
    input_is_process_matrix = false;
    if ~iscell(Wr_canonical)
        input_is_process_matrix = true;
        Wr_canonical = {Wr_canonical};
    end
    
    R = length(Wr_canonical);
    N = length(parties) - 2;
    
    d = prod(dims);
    
    for i = 1:R
       assert(isequal(size(Wr_canonical{i}),[d,d]),'Process size doesn''t agree with specified dimensions.'); 
    end
    
    % We just apply the inverse permutation of the systems specified by dim to the
    % order given by parties
    Wr = cell(1,R);
    perm = [parties{1}{1}];
    for n = 2:N+1
       perm = [perm, parties{n}{1}, parties{n}{2}]; 
    end
    perm = [perm, parties{N+2}{1}];
    assert(length(perm) == length(dims), 'Incompatibility between dimensions and parties.');
    for i = 1:R
        % Now we use the inverse permutation
       Wr{i} = PermuteSystems(Wr_canonical{i},perm,dims,0,1); 
    end
    
    % return a process matrix if input was a process matrix rather than a superinstrument
    if input_is_process_matrix
        Wr = Wr{1};
    end

end

