function [Wr_canonical,dims_canonical,parties_canonical] = superop_to_canonical_ordering(Wr,dims,parties)
%operator_to_canonical_ordering Reorders spaces so that they're canonically ordered
%   [Wr_canonical,dims_canonical,parties_canonical] = superop_to_canonical_ordering(Wr,dims,parties)
%
%   In general a process, witness or superinstrument can be specified on arbitrary spaces, with
%   arbitrary ordering (and a space can be broken into subspaces).
%   operator_to_canonical_ordering transforms any such operator into one with spaces ordered as
%   P,AI,AO,BI,BO,...,F, and all subspaces grouped together
%
% Requires QETLAB for PermuteSystems

% Written by Alastair Abbott, last modified 28 April 2021

    % Everything here works equally well for witnesses as for process matrices
    input_is_process_matrix = false;
    if ~iscell(Wr)
        input_is_process_matrix = true;
        Wr = {Wr};
    end
    
    R = length(Wr);
    N = length(parties) - 2;
    
    d = prod(dims);
    
    for i = 1:R
       assert(isequal(size(Wr{i}),[d,d]),'Process size doesn''t agree with specified dimensions.'); 
    end
    
    % We want some canonical dimensions [dP dAI dAO ... dF]
    dims_canonical = zeros(1,2*N+2);
    parties_canonical = cell(1,N+2);
    
    dims_canonical(1) = prod(dims(parties{1}{2}));
    parties_canonical{1} = {[],1};
    for n = 1:N
        dims_canonical(2*n) = prod(dims(parties{n+1}{1}));
        dims_canonical(2*n+1) = prod(dims(parties{n+1}{2}));
        parties_canonical{n+1} = {2*n,2*n+1};
    end
    dims_canonical(end) = prod(dims(parties{N+2}{1}));
    parties_canonical{end} = {2*N+2,[]};
    
    dI = prod(dims_canonical(2:2:end));
    dO = prod(dims_canonical(1:2:end));
    assert(d == dI*dO, 'Incompatibility between dimensions and parties.');
    
    % We will permute the process so that the systems are ordered canonically
    % This will make things much easier later, since we can assume the ordering
    Wr_canonical = cell(1,R);
    perm = [];
    for n = 1:N+2
       perm = [perm, parties{n}{1}, parties{n}{2}]; 
    end
    assert(length(perm) == length(dims), 'Incompatibility between dimensions and parties.');
    for i = 1:R
       Wr_canonical{i} = PermuteSystems(Wr{i},perm,dims); 
    end
    
    % return a process matrix if input was a process matrix rather than a superinstrument
    if input_is_process_matrix
        Wr_canonical = Wr_canonical{1};
    end

end

