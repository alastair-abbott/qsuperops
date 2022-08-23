function [Wr_canonical,dims_canonical,parties_canonical] = superop_to_canonical_ordering(Wr,dims,parties)
%operator_to_canonical_ordering Reorders spaces so that they're canonically ordered
%   [Wr_canonical,dims_canonical,parties_canonical] = superop_to_canonical_ordering(Wr,dims,parties)
%
%   In general a process, witness or superinstrument can be specified on arbitrary spaces, with
%   arbitrary ordering (and a space can be broken into subspaces).
%   operator_to_canonical_ordering transforms any such operator into one with spaces ordered as
%   P,AI,AO,BI,BO,...,F, and all subspaces grouped together
%
%   If parties is not specified, dims will be interpreted as specifying the dims of spaces
%   P,AI,AO,...,F. Using this shortcut is not recommended!
%
% Requires QETLAB for PermuteSystems

% Written by Alastair Abbott, last modified 28 April 2021

    % Check parties are properly specified.
    if ~exist('parties','var') || isempty(parties)
        N = length(dims)/2 - 1;
        assert(mod(N,1) == 0, 'Error: please specify parties or ensure dimensions of P and F specified (even if trivial)');
        disp(['Warning: parties not specified. Interpreting superoperator as ', num2str(N), '-partite operator.']);
        parties = cell(1,N+2);
        parties{1} = {1};
        for n = 1:N
            parties{n+1} = {2*n, 2*n+1};
        end
        parties{end} = {2*N+2};
    end

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
    
    % Catch case where P is empty by accident
    if isempty(parties{1})
        parties{1} = {[]};
    end
    dims_canonical(1) = prod(dims(parties{1}{1}));

    parties_canonical{1} = {1};
    for n = 1:N
        dims_canonical(2*n) = prod(dims(parties{n+1}{1}));
        dims_canonical(2*n+1) = prod(dims(parties{n+1}{2}));
        parties_canonical{n+1} = {2*n,2*n+1};
    end
    if isempty(parties{N+2})
        parties{N+2} = {};
    end
    dims_canonical(end) = prod(dims(parties{N+2}{1}));
    parties_canonical{end} = {2*N+2};
    
    dI = prod(dims_canonical(2:2:end));
    dO = prod(dims_canonical(1:2:end));
    assert(d == dI*dO, 'Incompatibility between dimensions and parties.');
    
    % We will permute the process so that the systems are ordered canonically
    % This will make things much easier later, since we can assume the ordering
    Wr_canonical = cell(1,R);
    perm = [parties{1}{1}];
    for n = 2:N+1
       perm = [perm, parties{n}{1}, parties{n}{2}]; 
    end
    perm = [perm, parties{N+2}{1}];
    assert(length(perm) == length(dims), 'Incompatibility between dimensions and parties.');
    for i = 1:R
       Wr_canonical{i} = PermuteSystems(Wr{i},perm,dims); 
    end
    
    % return a process matrix if input was a process matrix rather than a superinstrument
    if input_is_process_matrix
        Wr_canonical = Wr_canonical{1};
    end

end

