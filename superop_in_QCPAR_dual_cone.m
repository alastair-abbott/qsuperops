function cone_constraints = superop_in_QCPAR_dual_cone(Sr, dims, parties)
%superop_in_QCPAR_dual_cone checks whether a set of operators is in the dual cone of QC-PARs (parallel processes)
%   cone_constraints = superop_in_QCPAR_dual_cone(Sr, dims, parties)
%   Sr can be either a witness or a set of witness elements
%   If Sr is an sdpvar and cone membership not trivially true/false, this
%   returns the yalmip constraints for Sr to be in the dual cone of valid Quantum Circuits with parallel operations
%   Note: this doesn't check/enforce the normalisation of the superoperator
%   
%   Formulation based on: J. Wechs, H. Dourdent, A. A. Abbott, C. Branciard, PRX Quantum 2, 030335 (2021)
%   (See Propositions 2 and 10 therein.)
%
% Requires QETLAB for PermuteSystems, PartialTrace

% Written by Alastair Abbott (2022), last modified 30 August 2022

    % First put Sr in canonical ordering (this checks the input validity too)
    % The spaces P,AI,AO,...,F then correspond to dims 1,2,3,...,2*N+2
    if exist('parties','var') && ~isempty(parties)
        [Sr, dims, parties] = superop_to_canonical_ordering(Sr, dims, parties);
    else
        [Sr, dims, parties] = superop_to_canonical_ordering(Sr, dims);
    end

    % Treat a process matrix as a 1-element superinstrument
    if ~iscell(Sr)
        Sr = {Sr};
    end
    
    R = length(Sr);
    N = length(parties) - 2;

    d = prod(dims);

    %% For parallel circuits it's easy to do for arbitrary N
    % In canonical ordering, input spaces are 2,4,... and output 3,5,...

    cone_constraints = [];

    % All input and output spaces
    P = 1;
    AO = 3:2:(2*N+1);
    F = 2*N + 2;

    d_P = prod(dims(P));

    S = sdpvar(d,d,'hermitian','complex'); % the fixed part of the witness elements

    % Define the PSD parts of the witness
    T_r = cell(1,R);
    for r = 1:R
        T_r{r} = Sr{r} - S;
    end
    cone_constraints = [cone_constraints, superop_in_PSD_cone(T_r)];

    % Now we check that S is in the right space
    S_proj = S - (tr_replace(S,F,dims) - tr_replace(S,[AO,F],dims));

    if d_P ~= 1
        S_proj = S_proj - (tr_replace(S_proj,2:F,dims) - tr_replace(S_proj,1:F,dims));
    end

    cone_constraints = [cone_constraints, S_proj == 0];
end

