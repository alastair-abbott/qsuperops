function cone_constraints = superop_in_QCFO_dual_cone(Sr, dims, parties)
%superop_in_QCFO_dual_cone checks whether a set of operators is in the dual cone of QC-FOs (fixed order processes)
%   cone_constraints = superop_in_QCFO_dual_cone(Sr, dims, parties)
%   Sr can be either a witness or a set of witness elements
%   The order is assumed to be that in which the parties are specified
%   Note: this doesn't check/enforce the normalisation of the superoperator
%   
%   Formulation based on: J. Wechs, H. Dourdent, A. A. Abbott, C. Branciard, PRX Quantum 2, 030335 (2021)
%   (See Propositions 2 and 10 therein.)
%
% Requires QETLAB for PermuteSystems, PartialTrace

% Written by Alastair Abbott (2022), last modified 29 August 2022

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

    %% We do this generically for arbitrary N

    cone_constraints = [];

    P = 1;
    d_P = dims(P);

    S_k1kN = sdpvar(d,d,'hermitian','complex');

    % Define the PSD parts of the witness
    T_r = cell(1,R);
    for r = 1:R
        T_r{r} = Sr{r} - S_k1kN;
        cone_constraints = [cone_constraints, T_r{r} >= 0];
    end

    % Now we check that S_k1kN is in the right space
    S_proj = S_k1kN;
    for n = N:-1:1
        S_proj = S_proj - (tr_replace(S_proj,(2*n+2):(2*N+2),dims) - tr_replace(S_proj,(2*n+1):(2*N+2),dims));
    end

    if d_P ~= 1
        S_proj = S_proj - (tr_replace(S_proj,2:(2*N+2),dims) - tr_replace(S_proj,1:(2*N+2),dims));
    end

    cone_constraints = [cone_constraints, S_proj == 0];
end

