function cone_constraints = superop_in_QCFO_dual_cone(Sr, dims, parties)
%superop_in_QCFO_dual_cone checks whether a set of operators is in the dual cone of QC-FOs (fixed order processes)
%   cone_constraints = superop_in_QCFO_dual_cone(Sr, dims, parties)
%   Sr can be either a witness or a set of witness elements
%   The order is assumed to be that in which the parties are specified
%   Note: this doesn't check/enforce the normalisation of the superoperator
%   
%   Formulation based on: J. Wechs, A. A. Abbott, C. Branciard, New J. Phys. 21, 013027 (2019)
%
% Requires QETLAB for PermuteSystems, PartialTrace

% Written by Alastair Abbott (2022), last modified 24 November 2022

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

    % Keep the logical and yalmip constraints separate until the end
    constraints_logical = true;
    constraints_yalmip = true;

    S_k1kN = sdpvar(d,d,'hermitian','complex');

    % Define the PSD parts of the witness
    T_r = cell(1,R);
    for r = 1:R
        T_r{r} = Sr{r} - S_k1kN;
    end
    constraints_yalmip = [constraints_yalmip, superop_in_PSD_cone(T_r)];

    %% Cone constraints
    S_proj = project_onto_QCFOs(S_k1kN,dims,parties);
    
    if isa(S_proj,'sdpvar')
        constraints_yalmip = [constraints_yalmip, S_proj == 0];
    else
        constraints_logical = [constraints_logical, matrix_is_equal(S_proj,zeros(d),tol)];
    end

    %% Combine the two types of constraints
    constraints_logical = all(constraints_logical);
    if constraints_logical == false
        cone_constraints = false;
    else
        cone_constraints = constraints_yalmip;
    end
end

