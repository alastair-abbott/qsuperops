function cone_constraints = superop_in_valid_dual_cone(Sr, dims, parties, tol)
%superop_in_valid_dual_cone checks whether an operator is in the dual cone of valid superinstruments
%   superop_in_valid_dual_cone = superop_in_valid_dual_cone(Sr, dims, parties[, tol])
%   Sr can be either a witness or a set of witness elements
%   Note: this doesn't check/enforce the normalisation
%
% Requires QETLAB for PartialTrace

% Written by Alastair Abbott (2022), last modified 30 August 2022

    % default tolerance
    if ~exist('tol','var')
        tol = 1e-6;
    end

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
    
    %% First we set up the things common to all N

    % Keep the logical and yalmip constraints separate until the end
    constraints_logical = true;
    constraints_yalmip = true;

    Sv = sdpvar(d,d,'hermitian','complex'); % For the validity constraints

    % Define the PSD parts of the witness
    T_r = cell(1,R);
    for r = 1:R
        T_r{r} = Sr{r} - Sv;
    end
    constraints_yalmip = [constraints_yalmip, superop_in_PSD_cone(T_r)];

    %% Cone constraints
    Sv_proj = project_onto_valid_superops(Sv,dims,parties);
    
    if isa(Sv_proj,'sdpvar')
        constraints_yalmip = [constraints_yalmip, Sv_proj == 0];
    else
        constraints_logical = [constraints_logical, matrix_is_equal(Sv_proj,zeros(d),tol)];
    end

    %% Combine the two types of constraints
    constraints_logical = all(constraints_logical);
    if constraints_logical == false
        cone_constraints = false;
    else
        cone_constraints = constraints_yalmip;
    end
end

