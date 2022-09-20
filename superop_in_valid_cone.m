function in_valid_cone = superop_in_valid_cone(Wr, dims, parties,tol)
%superop_in_valid_cone checks whether an operator is in the cone of valid superinstruments
%   in_valid_cone = superop_in_valid_cone(Wr, dims, parties[,tol])
%   Wr can be either a superinstrument or a superoperator/process matrix
%   If Wr is an sdpvar and cone membership not trivially true/false, this
%   returns the yalmip constraints for Wr to be in the cone of valid superinstruments
%   Default tolerance is 1e-6, and irrelevant for sdpvar constraints
%   Note: this doesn't check/enforce the normalisation
%
% Requires QETLAB for PartialTrace

% Written by Alastair Abbott (2022), last modified 20 September 2022

    %% Setup and process the input

    % default tolerance
    if ~exist('tol','var')
        tol = 1e-6;
    end

    % First put Wr in canonical ordering (this checks the input validity too)
    % The spaces P,AI,AO,...,F then correspond to dims 1,2,3,...,2*N+2
    if exist('parties','var') && ~isempty(parties)
        [Wr, dims, parties] = superop_to_canonical_ordering(Wr, dims, parties);
    else
        [Wr, dims, parties] = superop_to_canonical_ordering(Wr, dims);
    end

    % Treat a process matrix as a 1-element superinstrument
    if ~iscell(Wr)
        Wr = {Wr};
    end
    
    % Keep the logical and yalmip constraints separate until the end
    constraints_logical = true;
    constraints_yalmip = true;

    R = length(Wr);
    N = length(parties) - 2;
    d = prod(dims);

    %% Cone constraints

    % First we check each Wr{r} >= 0
    constraints_temp = superop_in_PSD_cone(Wr,tol);
    if isa(constraints_temp,'logical')
        constraints_logical = [constraints_logical, constraints_temp];
    else
        constraints_yalmip = [constraints_yalmip, constraints_temp];
    end

    W = zeros(prod(dims));
    for r = 1:R
        assert(all(prod(dims) == size(Wr{r})), 'Error: W size doesn''t match provided dimensions.');
        W = W + Wr{r};
    end
    
    % Project W onto space of valid superoperators
    Wproj = project_onto_valid_superops(W,dims,parties);

    diff_valid_space = W - Wproj;
    if isa(diff_valid_space,'sdpvar')
        constraints_yalmip = [constraints_yalmip, diff_valid_space == 0];
    else
        constraints_logical = [constraints_logical, matrix_is_equal(diff_valid_space,zeros(d),tol)];
    end

    %% Combine the two types of constraints
    constraints_logical = all(constraints_logical);
    if constraints_logical == false
        in_valid_cone = false;
    else
        in_valid_cone = constraints_yalmip;
    end
end

