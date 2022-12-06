function in_QCFO_cone = superop_in_QCFO_cone(Wr, dims, parties,tol)
%superop_in_QCFO_cone checks whether a superinstrument is in the cone of QC-FOs (fixed order processes)
%   in_QCFO_cone = superop_in_QCFO_cone(Wr, dims, parties[, tol])
%   Wr can be either a superinstrument or a superoperator/process matrix
%   If Wr is an sdpvar and cone membership not trivially true/false, this
%   returns the yalmip constraints for Wr to be in the cone of valid Quantum Circuits with Fixed-Order operations
%   The order is assumed to be that in which the parties are specified
%   Default tolerance is 1e-6, and irrelevant for sdpvar constraints
%   Note: this doesn't check/enforce the normalisation of the superoperator
%
% Requires QETLAB for PermuteSystems, PartialTrace

% Written by Alastair Abbott (2022), last modified 24 November 2022

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
    
    R = length(Wr);
    N = length(parties) - 2;
    d = prod(dims);
    
    %% Setup and check the elements are PSD

    % Keep the logical and yalmip constraints separate until the end
    constraints_logical = [true];
    constraints_yalmip = [true];

    % First we check each Wr{r} >= 0
    constraints_temp = superop_in_PSD_cone(Wr,tol);
    if isa(constraints_temp,'logical')
        constraints_logical = [constraints_logical, constraints_temp];
    else
        constraints_yalmip = [constraints_yalmip, constraints_temp];
    end

    W = Wr{1};
    for r = 2:R
       W = W + Wr{r}; 
    end

    % Project W onto space of QCFOs
    Wproj = project_onto_QCFOs(W,dims,parties);

    diff_valid_space = W - Wproj;
    if isa(diff_valid_space,'sdpvar')
        constraints_yalmip = [constraints_yalmip, diff_valid_space == 0];
    else
        constraints_logical = [constraints_logical, matrix_is_equal(diff_valid_space,zeros(d),tol)];
    end

    % Combine the two types of constraints
    constraints_logical = all(constraints_logical);
    if constraints_logical == false
        in_QCFO_cone = false;
    else
        in_QCFO_cone = constraints_yalmip;
    end
end

