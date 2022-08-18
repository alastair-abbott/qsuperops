function in_QCPAR_cone = superop_in_QCPAR_cone(Wr, dims, parties,tol)
%superop_in_QCPAR_cone checks whether a superinstrument is in the cone of QC-PARs
%   in_QCPAR_cone = superop_in_QCPAR_cone(Wr, dims, parties[, tol])
%   Wr can be either a superinstrument or a superoperator/process matrix
%   If Wr is an sdpvar and cone membership not trivially true/false, this
%   returns the yalmip constraints for Wr to be in the cone of valid Quantum Circuits with Parallel operations
%   Default tolerance is 1e-6, and irrelevant for sdpvar constraints
%   Note: this doesn't check/enforce the normalisation of the superoperator
%   
%   Formulation based on: J. Wechs, H. Dourdent, A. A. Abbott, C. Branciard, PRX Quantum 2, 030335 (2021)
%   (See Propositions 3 and 11 therein.)
%
% Requires QETLAB for PermuteSystems, PartialTrace

% Written by Alastair Abbott (2022), last modified 18 August 2022

    % default tolerance
    if nargin == 3
        tol = 1e-6;
    end

    % First put Wr in canonical ordering (this checks the input validity too)
    % The spaces P,AI,AO,...,F then correspond to dims 1,2,3,...,2*N+2
    [Wr, dims, parties] = superop_to_canonical_ordering(Wr, dims, parties);

    % Treat a process matrix as a 1-element superinstrument
    if ~iscell(Wr)
        Wr = {Wr};
    end
    
    R = length(Wr);
    N = length(parties) - 2;
    
    %% For parallel circuits it's easy to do for arbitrary N
    % In canonical ordering, input spaces are 2,4,... and output 3,5,...

    % All input and output spaces
    P = 1;
    AI = 2:2:2*N;
    AO = 3:2:(2*N+1);
    F = 2*N + 2;

    d_P = prod(dims(P));
    d_AO = prod(dims(AO));

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

    W_I = 1/d_AO*PartialTrace(W,[AO,F],dims);

    diff_F_last = PartialTrace(W,F,dims) - PermuteSystems(tensor_id(W_I,d_AO),[P,AI,AO],dims([P,AI,AO]),0,1);
    if isa(diff_F_last,'sdpvar')
        constraints_yalmip = [constraints_yalmip, diff_F_last == 0];
    else
        constraints_logical = [constraints_logical, matrix_is_equal(diff_F_last,zeros(size(diff_F_last)),tol)];
    end

    if d_P ~= 1 % Otherwise this is enforced by the normalisation of the input superinstrument
        % We only want it to be proportional to the identity, since we are only checking the cone
        diff_P_first = PartialTrace(W_I,2:(N+1),dims([P,AO])) - trace(W_I)/d_P*eye(d_P);
        if isa(diff_F_last,'sdpvar')
            constraints_yalmip = [constraints_yalmip, diff_P_first == 0];
        else
            constraints_logical = [constraints_logical, matrix_is_equal(diff_P_first,zeros(d_P),tol)];
        end
    end

    % Combine the two types of constraints
    constraints_logical = all(constraints_logical);
    if constraints_logical == false
        in_QCPAR_cone = false;
    else
        in_QCPAR_cone = constraints_yalmip;
    end
end

