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
    d = prod(dims);
    
    %% For parallel circuits it's easy to do for arbitrary N
    % In canonical ordering, input spaces are 2,4,... and output 3,5,...
    
    % If Wr is an sdpvar then we want to create the sdpvar constraints instead
    input_is_sdpvar = false;
    for r = 1:R
        if isa(Wr{r},'sdpvar')
            input_is_sdpvar = true;
            break
        end
    end

    constraints_PSD = [];
    W = zeros(d,d);
    for r = 1:R
        assert(all(d == size(Wr{r})), 'Error: W size doesn''t match provided dimensions.');
        % Each element of superinstrument should be PSD
        if input_is_sdpvar
            constraints_PSD = [constraints_PSD, Wr{r} >= 0];
        else
            constraints_PSD = [constraints_PSD, all(eig(Wr{r}) > -tol)];
        end
        W = W + Wr{r};
    end

    % If Wr is an sdpvar then we want to create the sdpvar constraints instead
    input_is_sdpvar = false;
    if isa(W,'sdpvar')
        input_is_sdpvar = true;
    end

    % All input and output spaces
    P = 1;
    AI = 2:2:2*N;
    AO = 3:2:(2*N+1);
    F = 2*N + 2;

    d_P = (dims(P));
    d_O = prod(dims(AO));

    if input_is_sdpvar
        in_QCPAR_cone = [W >= 0];
    else
        in_QCPAR_cone = [all(eig(Wr{r}) > -tol)];
    end


    W_I = 1/d_O*PartialTrace(W,[AO,F],dims);

    in_QCPAR_cone = [in_QCPAR_cone, PartialTrace(W,F,dims) == PermuteSystems(tensor_id(W_I,d_O),[P,AI,AO],dims([P,AI,AO]),0,1)];
    if d_P ~= 1 % Otherwise this is enforced by the normalisation of the input superinstrument
        % We only want it to be proportional to the identity, since we are only checking the cone
        in_QCPAR_cone = [in_QCPAR_cone, PartialTrace(W_I,2:(N+1),dims([P,AO])) == trace(W_I)/d_P*eye(d_P)];
    end

end

