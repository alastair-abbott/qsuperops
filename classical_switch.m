function [W,dims,parties] = classical_switch(N,d_t,trace_F_t,input_state)
%classical_switch Generates the process matrix for the classical N-switch
%   [W,dims,parties] = classical_switch(N,d_t,trace_F_t,input_state)
%   N: number of parties
%   d_t: target dimension. Default: 2
%   trace_F_t: if true, throw away final target state. Default: false
%   input_state: an optional initial control-target input state.
%   
%   Ordering of systems is taken to be P_c,P_t,AI,AO,BI,BO,...,F_t,F_c
%
% Requires QETLAB for PermuteSystems, PartialTrace

% Written by Alastair Abbott 2022, last modified 19 August 2022

    %% Process input
    assert(N >= 2,'Error: Need to specify N >= 2 as input.');

    if nargin < 2
        d_t = 2; % default is qubit target
    end

    if nargin < 3
        trace_F_t = false;
    end

    %% 
    d_C = factorial(N);
    dims = [d_C, d_t*ones(1,2*N+2), d_C];
    parties = cell(1,N+2);
    parties{1} = {[],[1,2]}; % P = P_c P_t
    parties{end} = {[2*N+3, 2*N+4],[]}; % F = F_t F_c
    for i = 1:N
        parties{i+1} = {2*i+1, 2*i+2};
    end

    control_basis = eye(d_C);
    id_CJ = pure_CJ(eye(d_t)); % pure CJ of identity channel

    allPerms = reshape([2*perms(1:N)-1; 2*perms(1:N)],factorial(N),[]); % All permutaion of parties spaces
    allPerms = flip(allPerms); % More natural ordering of permutations, one per row

    % We build up the process vector |w>> one permutation at a time
    W = zeros(d_t^(2*N+2)*d_C^2);
    for i = 1:factorial(N)
        w_element = Tensor(control_basis(:,i),Tensor(id_CJ,N+1),control_basis(:,i));
        w_element = PermuteSystems(w_element,[1,2,allPerms(i,:)+2,2*N+3,2*N+4],dims,0,1);
        W = W + w_element*w_element';
    end

    if trace_F_t == true
        [W, dims, parties] = trace_superop_output(W,dims,parties,1);
    end

    if nargin == 4
        [W, dims, parties] = insert_superop_input(W,dims,parties,input_state,[1,2]);
    end
end