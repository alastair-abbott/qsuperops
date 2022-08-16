function [W_new, dims_new, parties_new] = insert_superop_input(W,dims,parties,rho,P_sys)
%insert_superop_input New superoperator from old by plugging in inputs
%   [W_new, dims_new, parties_new] = insert_superop_input(W,dims,parties,rho,P_sys) inserts rho into specified systems of P
%
%   Given a superoperator W, generate a new one by inserting input rho
%   into the spaces specified in P_sys (as indices *within* P).
%   Also works for superinstruments, in which case state inserted into each probabilistic process
%
% Require QETLAB for PermuteSystems, PartialTrace

% Written by Alastair Abbott, last modified 28 April 2021

    % if input a process matrix, treat as trivial superinstrument
    input_is_process_matrix = false;
    if ~iscell(W)
        input_is_process_matrix = true;
        W = {W};
    end
    R = length(W); % number of superinstrument elements

    P = parties{1}{2};
    assert(length(P_sys) <= length(P), 'Psys should specify which systems in P state to be inserted in.');
   
    d = prod(dims);
    dPsys = prod(dims(P(P_sys)));
    assert(all(size(rho) == dPsys), 'rho has incorrect dimensions.');
    
    % First construct rho \otimes Id, and permute spaces to correct place
    permutation = 1:length(dims);
    permutation(P(P_sys)) = [];
    permutation = [P(P_sys), permutation];
    rhoId = PermuteSystems(Tensor(transpose(rho),eye(d/dPsys)), permutation, dims(permutation), 0, 1);
    W_new = cell(1,R);
    for r = 1:R
        W_new{r} = PartialTrace(rhoId*W{r}, P(P_sys), dims);
    end
    
    dims_new = dims;
    dims_new(P(P_sys)) = [];
    
    parties_new = parties;
    parties_new{1}{2}(P_sys) = [];
    % Need to readjust the space labels
    for k = 1:length(parties_new)
        for j = 1:2
           for i = 1:length(parties_new{k}{j})
               x = parties_new{k}{j}(i);
               parties_new{k}{j}(i) = x - sum(x > P(P_sys));
           end
        end
    end
    
    if input_is_process_matrix
        W_new = W_new{1};
    end
end

