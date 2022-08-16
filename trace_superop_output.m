function [Wr_new, dims_new, parties_new] = trace_superop_output(Wr,dims,parties,F_sys)
%trace_superop_output Traces out output spaces of a process matrix or superinstrument
%   [Wr_new, dims_new, parties_new] = trace_superop_output(Wr,dims,parties,F_sys) traces out the systems in F_sys
%
%   Can be applied either to a superoperator or a superinstrrument, in which case spaces traced
%   from each probabilistic circuit making up the superinstrument
%
% Require QETLAB for PermuteSystems, PartialTrace

% Written by Alastair Abbott, last modified 23 April 2021

    input_is_process_matrix = false;
    if ~iscell(Wr)
       Wr = {Wr};
       input_is_process_matrix = true;
    end
    
    F = parties{end}{1};
    assert(length(F_sys) <= length(F),'Fsys should specify which systems in F to be measured.');
    
    d_F_sys = prod(dims(F(F_sys)));
    
    % We actually treat this as "measuring" the specified systems with deterministic POVM
    M = eye(d_F_sys,d_F_sys);
    [Wr_new, dims_new, parties_new] = measureOutputOfProcess(Wr,dims,parties,M,F_sys);
    
    if input_is_process_matrix
        Wr_new = Wr_new{1};
    end
end

