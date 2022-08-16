function [Wr_new, dims_new, parties_new] = measure_superop_output(Wr,dims,parties,M,F_sys)
%measureOutputOfProcess Generates superinstrument by measuring output of a process
%   [Wr_new, dims_new, parties_new] = measure_superop_output(Wr,dims,parties,M,F_sys) measures M on systems in F_sys
%
%   Applies the POVM M to the systems F_sys of the superinstrument/operator Wr, generating a new
%   superinstrrument. 
%
% Require QETLAB for PermuteSystems, PartialTrace

% Written by Alastair Abbott, last modified 28 April 2021
   
    % convert to trivial superinstrrument if input a process
    if ~iscell(Wr)
        Wr = {Wr};
    end
    r_init = length(Wr);
    
    F = parties{end}{1};
    assert(length(F_sys) <= length(F),'Fsys should specify which systems in F to be measured.');
    
    d = prod(dims);
    d_F_sys = prod(dims(F(F_sys)));
    assert(all(size(M,[1,2]) == d_F_sys), 'M has incorrect dimensions.');
    
    % number of measurement outcomes
    R = size(M,3);
    
    M_total = zeros(d_F_sys,d_F_sys);
    for r = 1:R
        M_total = M_total + M(:,:,r);
        assert(min(eig(M(:,:,r))) > -1e-6,['POVM element M(:,:,', num2str(r),') not PSD.']);
    end
    assert(matrix_is_equal(M_total,eye(d_F_sys),1e-6),'POVM M does not sum to identity.');
    
    Wr_new = cell(1,R*r_init);
    % First construct Id \otimes M_r, and permute spaces to correct place
    for r = 1:R
        permutation = 1:length(dims);
        permutation(F(F_sys)) = [];
        permutation = [permutation, F(F_sys)];
        idMr = PermuteSystems(Tensor(eye(d/d_F_sys),transpose(M(:,:,r))), permutation, dims(permutation), 0, 1);
        for i = 1:r_init
            Wr_new{(i-1)*R+r} = PartialTrace(idMr*Wr{i}, F(F_sys), dims);
        end
    end
    
    dims_new = dims;
    dims_new(F(F_sys)) = [];
    
    parties_new = parties;
    parties_new{end}{1}(F_sys) = [];
    % Need to readjust the space labels
    % If spaces are give in same order as dim, this shouldn't do anything
    for k = 1:length(parties_new)
        for j = 1:2
           for i = 1:length(parties_new{k}{j})
               x = parties_new{k}{j}(i);
               parties_new{k}{j}(i) = x - sum(x > F(F_sys));
           end
        end
    end
end

