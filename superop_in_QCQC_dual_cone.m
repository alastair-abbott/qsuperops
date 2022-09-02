function cone_constraints = superop_in_QCQC_dual_cone(Sr, dims, parties)
%superop_in_QCQC_dual_cone Yalmip constraints for a set of matrices to be in the cone of witnesses for QC-QCs
%   cone_constraints = superop_in_QCQC_dual_cone(Sr, dims, parties) 
%   Sr can be either a single witness S, or a superinstrument witness Sr
%   Returns the yalmip constraints, i.e. for the dual cone
%
% Requires QETLAB for PermuteSystems, PartialTrace

% Written by Alastair Abbott 2022, last modified 30 August 2022

    % First put Sr in canonical ordering (this checks the input validity too)
    % The spaces P,AI,AO,...,F then correspond to dims 1,2,3,...,2*N+2
    if exist('parties','var') && ~isempty(parties)
        [Sr, dims, parties] = superop_to_canonical_ordering(Sr, dims, parties);
    else
        [Sr, dims, parties] = superop_to_canonical_ordering(Sr, dims);
    end

    % treat a process matrix witness as a single element superinstrument witness
    if ~iscell(Sr)
        Sr = {Sr};
    end
    
    R = length(Sr);
    N = length(parties) - 2;
    
    d = prod(dims);
    
    %% We treat each N separately
    
    cone_constraints = [];
    
    switch N
        case 1
            %%
            % For N=1, QCCC = QCQC = QCFO = QCPAR
            cone_constraints = [cone_constraints, superop_in_QCCC_dual_cone(Sr,dims,parties)];
            
        case 2
            %%
            % For N=2, W is a QCQC if tr_F(W) is a QC-CC
            % Implement this in the dual picture
            F = 2*N+2;

            S_tot = Sr{1};
            for r = 2:R
                S_tot = S_tot + Sr{r};
            end
            cone_constraints = [cone_constraints, S_tot - tr_replace(S_tot,F,dims) == 0];

            [Sr_F, dims_F, parties_F] = trace_superop_output(Sr,dims,parties,1);
            cone_constraints = [cone_constraints, superop_in_QCCC_dual_cone(Sr_F,dims_F,parties_F)];
            
        case 3
            %%
            
            
        case 4
            %%
            
            
        otherwise
            error('Currently only implemented up to N=4');
    end
end

