function cone_constraints = superop_in_QCSup_cone(Wr, dims, parties)
%superop_in_QCSup_cone Yalmip constraints for a superinstrument to be in the cone of Quantum circuits with Superposition of orders
%   cone_constraints = superop_in_QCSup_cone(Wr, dims, parties)
%   Wr can be either a superinstrument or a superoperator/process matrix
%   Returns the yalmip SDP constraints for Wr to be in the convex hull of cones of valid Quantum circuits with Superposition of orders
%   This set is effectively the purification of convQCQC and was described in
%   Qiushi Liu, Zihao Hu, Haidong Yuan and Yuxiang Yang, arXiv:2203.09758 [quant-ph]
%   Note: this doesn't enforce the normalisation of the superoperator

% Written by Alastair Abbott 2022, last modified 31 August 2022

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
    
    %% We can do this for arbitrary N     

    % Need every Wr{i} >= 0
    cone_constraints = superop_in_PSD_cone(Wr);

    W = Wr{1};
    for r = 2:R
       W = W + Wr{r}; 
    end

    % Trace out F
    [W_reduced, dims_reduced, parties_reduced] = trace_superop_output(W,dims,parties,1);

    cone_constraints = [cone_constraints, superop_in_convQCFO_cone(W_reduced, dims_reduced, parties_reduced)];

end

