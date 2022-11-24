function W_proj = project_onto_QCFOs(W, dims_raw, parties_raw)
%project_onto_QCFOs Projects a superoperator onto the subspace of QCFOs
%   W_projected = project_onto_QCFOs(W, dims, parties)
%   W is taken to be a superoperator (not a superinstrument)
%
% Requires QETLAB for PartialTrace
%
% Rather than the definition of QC-FOs given in J. Wechs, H. Dourdent, A. A. Abbott, C. Branciard, PRX Quantum 2, 030335 (2021) 
% (see Propositions 2 and 10 therein), we make use of the conditions specified on the whole process matrix as defined in
% J. Wechs, A. A. Abbott, C. Branciard, New J. Phys. 21, 013027 (2019)

% Written by Alastair Abbott (2022), last modified 24 November 2022

    %% Setup and process the input

    % First put W in canonical ordering (this checks the input validity too)
    % The spaces P,AI,AO,...,F then correspond to dims 1,2,3,...,2*N+2
    if exist('parties_raw','var') && ~isempty(parties_raw)
        [W, dims, parties] = superop_to_canonical_ordering(W, dims_raw, parties_raw);
    else
        [W, dims, parties] = superop_to_canonical_ordering(W, dims_raw);
        parties_raw = parties;
    end

    N = length(parties) - 2;
    
    %% We do this generically for arbitrary N
    F = 2*N+2;

    W_proj = W; % Start with the full W

    % Then project on each constraint (1-An_O)A_rest,F
    for n = N:-1:0 % We go to party 0, which is P
        W_proj = W_proj - (tr_replace(W_proj,2*n+2:F,dims) - tr_replace(W_proj,2*n+1:F,dims));
    end

    % Put W back in its original ordering
    W_proj = superop_from_canonical_ordering(W_proj,dims_raw,parties_raw);
end


