function W_proj = project_onto_valid_superops(W, dims_raw, parties_raw)
%project_onto_valid_superops Projects a superoperator onto the subspace of valid superoperators
%   W_projected = project_onto_superop_subspace(W, dims, parties)
%   W is taken to be a superoperator (not a superinstrument)
%
% Requires QETLAB for PartialTrace

% Written by Alastair Abbott (2022), last modified 20 September 2022

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
    
    %% Implement projection separately and explicitly for each N
    switch N
        case 1
            %%
            P = 1;
            AI = 2;
            AO = 3;
            F = 4;
            
            W_proj = W - (tr_replace(W,F,dims) - tr_replace(W,[AO,F],dims)); % (1-AO)F
            W_proj = W_proj - (tr_replace(W_proj,[AI,AO,F],dims) - tr_replace(W_proj,[P,AI,AO,F],dims)); % (1-P)AF
              
        case 2
            %%
            P = 1;
            AI = 2;
            AO = 3;
            A = [AI, AO];
            BI = 4;
            BO = 5;
            B = [BI, BO];
            F = 6;
            
            W_proj = W - (tr_replace(W,[B,F],dims) - tr_replace(W,[AO,B,F],dims)); % (1-AO)BF
            W_proj = W_proj - (tr_replace(W_proj,[A,F],dims) - tr_replace(W_proj,[A,BO,F],dims)); % (1-BO)AF
            W_proj = W_proj - (tr_replace(W_proj,F,dims) - tr_replace(W_proj,[AO,F],dims) ...
                             - tr_replace(W_proj,[BO,F],dims) + tr_replace(W_proj,[AO,BO,F],dims)); % (1-AO)(1-BO)F
            W_proj = W_proj - (tr_replace(W_proj,[A,B,F],dims) - tr_replace(W_proj,[P,A,B,F],dims)); % (1-P)ABF
           
        case 3
            %%
            P = 1;
            AI = 2;
            AO = 3;
            A = [AI, AO];
            BI = 4;
            BO = 5;
            B = [BI, BO];
            CI = 6;
            CO = 7;
            C = [CI, CO];
            F = 8;
            
            W_proj = W - (tr_replace(W,[B,C,F],dims) - tr_replace(W,[AO,B,C,F],dims)); % (1-AO)BCF
            W_proj = W_proj - (tr_replace(W_proj,[A,C,F],dims) - tr_replace(W_proj,[A,BO,C,F],dims)); % (1-BO)ACF
            W_proj = W_proj - (tr_replace(W_proj,[A,B,F],dims) - tr_replace(W_proj,[A,B,CO,F],dims)); % (1-CO)ABF
            
            W_proj = W_proj - (tr_replace(W_proj,[C,F],dims) - tr_replace(W_proj,[AO,C,F],dims) ...
                             - tr_replace(W_proj,[BO,C,F],dims) + tr_replace(W_proj,[AO,BO,C,F],dims)); % (1-AO)(1-BO)CF
            W_proj = W_proj - (tr_replace(W_proj,[B,F],dims) - tr_replace(W_proj,[AO,B,F],dims) ...
                             - tr_replace(W_proj,[B,CO,F],dims) + tr_replace(W_proj,[AO,B,CO,F],dims)); % (1-AO)B(1-CO)F
            W_proj = W_proj - (tr_replace(W_proj,[A,F],dims) - tr_replace(W_proj,[A,BO,F],dims) ...
                             - tr_replace(W_proj,[A,CO,F],dims) + tr_replace(W_proj,[A,BO,CO,F],dims)); % A(1-BO)(1-CO)F
            
            W_proj = W_proj - (tr_replace(W_proj,F,dims) - tr_replace(W_proj,[AO,F],dims) - tr_replace(W_proj,[BO,F],dims) - tr_replace(W_proj,[CO,F],dims) ...
                             + tr_replace(W_proj,[AO,BO,F],dims) + tr_replace(W_proj,[AO,CO,F],dims) + tr_replace(W_proj,[BO,CO,F],dims) ...
                             - tr_replace(W_proj,[AO,BO,CO,F],dims)); % (1-AO)(1-BO)(1-CO)F
            
            W_proj = W_proj - (tr_replace(W_proj,[A,B,C,F],dims) - tr_replace(W_proj,[P,A,B,C,F],dims)); % (1-P)ABCF
            
        case 4
            %%
            P = 1;
            AI = 2;
            AO = 3;
            A = [AI, AO];
            BI = 4;
            BO = 5;
            B = [BI, BO];
            CI = 6;
            CO = 7;
            C = [CI, CO];
            DI = 8;
            DO = 9;
            D = [DI, DO];
            F = 10;
            
            W_proj = W - (tr_replace(W,[B,C,D,F],dims) - tr_replace(W,[AO,B,C,D,F],dims)); % (1-AO)BCDF
            W_proj = W_proj - (tr_replace(W_proj,[A,C,D,F],dims) - tr_replace(W_proj,[A,BO,C,D,F],dims)); % (1-BO)ACDF
            W_proj = W_proj - (tr_replace(W_proj,[A,B,D,F],dims) - tr_replace(W_proj,[A,B,CO,D,F],dims)); % (1-CO)ABDF
            W_proj = W_proj - (tr_replace(W_proj,[A,B,C,F],dims) - tr_replace(W_proj,[A,B,C,DO,F],dims)); % (1-DO)ABCF
            
            W_proj = W_proj - (tr_replace(W_proj,[C,D,F],dims) - tr_replace(W_proj,[AO,C,D,F],dims) ...
                             - tr_replace(W_proj,[BO,C,D,F],dims) + tr_replace(W_proj,[AO,BO,C,D,F],dims)); % (1-AO)(1-BO)CDF
            W_proj = W_proj - (tr_replace(W_proj,[B,D,F],dims) - tr_replace(W_proj,[AO,B,D,F],dims) ...
                             - tr_replace(W_proj,[B,CO,D,F],dims) + tr_replace(W_proj,[AO,B,CO,D,F],dims)); % (1-AO)B(1-CO)DF
            W_proj = W_proj - (tr_replace(W_proj,[B,C,F],dims) - tr_replace(W_proj,[AO,B,C,F],dims) ...
                             - tr_replace(W_proj,[B,C,DO,F],dims) + tr_replace(W_proj,[AO,B,C,DO,F],dims)); % (1-AO)BC(1-DO)F
            W_proj = W_proj - (tr_replace(W_proj,[A,D,F],dims) - tr_replace(W_proj,[A,BO,D,F],dims) ...
                             - tr_replace(W_proj,[A,CO,D,F],dims) + tr_replace(W_proj,[A,BO,CO,D,F],dims)); % A(1-BO)(1-CO)DF
            W_proj = W_proj - (tr_replace(W_proj,[A,C,F],dims) - tr_replace(W_proj,[A,BO,C,F],dims) ...
                             - tr_replace(W_proj,[A,C,DO,F],dims) + tr_replace(W_proj,[A,BO,C,DO,F],dims)); % A(1-BO)C(1-DO)F
            W_proj = W_proj - (tr_replace(W_proj,[A,B,F],dims) - tr_replace(W_proj,[A,B,CO,F],dims) ...
                             - tr_replace(W_proj,[A,B,DO,F],dims) + tr_replace(W_proj,[A,B,CO,D0,F],dims)); % AB(1-CO)(1-DO)F
            
            W_proj = W_proj - (tr_replace(W_proj,[D,F],dims) - tr_replace(W_proj,[AO,D,F],dims) - tr_replace(W_proj,[BO,D,F],dims) - tr_replace(W_proj,[CO,D,F],dims) ...
                             + tr_replace(W_proj,[AO,BO,D,F],dims) + tr_replace(W_proj,[AO,CO,D,F],dims) + tr_replace(W_proj,[BO,CO,D,F],dims) ...
                             - tr_replace(W_proj,[AO,BO,CO,D,F],dims)); % (1-AO)(1-BO)(1-CO)DF
            W_proj = W_proj - (tr_replace(W_proj,[C,F],dims) - tr_replace(W_proj,[AO,C,F],dims) - tr_replace(W_proj,[BO,C,F],dims) - tr_replace(W_proj,[C,DO,F],dims) ...
                             + tr_replace(W_proj,[AO,BO,C,F],dims) + tr_replace(W_proj,[AO,C,DO,F],dims) + tr_replace(W_proj,[BO,C,DO,F],dims) ...
                             - tr_replace(W_proj,[AO,BO,C,DO,F],dims)); % (1-AO)(1-BO)C(1-DO)F
            W_proj = W_proj - (tr_replace(W_proj,[B,F],dims) - tr_replace(W_proj,[AO,B,F],dims) - tr_replace(W_proj,[B,CO,F],dims) - tr_replace(W_proj,[B,DO,F],dims) ...
                             + tr_replace(W_proj,[AO,B,CO,F],dims) + tr_replace(W_proj,[AO,B,DO,F],dims) + tr_replace(W_proj,[B,CO,DO,F],dims) ...
                             - tr_replace(W_proj,[AO,B,CO,DO,F],dims)); % (1-AO)B(1-CO)(1-DO)F
            W_proj = W_proj - (tr_replace(W_proj,[A,F],dims) - tr_replace(W_proj,[A,BO,F],dims) - tr_replace(W_proj,[A,CO,F],dims) - tr_replace(W_proj,[A,DO,F],dims) ...
                             + tr_replace(W_proj,[A,BO,CO,F],dims) + tr_replace(W_proj,[A,BO,DO,F],dims) + tr_replace(W_proj,[A,CO,DO,F],dims) ...
                             - tr_replace(W_proj,[A,BO,CO,DO,F],dims)); % A(1-BO)(1-CO)(1-DO)F

            W_proj = W_proj - (tr_replace(W_proj,[F],dims) - tr_replace(W_proj,[AO,F],dims) - tr_replace(W_proj,[BO,F],dims) - tr_replace(W_proj,[CO,F],dims) - tr_replace(W_proj,[DO,F],dims) ...
                             + tr_replace(W_proj,[AO,BO,F],dims) + tr_replace(W_proj,[AO,CO,F],dims) + tr_replace(W_proj,[AO,DO,F],dims) + tr_replace(W_proj,[BO,CO,F],dims) + tr_replace(W_proj,[BO,DO,F],dims) + tr_replace(W_proj,[CO,DO,F],dims) ...
                             - tr_replace(W_proj,[AO,BO,CO,F],dims) - tr_replace(W_proj,[AO,BO,DO,F],dims) - tr_replace(W_proj,[AO,CO,DO,F],dims) - tr_replace(W_proj,[BO,CO,DO,F],dims) ...
                             + tr_replace(W_proj,[AO,B0,CO,DO,F],dims)); % (1-AO)(1-BO)(1-CO)(1-DO)F
            
            W_proj = W_proj - (tr_replace(W_proj,[A,B,C,D,F],dims) - tr_replace(W_proj,[P,A,B,C,D,F],dims)); % (1-P)ABCDF
        otherwise
            error('Projector onto valid superoperators currently only implemented up to N=4.');
    end

    % Put W back in its original ordering
    W_proj = superop_from_canonical_ordering(W_proj,dims_raw,parties_raw);

end


