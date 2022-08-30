function cone_constraints = superop_in_valid_dual_cone(Sr, dims, parties, tol)
%superop_in_valid_dual_cone checks whether an operator is in the dual cone of valid superinstruments
%   superop_in_valid_dual_cone = superop_in_valid_dual_cone(Sr, dims, parties[, tol])
%   Sr can be either a witness or a set of witness elements
%   Note: this doesn't check/enforce the normalisation
%
% Requires QETLAB for PartialTrace

% Written by Alastair Abbott (2022), last modified 30 August 2022

    % default tolerance
    if ~exist('tol','var')
        tol = 1e-6;
    end

    % First put Sr in canonical ordering (this checks the input validity too)
    % The spaces P,AI,AO,...,F then correspond to dims 1,2,3,...,2*N+2
    if exist('parties','var') && ~isempty(parties)
        [Sr, dims, parties] = superop_to_canonical_ordering(Sr, dims, parties);
    else
        [Sr, dims, parties] = superop_to_canonical_ordering(Sr, dims);
    end

    % Treat a process matrix as a 1-element superinstrument
    if ~iscell(Sr)
        Sr = {Sr};
    end
    
    R = length(Sr);
    N = length(parties) - 2;

    d = prod(dims);
    
    %% First we set up the things common to all N

    % Keep the logical and yalmip constraints separate until the end
    constraints_logical = true;
    constraints_yalmip = true;

    Sv = sdpvar(d,d,'hermitian','complex'); % For the validity constraints

    % Define the PSD parts of the witness
    T_r = cell(1,R);
    for r = 1:R
        T_r{r} = Sr{r} - Sv;
    end
    constraints_yalmip = [constraints_yalmip, superop_in_PSD_cone(T_r)];

    %% Then the specific constraints on Sv for each N

    switch N
        case 1
            %%
            P = 1;
            AI = 2;
            AO = 3;
            F = 4;
            
            % Project Sv onto space of valid processes
            Sv_proj = Sv - (tr_replace(Sv,F,dims) - tr_replace(Sv,[AO,F],dims)); % (1-AO)F
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[AI,AO,F],dims) - tr_replace(Sv_proj,[P,AI,AO,F],dims)); % (1-P)AF
              
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
            
            % Project Sv onto the space of valid processes
            Sv_proj = Sv - (tr_replace(Sv,[B,F],dims) - tr_replace(Sv,[AO,B,F],dims)); % (1-AO)BF
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[A,F],dims) - tr_replace(Sv_proj,[A,BO,F],dims)); % (1-BO)AF
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,F,dims) - tr_replace(Sv_proj,[AO,F],dims) ...
                             - tr_replace(Sv_proj,[BO,F],dims) + tr_replace(Sv_proj,[AO,BO,F],dims)); % (1-AO)(1-BO)F
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[A,B,F],dims) - tr_replace(Sv_proj,[P,A,B,F],dims)); % (1-P)ABF
           
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
            
            Sv_proj = Sv - (tr_replace(Sv,[B,C,F],dims) - tr_replace(Sv,[AO,B,C,F],dims)); % (1-AO)BCF
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[A,C,F],dims) - tr_replace(Sv_proj,[A,BO,C,F],dims)); % (1-BO)ACF
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[A,B,F],dims) - tr_replace(Sv_proj,[A,B,CO,F],dims)); % (1-CO)ABF
            
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[C,F],dims) - tr_replace(Sv_proj,[AO,C,F],dims) ...
                             - tr_replace(Sv_proj,[BO,C,F],dims) + tr_replace(Sv_proj,[AO,BO,C,F],dims)); % (1-AO)(1-BO)CF
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[B,F],dims) - tr_replace(Sv_proj,[AO,B,F],dims) ...
                             - tr_replace(Sv_proj,[B,CO,F],dims) + tr_replace(Sv_proj,[AO,B,CO,F],dims)); % (1-AO)B(1-CO)F
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[A,F],dims) - tr_replace(Sv_proj,[A,BO,F],dims) ...
                             - tr_replace(Sv_proj,[A,CO,F],dims) + tr_replace(Sv_proj,[A,BO,CO,F],dims)); % A(1-BO)(1-CO)F
            
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,F,dims) - tr_replace(Sv_proj,[AO,F],dims) - tr_replace(Sv_proj,[BO,F],dims) - tr_replace(Sv_proj,[CO,F],dims) ...
                             + tr_replace(Sv_proj,[AO,BO,F],dims) + tr_replace(Sv_proj,[AO,CO,F],dims) + tr_replace(Sv_proj,[BO,CO,F],dims) ...
                             - tr_replace(Sv_proj,[AO,BO,CO,F],dims)); % (1-AO)(1-BO)(1-CO)F
            
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[A,B,C,F],dims) - tr_replace(Sv_proj,[P,A,B,C,F],dims)); % (1-P)ABCF
            
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
            
            Sv_proj = Sv - (tr_replace(Sv,[B,C,D,F],dims) - tr_replace(Sv,[AO,B,C,D,F],dims)); % (1-AO)BCDF
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[A,C,D,F],dims) - tr_replace(Sv_proj,[A,BO,C,D,F],dims)); % (1-BO)ACDF
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[A,B,D,F],dims) - tr_replace(Sv_proj,[A,B,CO,D,F],dims)); % (1-CO)ABDF
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[A,B,C,F],dims) - tr_replace(Sv_proj,[A,B,C,DO,F],dims)); % (1-DO)ABCF
            
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[C,D,F],dims) - tr_replace(Sv_proj,[AO,C,D,F],dims) ...
                             - tr_replace(Sv_proj,[BO,C,D,F],dims) + tr_replace(Sv_proj,[AO,BO,C,D,F],dims)); % (1-AO)(1-BO)CDF
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[B,D,F],dims) - tr_replace(Sv_proj,[AO,B,D,F],dims) ...
                             - tr_replace(Sv_proj,[B,CO,D,F],dims) + tr_replace(Sv_proj,[AO,B,CO,D,F],dims)); % (1-AO)B(1-CO)DF
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[B,C,F],dims) - tr_replace(Sv_proj,[AO,B,C,F],dims) ...
                             - tr_replace(Sv_proj,[B,C,DO,F],dims) + tr_replace(Sv_proj,[AO,B,C,DO,F],dims)); % (1-AO)BC(1-DO)F
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[A,D,F],dims) - tr_replace(Sv_proj,[A,BO,D,F],dims) ...
                             - tr_replace(Sv_proj,[A,CO,D,F],dims) + tr_replace(Sv_proj,[A,BO,CO,D,F],dims)); % A(1-BO)(1-CO)DF
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[A,C,F],dims) - tr_replace(Sv_proj,[A,BO,C,F],dims) ...
                             - tr_replace(Sv_proj,[A,C,DO,F],dims) + tr_replace(Sv_proj,[A,BO,C,DO,F],dims)); % A(1-BO)C(1-DO)F
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[A,B,F],dims) - tr_replace(Sv_proj,[A,B,CO,F],dims) ...
                             - tr_replace(Sv_proj,[A,B,DO,F],dims) + tr_replace(Sv_proj,[A,B,CO,D0,F],dims)); % AB(1-CO)(1-DO)F
            
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[D,F],dims) - tr_replace(Sv_proj,[AO,D,F],dims) - tr_replace(Sv_proj,[BO,D,F],dims) - tr_replace(Sv_proj,[CO,D,F],dims) ...
                             + tr_replace(Sv_proj,[AO,BO,D,F],dims) + tr_replace(Sv_proj,[AO,CO,D,F],dims) + tr_replace(Sv_proj,[BO,CO,D,F],dims) ...
                             - tr_replace(Sv_proj,[AO,BO,CO,D,F],dims)); % (1-AO)(1-BO)(1-CO)DF
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[C,F],dims) - tr_replace(Sv_proj,[AO,C,F],dims) - tr_replace(Sv_proj,[BO,C,F],dims) - tr_replace(Sv_proj,[C,DO,F],dims) ...
                             + tr_replace(Sv_proj,[AO,BO,C,F],dims) + tr_replace(Sv_proj,[AO,C,DO,F],dims) + tr_replace(Sv_proj,[BO,C,DO,F],dims) ...
                             - tr_replace(Sv_proj,[AO,BO,C,DO,F],dims)); % (1-AO)(1-BO)C(1-DO)F
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[B,F],dims) - tr_replace(Sv_proj,[AO,B,F],dims) - tr_replace(Sv_proj,[B,CO,F],dims) - tr_replace(Sv_proj,[B,DO,F],dims) ...
                             + tr_replace(Sv_proj,[AO,B,CO,F],dims) + tr_replace(Sv_proj,[AO,B,DO,F],dims) + tr_replace(Sv_proj,[B,CO,DO,F],dims) ...
                             - tr_replace(Sv_proj,[AO,B,CO,DO,F],dims)); % (1-AO)B(1-CO)(1-DO)F
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[A,F],dims) - tr_replace(Sv_proj,[A,BO,F],dims) - tr_replace(Sv_proj,[A,CO,F],dims) - tr_replace(Sv_proj,[A,DO,F],dims) ...
                             + tr_replace(Sv_proj,[A,BO,CO,F],dims) + tr_replace(Sv_proj,[A,BO,DO,F],dims) + tr_replace(Sv_proj,[A,CO,DO,F],dims) ...
                             - tr_replace(Sv_proj,[A,BO,CO,DO,F],dims)); % A(1-BO)(1-CO)(1-DO)F

            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[F],dims) - tr_replace(Sv_proj,[AO,F],dims) - tr_replace(Sv_proj,[BO,F],dims) - tr_replace(Sv_proj,[CO,F],dims) - tr_replace(Sv_proj,[DO,F],dims) ...
                             + tr_replace(Sv_proj,[AO,BO,F],dims) + tr_replace(Sv_proj,[AO,CO,F],dims) + tr_replace(Sv_proj,[AO,DO,F],dims) + tr_replace(Sv_proj,[BO,CO,F],dims) + tr_replace(Sv_proj,[BO,DO,F],dims) + tr_replace(Sv_proj,[CO,DO,F],dims) ...
                             - tr_replace(Sv_proj,[AO,BO,CO,F],dims) - tr_replace(Sv_proj,[AO,BO,DO,F],dims) - tr_replace(Sv_proj,[AO,CO,DO,F],dims) - tr_replace(Sv_proj,[BO,CO,DO,F],dims) ...
                             + tr_replace(Sv_proj,[AO,B0,CO,DO,F],dims)); % (1-AO)(1-BO)(1-CO)(1-DO)F
            
            Sv_proj = Sv_proj - (tr_replace(Sv_proj,[A,B,C,D,F],dims) - tr_replace(Sv_proj,[P,A,B,C,D,F],dims)); % (1-P)ABCDF
        otherwise
            error('Currently only implemented up to N=4');
    end

    if isa(Sv_proj,'sdpvar')
        constraints_yalmip = [constraints_yalmip, Sv_proj == 0];
    else
        constraints_logical = [constraints_logical, matrix_is_equal(Sv_proj,zeros(d),tol)];
    end

    %% Combine the two types of constraints
    constraints_logical = all(constraints_logical);
    if constraints_logical == false
        cone_constraints = false;
    else
        cone_constraints = constraints_yalmip;
    end
end

