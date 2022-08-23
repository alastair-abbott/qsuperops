function in_valid_cone = superop_in_valid_cone(Wr, dims, parties,tol)
%superop_in_valid_cone checks whether an operator is in the cone of valid superinstruments
%   in_valid_cone = superop_in_valid_cone(Wr, dims, parties[,tol])
%   Wr can be either a superinstrument or a superoperator/process matrix
%   If Wr is an sdpvar and cone membership not trivially true/false, this
%   returns the yalmip constraints for Wr to be in the cone of valid superinstruments
%   Default tolerance is 1e-6, and irrelevant for sdpvar constraints
%   Note: this doesn't check/enforce the normalisation
%
% Requires QETLAB for PartialTrace

% Written by Alastair Abbott (2022), last modified 18 August 2022

    %% Setup and process the input

    % default tolerance
    if ~exist('tol','var')
        tol = 1e-6;
    end

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
    
    % Keep the logical and yalmip constraints separate until the end
    constraints_logical = true;
    constraints_yalmip = true;

    %% Common for all N

    R = length(Wr);
    N = length(parties) - 2;
    d = prod(dims);

    % First we check each Wr{r} >= 0
    constraints_temp = superop_in_PSD_cone(Wr,tol);
    if isa(constraints_temp,'logical')
        constraints_logical = [constraints_logical, constraints_temp];
    else
        constraints_yalmip = [constraints_yalmip, constraints_temp];
    end

    W = zeros(prod(dims));
    for r = 1:R
        assert(all(prod(dims) == size(Wr{r})), 'Error: W size doesn''t match provided dimensions.');
        W = W + Wr{r};
    end
    
    %% Now the specific constraints for each N seperately

    switch N
        case 1
            %%
            P = 1;
            AI = 2;
            AO = 3;
            F = 4;
            
            % Project W onto the space of valid processes
            Wproj = W - (tr_replace(W,F,dims) - tr_replace(W,[AO,F],dims)); % (1-AO)F
            Wproj = Wproj - (tr_replace(Wproj,[AI,AO,F],dims) - tr_replace(Wproj,[P,AI,AO,F],dims)); % (1-P)AF
              
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
            
            % Project W onto the space of valid processes
            Wproj = W - (tr_replace(W,[B,F],dims) - tr_replace(W,[AO,B,F],dims)); % (1-AO)BF
            Wproj = Wproj - (tr_replace(Wproj,[A,F],dims) - tr_replace(Wproj,[A,BO,F],dims)); % (1-BO)AF
            Wproj = Wproj - (tr_replace(Wproj,F,dims) - tr_replace(Wproj,[AO,F],dims) ...
                             - tr_replace(Wproj,[BO,F],dims) + tr_replace(Wproj,[AO,BO,F],dims)); % (1-AO)(1-BO)F
            Wproj = Wproj - (tr_replace(Wproj,[A,B,F],dims) - tr_replace(Wproj,[P,A,B,F],dims)); % (1-P)ABF
           
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
            
            Wproj = W - (tr_replace(W,[B,C,F],dims) - tr_replace(W,[AO,B,C,F],dims)); % (1-AO)BCF
            Wproj = Wproj - (tr_replace(Wproj,[A,C,F],dims) - tr_replace(Wproj,[A,BO,C,F],dims)); % (1-BO)ACF
            Wproj = Wproj - (tr_replace(Wproj,[A,B,F],dims) - tr_replace(Wproj,[A,B,CO,F],dims)); % (1-CO)ABF
            
            Wproj = Wproj - (tr_replace(Wproj,[C,F],dims) - tr_replace(Wproj,[AO,C,F],dims) ...
                             - tr_replace(Wproj,[BO,C,F],dims) + tr_replace(Wproj,[AO,BO,C,F],dims)); % (1-AO)(1-BO)CF
            Wproj = Wproj - (tr_replace(Wproj,[B,F],dims) - tr_replace(Wproj,[AO,B,F],dims) ...
                             - tr_replace(Wproj,[B,CO,F],dims) + tr_replace(Wproj,[AO,B,CO,F],dims)); % (1-AO)B(1-CO)F
            Wproj = Wproj - (tr_replace(Wproj,[A,F],dims) - tr_replace(Wproj,[A,BO,F],dims) ...
                             - tr_replace(Wproj,[A,CO,F],dims) + tr_replace(Wproj,[A,BO,CO,F],dims)); % A(1-BO)(1-CO)F
            
            Wproj = Wproj - (tr_replace(Wproj,F,dims) - tr_replace(Wproj,[AO,F],dims) - tr_replace(Wproj,[BO,F],dims) - tr_replace(Wproj,[CO,F],dims) ...
                             + tr_replace(Wproj,[AO,BO,F],dims) + tr_replace(Wproj,[AO,CO,F],dims) + tr_replace(Wproj,[BO,CO,F],dims) ...
                             - tr_replace(Wproj,[AO,BO,CO,F],dims)); % (1-AO)(1-BO)(1-CO)F
            
            Wproj = Wproj - (tr_replace(Wproj,[A,B,C,F],dims) - tr_replace(Wproj,[P,A,B,C,F],dims)); % (1-P)ABCF
            
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
            
            Wproj = W - (tr_replace(W,[B,C,D,F],dims) - tr_replace(W,[AO,B,C,D,F],dims)); % (1-AO)BCDF
            Wproj = Wproj - (tr_replace(Wproj,[A,C,D,F],dims) - tr_replace(Wproj,[A,BO,C,D,F],dims)); % (1-BO)ACDF
            Wproj = Wproj - (tr_replace(Wproj,[A,B,D,F],dims) - tr_replace(Wproj,[A,B,CO,D,F],dims)); % (1-CO)ABDF
            Wproj = Wproj - (tr_replace(Wproj,[A,B,C,F],dims) - tr_replace(Wproj,[A,B,C,DO,F],dims)); % (1-DO)ABCF
            
            Wproj = Wproj - (tr_replace(Wproj,[C,D,F],dims) - tr_replace(Wproj,[AO,C,D,F],dims) ...
                             - tr_replace(Wproj,[BO,C,D,F],dims) + tr_replace(Wproj,[AO,BO,C,D,F],dims)); % (1-AO)(1-BO)CDF
            Wproj = Wproj - (tr_replace(Wproj,[B,D,F],dims) - tr_replace(Wproj,[AO,B,D,F],dims) ...
                             - tr_replace(Wproj,[B,CO,D,F],dims) + tr_replace(Wproj,[AO,B,CO,D,F],dims)); % (1-AO)B(1-CO)DF
            Wproj = Wproj - (tr_replace(Wproj,[B,C,F],dims) - tr_replace(Wproj,[AO,B,C,F],dims) ...
                             - tr_replace(Wproj,[B,C,DO,F],dims) + tr_replace(Wproj,[AO,B,C,DO,F],dims)); % (1-AO)BC(1-DO)F
            Wproj = Wproj - (tr_replace(Wproj,[A,D,F],dims) - tr_replace(Wproj,[A,BO,D,F],dims) ...
                             - tr_replace(Wproj,[A,CO,D,F],dims) + tr_replace(Wproj,[A,BO,CO,D,F],dims)); % A(1-BO)(1-CO)DF
            Wproj = Wproj - (tr_replace(Wproj,[A,C,F],dims) - tr_replace(Wproj,[A,BO,C,F],dims) ...
                             - tr_replace(Wproj,[A,C,DO,F],dims) + tr_replace(Wproj,[A,BO,C,DO,F],dims)); % A(1-BO)C(1-DO)F
            Wproj = Wproj - (tr_replace(Wproj,[A,B,F],dims) - tr_replace(Wproj,[A,B,CO,F],dims) ...
                             - tr_replace(Wproj,[A,B,DO,F],dims) + tr_replace(Wproj,[A,B,CO,D0,F],dims)); % AB(1-CO)(1-DO)F
            
            Wproj = Wproj - (tr_replace(Wproj,[D,F],dims) - tr_replace(Wproj,[AO,D,F],dims) - tr_replace(Wproj,[BO,D,F],dims) - tr_replace(Wproj,[CO,D,F],dims) ...
                             + tr_replace(Wproj,[AO,BO,D,F],dims) + tr_replace(Wproj,[AO,CO,D,F],dims) + tr_replace(Wproj,[BO,CO,D,F],dims) ...
                             - tr_replace(Wproj,[AO,BO,CO,D,F],dims)); % (1-AO)(1-BO)(1-CO)DF
            Wproj = Wproj - (tr_replace(Wproj,[C,F],dims) - tr_replace(Wproj,[AO,C,F],dims) - tr_replace(Wproj,[BO,C,F],dims) - tr_replace(Wproj,[C,DO,F],dims) ...
                             + tr_replace(Wproj,[AO,BO,C,F],dims) + tr_replace(Wproj,[AO,C,DO,F],dims) + tr_replace(Wproj,[BO,C,DO,F],dims) ...
                             - tr_replace(Wproj,[AO,BO,C,DO,F],dims)); % (1-AO)(1-BO)C(1-DO)F
            Wproj = Wproj - (tr_replace(Wproj,[B,F],dims) - tr_replace(Wproj,[AO,B,F],dims) - tr_replace(Wproj,[B,CO,F],dims) - tr_replace(Wproj,[B,DO,F],dims) ...
                             + tr_replace(Wproj,[AO,B,CO,F],dims) + tr_replace(Wproj,[AO,B,DO,F],dims) + tr_replace(Wproj,[B,CO,DO,F],dims) ...
                             - tr_replace(Wproj,[AO,B,CO,DO,F],dims)); % (1-AO)B(1-CO)(1-DO)F
            Wproj = Wproj - (tr_replace(Wproj,[A,F],dims) - tr_replace(Wproj,[A,BO,F],dims) - tr_replace(Wproj,[A,CO,F],dims) - tr_replace(Wproj,[A,DO,F],dims) ...
                             + tr_replace(Wproj,[A,BO,CO,F],dims) + tr_replace(Wproj,[A,BO,DO,F],dims) + tr_replace(Wproj,[A,CO,DO,F],dims) ...
                             - tr_replace(Wproj,[A,BO,CO,DO,F],dims)); % A(1-BO)(1-CO)(1-DO)F

            Wproj = Wproj - (tr_replace(Wproj,[F],dims) - tr_replace(Wproj,[AO,F],dims) - tr_replace(Wproj,[BO,F],dims) - tr_replace(Wproj,[CO,F],dims) - tr_replace(Wproj,[DO,F],dims) ...
                             + tr_replace(Wproj,[AO,BO,F],dims) + tr_replace(Wproj,[AO,CO,F],dims) + tr_replace(Wproj,[AO,DO,F],dims) + tr_replace(Wproj,[BO,CO,F],dims) + tr_replace(Wproj,[BO,DO,F],dims) + tr_replace(Wproj,[CO,DO,F],dims) ...
                             - tr_replace(Wproj,[AO,BO,CO,F],dims) - tr_replace(Wproj,[AO,BO,DO,F],dims) - tr_replace(Wproj,[AO,CO,DO,F],dims) - tr_replace(Wproj,[BO,CO,DO,F],dims) ...
                             + tr_replace(Wproj,[AO,B0,CO,DO,F],dims)); % (1-AO)(1-BO)(1-CO)(1-DO)F
            
            Wproj = Wproj - (tr_replace(Wproj,[A,B,C,D,F],dims) - tr_replace(Wproj,[P,A,B,C,D,F],dims)); % (1-P)ABCDF
        otherwise
            disp('Check currently only implemented up to N=4.');
            in_valid_cone = false;
            return
    end

    diff_valid_space = W - Wproj;
    if isa(diff_valid_space,'sdpvar')
        constraints_yalmip = [constraints_yalmip, diff_valid_space == 0];
    else
        constraints_logical = [constraints_logical, matrix_is_equal(diff_valid_space,zeros(d),tol)];
    end

    %% Combine the two types of constraints
    constraints_logical = all(constraints_logical);
    if constraints_logical == false
        in_valid_cone = false;
    else
        in_valid_cone = constraints_yalmip;
    end
end

