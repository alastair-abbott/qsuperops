function is_valid = is_valid_superop(Wr, dims, parties, tol)
%is_valid_process determines if input is a valid process matrix or superinstrument
%   is_valid = is_valid_superop(W, dims, parties, tol) is true if W is a valid process matrix or superinstrument.
%   If W is an sdpvar, then the function returns the yalmip constraints for Wr to be valid
%   
%   The arguments are specified as follows:
%       - W: the (potential) process matrix, a d-d matrix, or a cell of d-d matrices (may be sdpvars)
%       - dims: a vector specifying the dimensions of individual spaces and satisfying prod(dims) == d
%       - parties: a cell-array specifying the spaces that correspond to each party
%       - tol: the numerical tolerance for the validity check (default: 1e-6)
%
%   Note: parties{1}{2} specifies P, parties{n+1}{1/2} specifies A_n^{I/O}, parties{N+2}{1} specifies F

% Written by Alastair Abbott, last modified 19 April 2021

    % default tolerance
    if nargin == 3
        tol = 1e-6;
    end
    
    % Treat process matrices as single element superinstrument
    if ~iscell(Wr)
        Wr = {Wr};
    end
    R = length(Wr); % number of superinstrument elements

    % If Wr is an sdpvar then we want to create the sdpvar constraints instead
    input_is_sdpvar = false;
    for r = 1:R
        if isa(Wr{r},'sdpvar')
            input_is_sdpvar = true;
            break
        end
    end
    
    constraints_PSD = [];
    W = zeros(prod(dims),prod(dims));
    for r = 1:R
        assert(all(prod(dims) == size(Wr{r})), 'Error: W size doesn''t match provided dimensions.');
        % Each element of superinstrument should be PSD
        if input_is_sdpvar
            constraints_PSD = [constraints_PSD, Wr{r} >= 0];
        else
            constraints_PSD = [constraints_PSD, all(eig(Wr{r}) > -tol)];
        end
        W = W + Wr{r};
    end

    N = length(parties) - 2; % number of parties/operations
    
    % We do for each N seperately

    switch N
        case 1
            %%
            P = parties{1}{2};
            AI = parties{2}{1};
            AO = parties{2}{2};
            F = parties{3}{1};
            
            d_O = prod(dims([P, AO])); % product of output dims
            
            % Project W onto the space of valid processes
            Wproj = W - (tr_replace(W,F,dims) - tr_replace(W,[AO,F],dims)); % (1-AO)F
            Wproj = Wproj - (tr_replace(Wproj,[AI,AO,F],dims) - tr_replace(Wproj,[P,AI,AO,F],dims)); % (1-P)AF
              
        case 2
            %%
            P = parties{1}{2};
            AI = parties{2}{1};
            AO = parties{2}{2};
            A = [AI, AO];
            BI = parties{3}{1};
            BO = parties{3}{2};
            B = [BI, BO];
            F = parties{4}{1};
            
            d_O = prod(dims([P, AO, BO])); % product of output dims
            
            % Project W onto the space of valid processes
            Wproj = W - (tr_replace(W,[B,F],dims) - tr_replace(W,[AO,B,F],dims)); % (1-AO)BF
            Wproj = Wproj - (tr_replace(Wproj,[A,F],dims) - tr_replace(Wproj,[A,BO,F],dims)); % (1-BO)AF
            Wproj = Wproj - (tr_replace(Wproj,F,dims) - tr_replace(Wproj,[AO,F],dims) ...
                             - tr_replace(Wproj,[BO,F],dims) + tr_replace(Wproj,[AO,BO,F],dims)); % (1-AO)(1-BO)F
            Wproj = Wproj - (tr_replace(Wproj,[A,B,F],dims) - tr_replace(Wproj,[P,A,B,F],dims)); % (1-P)ABF
           
        case 3
            %%
            P = parties{1}{2};
            AI = parties{2}{1};
            AO = parties{2}{2};
            A = [AI, AO];
            BI = parties{3}{1};
            BO = parties{3}{2};
            B = [BI, BO];
            CI = parties{4}{1};
            CO = parties{4}{2};
            C = [CI, CO];
            F = parties{5}{1};
            
            d_O = prod(dims([P,AO,BO,CO])); % product of output dims
            
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
            P = parties{1}{2};
            AI = parties{2}{1};
            AO = parties{2}{2};
            A = [AI, AO];
            BI = parties{3}{1};
            BO = parties{3}{2};
            B = [BI, BO];
            CI = parties{4}{1};
            CO = parties{4}{2};
            C = [CI, CO];
            DI = parties{5}{1};
            DO = parties{5}{2};
            D = [DI, DO];
            F = parties{6}{1};
            
            d_O = prod(dims([P,AO,BO,CO,DO])); % product of output dims
            
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
            is_valid = false;
            return
    end
    
    if input_is_sdpvar
        % In some cases we might be able to know even an SDP var isn't a valid process
        % e.g., if it's not normalised properly
        if isa(trace(W) - d_O,'double') && abs(trace(W) - d_O) > tol
            is_valid = false;
        elseif isa(W - Wproj,'double') && ~matrix_is_equal(W-Wproj,zeros(prod(dims)),tol)
            is_valid = false;
        else 
            is_valid = [constraints_PSD, trace(W) == d_O, W == Wproj];
        end
    else
        is_PSD = all(constraints_PSD); % every W{r} is PSD
        is_normalised = abs(trace(W) - d_O) < tol;
        is_in_valid_space = matrix_is_equal(W,Wproj,tol);
        is_valid = is_PSD && is_normalised && is_in_valid_space;
    end
end

