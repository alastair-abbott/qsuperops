function is_valid = is_valid_process(Wr, dims, parties, tol)
%is_valid_process determines if input is a valid process matrix or superinstrument
%   is_valid = is_valid_process(W, dims, parties, tol) is true if W is a valid process matrix or superinstrument.
%   
%   The arguments are specified as follows:
%       - W: the (potential) process matrix, a d-d matrix, or a cell of d-d matrices
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
    
    is_PSD = true;
    W = zeros(prod(dims),prod(dims));
    for r = 1:R
        assert(all(prod(dims) == size(Wr{r})), 'Error: W size doesn''t match provided dimensions.');
        % Each element of superinstrument should be PSD
        is_PSD = is_PSD && all(eig(Wr{r}) > -tol);
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
            
            dO = prod(dims([P, AO])); % product of output dims
            
            % Project W onto the space of valid processes
            Wproj = W - (tr_replace(W,F,dims) - tr_replace(W,[AO,F],dims)); % (1-AO)F
            Wproj = Wproj - (tr_replace(Wproj,[AI,AO,F],dims) - tr_replace(Wproj,[P,AI,AO,F],dims)); % (1-P)AF
            
            is_normalised = abs(trace(W) - dO) < tol;
            is_in_valid_space = matrix_is_equal(W,Wproj,tol);
              
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
            
            dO = prod(dims([P, AO, BO])); % product of output dims
            
            % Project W onto the space of valid processes
            Wproj = W - (tr_replace(W,[B,F],dims) - tr_replace(W,[AO,B,F],dims)); % (1-AO)BF
            Wproj = Wproj - (tr_replace(Wproj,[A,F],dims) - tr_replace(Wproj,[A,BO,F],dims)); % (1-BO)AF
            Wproj = Wproj - (tr_replace(Wproj,F,dims) - tr_replace(Wproj,[AO,F],dims) ...
                             - tr_replace(Wproj,[BO,F],dims) + tr_replace(Wproj,[AO,BO,F],dims)); % (1-AO)(1-BO)F
            Wproj = Wproj - (tr_replace(Wproj,[A,B,F],dims) - tr_replace(Wproj,[P,A,B,F],dims)); % (1-P)ABF
            
            is_normalised = abs(trace(W) - dO) < tol;
            is_in_valid_space = matrix_is_equal(W,Wproj,tol);
           
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
            
            dO = prod(dims([P,AO,BO,CO])); % product of output dims
            
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
            
            is_normalised = abs(trace(W) - dO) < tol;
            is_in_valid_space = matrix_is_equal(W,Wproj,tol);
            
        otherwise
            disp('Check currently not implemented for this number of parties.');
            is_in_valid_space = false;
    end
    
    is_valid = is_PSD && is_normalised && is_in_valid_space;
end

