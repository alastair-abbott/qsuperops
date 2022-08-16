function isValid = isValidProcess(Wr, dim, parties, tol)
%isValidProcess determines if input is a valid process matrix or superinstrument
%   isValid = isValidProcess(W, dim, parties, tol) is true if W is a valid process matrix or superinstrument.
%   
%   The arguments are specified as follows:
%       - W: the (potential) process matrix, a d-d matrix, or a cell of d-d matrices
%       - dim: a vector specifying the dimensions of individual spaces and satisfying prod(dim) == d
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
    
    isPSD = true;
    W = zeros(prod(dim),prod(dim));
    for r = 1:R
        assert(all(prod(dim) == size(Wr{r})), 'Error: W size doesn''t match provided dimensions.');
        % Each element of superinstrument should be PSD
        isPSD = isPSD && all(eig(Wr{r}) > -tol);
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
            
            dO = prod(dim([P, AO])); % product of output dims
            
            % Project W onto the space of valid processes
            Wproj = W - (trRep(W,F,dim) - trRep(W,[AO,F],dim)); % (1-AO)F
            Wproj = Wproj - (trRep(Wproj,[AI,AO,F],dim) - trRep(Wproj,[P,AI,AO,F],dim)); % (1-P)AF
            
            isNormalised = abs(trace(W) - dO) < tol;
            isInValidSpace = matrixIsEqual(W,Wproj,tol);
              
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
            
            dO = prod(dim([P, AO, BO])); % product of output dims
            
            % Project W onto the space of valid processes
            Wproj = W - (trRep(W,[B,F],dim) - trRep(W,[AO,B,F],dim)); % (1-AO)BF
            Wproj = Wproj - (trRep(Wproj,[A,F],dim) - trRep(Wproj,[A,BO,F],dim)); % (1-BO)AF
            Wproj = Wproj - (trRep(Wproj,F,dim) - trRep(Wproj,[AO,F],dim) ...
                             - trRep(Wproj,[BO,F],dim) + trRep(Wproj,[AO,BO,F],dim)); % (1-AO)(1-BO)F
            Wproj = Wproj - (trRep(Wproj,[A,B,F],dim) - trRep(Wproj,[P,A,B,F],dim)); % (1-P)ABF
            
            isNormalised = abs(trace(W) - dO) < tol;
            isInValidSpace = matrixIsEqual(W,Wproj,tol);
           
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
            
            dO = prod(dim([P,AO,BO,CO])); % product of output dims
            
            Wproj = W - (trRep(W,[B,C,F],dim) - trRep(W,[AO,B,C,F],dim)); % (1-AO)BCF
            Wproj = Wproj - (trRep(Wproj,[A,C,F],dim) - trRep(Wproj,[A,BO,C,F],dim)); % (1-BO)ACF
            Wproj = Wproj - (trRep(Wproj,[A,B,F],dim) - trRep(Wproj,[A,B,CO,F],dim)); % (1-CO)ABF
            
            Wproj = Wproj - (trRep(Wproj,[C,F],dim) - trRep(Wproj,[AO,C,F],dim) ...
                             - trRep(Wproj,[BO,C,F],dim) + trRep(Wproj,[AO,BO,C,F],dim)); % (1-AO)(1-BO)CF
            Wproj = Wproj - (trRep(Wproj,[B,F],dim) - trRep(Wproj,[AO,B,F],dim) ...
                             - trRep(Wproj,[B,CO,F],dim) + trRep(Wproj,[AO,B,CO,F],dim)); % (1-AO)B(1-CO)F
            Wproj = Wproj - (trRep(Wproj,[A,F],dim) - trRep(Wproj,[A,BO,F],dim) ...
                             - trRep(Wproj,[A,CO,F],dim) + trRep(Wproj,[A,BO,CO,F],dim)); % A(1-BO)(1-CO)F
            
            Wproj = Wproj - (trRep(Wproj,F,dim) - trRep(Wproj,[AO,F],dim) - trRep(Wproj,[BO,F],dim) - trRep(Wproj,[CO,F],dim) ...
                             + trRep(Wproj,[AO,BO,F],dim) + trRep(Wproj,[AO,CO,F],dim) + trRep(Wproj,[BO,CO,F],dim) ...
                             - trRep(Wproj,[AO,BO,CO,F],dim)); % (1-AO)(1-BO)(1-CO)F
            
            Wproj = Wproj - (trRep(Wproj,[A,B,C,F],dim) - trRep(Wproj,[P,A,B,C,F],dim)); % (1-P)ABCF
            
            isNormalised = abs(trace(W) - dO) < tol;
            isInValidSpace = matrixIsEqual(W,Wproj,tol);
            
        otherwise
            disp('Check currently not implemented for this number of parties.');
            isInValidSpace = false;
    end
    
    isValid = isPSD && isNormalised && isInValidSpace;
end

