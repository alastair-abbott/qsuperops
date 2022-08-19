function is_valid_convQCFO = is_convQCFO(Wr, dims, parties, tol)
%is_convQCFO determines if input is a convex mixture of QC-FO process with different orders
%   is_valid_convQCFO = is_convQCFO(Wr, dims, parties, tol) is true if W is a valid conv(QC-FO) (including superinstrument)
%   If W is an sdpvar, then the function returns the yalmip constraints for Wr to be valid
%   
%   The arguments are specified as follows:
%       - W: the (potential) process matrix, a d-d matrix, or a cell of d-d matrices (may be sdpvars)
%       - dims: a vector specifying the dimensions of individual spaces and satisfying prod(dims) == d
%       - parties: a cell-array specifying the spaces that correspond to each party
%       - tol: the numerical tolerance for the validity check (default: 1e-6)

% Written by Alastair Abbott, last modified 19 August 2022

    % default tolerance
    if ~exist('tol','var')
        tol = 1e-6;
    end

    % First put Wr in canonical ordering (this checks the input validity too)
    % The spaces P,AI,AO,...,F then correspond to dims 1,2,3,...,2*N+2
    if exist('parties','var')
        [Wr, dims, parties] = superop_to_canonical_ordering(Wr, dims, parties);
    else
        [Wr, dims, parties] = superop_to_canonical_ordering(Wr, dims);
    end
    
    % Treat process matrices as single element superinstrument
    if ~iscell(Wr)
        Wr = {Wr};
    end
    R = length(Wr); % number of superinstrument elements
    N = length(parties) - 2;
   
    d_O = prod(dims(1:2:(2*N+1)));
    
    W = zeros(prod(dims));
    for r = 1:R
        assert(all(prod(dims) == size(Wr{r})), 'Error: W size doesn''t match provided dimensions.');
        W = W + Wr{r};
    end

    if isa(W,'sdpvar')
        is_valid_convQCFO = [trace(W) == d_O, superop_in_convQCFO_cone(Wr,dims,parties)];
    else
        % input is an explicit W; check the random robustness
        is_valid_convQCFO = (abs(trace(W) - d_O) <= tol) ...
            && superop_random_robustness(Wr,dims,parties,'convQCFO',sdpsettings('solver','mosek','verbose',0)) < tol;
    end
end

