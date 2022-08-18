function is_QCPAR = is_QCPAR(Wr, dims, parties, tol)
%is_QCPAR determines if input is a QC-PAR process
%   is_QCPAR = is_valid_QCPAR(W, dims, parties, tol) is true if W is a valid QC-PAR (including superinstrument)
%   If W is an sdpvar, then the function returns the yalmip constraints for Wr to be valid
%   
%   The arguments are specified as follows:
%       - W: the (potential) process matrix, a d-d matrix, or a cell of d-d matrices (may be sdpvars)
%       - dims: a vector specifying the dimensions of individual spaces and satisfying prod(dims) == d
%       - parties: a cell-array specifying the spaces that correspond to each party
%       - tol: the numerical tolerance for the validity check (default: 1e-6)

% Written by Alastair Abbott, last modified 18 August 2022

    % default tolerance
    if nargin == 3
        tol = 1e-6;
    end

    % First put Wr in canonical ordering (this checks the input validity too)
    % The spaces P,AI,AO,...,F then correspond to dims 1,2,3,...,2*N+2
    [Wr, dims, parties] = superop_to_canonical_ordering(Wr, dims, parties);
    
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
        is_QCPAR = [trace(W) == d_O, superop_in_QCPAR_cone(Wr,dims,parties,tol)];
    else
        is_QCPAR = (abs(trace(W) - d_O) <= tol) && superop_in_QCPAR_cone(Wr,dims,parties,tol);
    end
end

