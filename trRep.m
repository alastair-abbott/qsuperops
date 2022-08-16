function XW = trRep(W, sys, dim)
%trRep Computes the "trace-and-replace" operation
%   XW = trRep(W,sys,dim) Traces out the systems specified by sys
%   and replaces them with the normalised identity
%   See, e.g., arXiv:1506.03776 [quant-ph].
%	
%   dim = [d1, d2, ...] specifies the dimensions of the spaces.
%   sys is a vector indicating which of these spaces to trace out.

% Written by Alastair Abbott, last modified 19 April 2021

 	d = prod(dim); % Size of the full space (i.e., matrix size)
    assert(all(d == size(W)), 'Error: W size doesn''t match supplied dimensions');
    assert(all(sys > 0) && all(sys <= length(dim)), 'Error: Invalid systems specified');
	
    dX = prod(dim(sys)); % Size of spaces to trace and replace
	% If all sys are dim 1 or sys = [], do nothing
    if dX == 1
        XW = W;
        return
    end
	
 	idX = 1/dX*eye(dX);
	% If tracing out all systems, we can just return the identity
    if dX == d
        XW = trace(W)*idX;
        return
    end
	
	% Trace out the systems in sys and tensor idX at the start
    Wtraced = PartialTrace(W, sys, dim);
    % W = Tensor(idX, Wtraced);
    % Faster approach without explicitly computing Tensor (otherwise it's very slow with sdp vars)
    XW = tensorID(Wtraced,dX,false)/dX;
    
    % PermuteSystems doesn't play well when some systems have dimension 1,
    % so we pretend these systems don't exist and adjust the indices accordingly
    if sum(dim(2:end) == 1) > 0 % (if first system has dim 1, it works ok)
        j = 1;
        sysKey = zeros(1,length(dim));
        for i = 1:length(dim)
            if dim(i) > 1
                sysKey(i) = j;
                j = j+1;
            else
                sysKey(i) = 0;
            end
        end
        % Just keep the systems that had dimension > 1
        sys = sysKey(sys);
        sys = sys(sys>0);
        dim = dim(dim>1);
    end
    
    % Get the permutation to put the identities in the right places
    rest = 1:length(dim);
    rest(sys) = [];
    XW = PermuteSystems(XW,[sys rest],dim([sys rest]),0,1);
end

