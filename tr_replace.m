function XW = tr_replace(W, sys, dims)
%tr_replace Computes the "trace-and-replace" operation
%   XW = tr_replace(W,sys,dims) Traces out the systems specified by sys
%   and replaces them with the normalised identity
%   See, e.g., arXiv:1506.03776 [quant-ph].
%	
%   dims = [d1, d2, ...] specifies the dimensions of the spaces.
%   sys is a vector indicating which of these spaces to trace out.
%
% Requires QETLAB for PartialTrace and PermuteSystems

% Written by Alastair Abbott, last modified 19 April 2021

 	d = prod(dims); % Size of the full space (i.e., matrix size)
    assert(all(d == size(W)), 'Error: W size doesn''t match supplied dimensions');
    assert(all(sys > 0) && all(sys <= length(dims)), 'Error: Invalid systems specified');
	
    dX = prod(dims(sys)); % Size of spaces to trace and replace
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
    Wtraced = PartialTrace(W, sys, dims);
    % W = kron(idX, Wtraced);
    % Faster approach without explicitly computing kron (otherwise it's very slow with sdp vars)
    XW = tensor_id(Wtraced,dX,false)/dX;
    
    % PermuteSystems doesn't play well when some systems have dimension 1,
    % so we pretend these systems don't exist and adjust the indices accordingly
    if sum(dims(2:end) == 1) > 0 % (if first system has dim 1, it works ok)
        j = 1;
        sysKey = zeros(1,length(dims));
        for i = 1:length(dims)
            if dims(i) > 1
                sysKey(i) = j;
                j = j+1;
            else
                sysKey(i) = 0;
            end
        end
        % Just keep the systems that had dimension > 1
        sys = sysKey(sys);
        sys = sys(sys>0);
        dims = dims(dims>1);
    end
    
    % Get the permutation to put the identities in the right places
    rest = 1:length(dims);
    rest(sys) = [];
    XW = PermuteSystems(XW,[sys rest],dims([sys rest]),0,1);
end

