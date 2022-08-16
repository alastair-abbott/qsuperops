function AI = tensorID(A,d,idAtEnd)
%tensorEye Efficient tensoring of a matrix with the identity
%   AI = tensorEye(A,d) computes Tensor(A,eye(d))
%   AI = tensorEye(A,d,idAtEnd) computes Tensor(A,eye(d)) if idAtEnd true, otherwise Tensor(eye(d),A)
%
%   This method is much faster when using overloaded kron to tensor sdp variables with Yalmip or cvx.
%   Here we make use of fact that Tensor(eye(d), A) is just a block-diagonal matrix with several
%   copies of A as the blocks.

% Written by Alastair Abbott, last modified 19 April 2021

    % By default, we tensor eye(d) at end
    if nargin == 2
        idAtEnd = true;
    end

    switch d
        case 1
            AI = A;
        case 2 % d = 2 case faster
            AI = blkdiag(A,A);
        otherwise
            AI = repmat({A},1,d);
            AI = blkdiag(AI{:});
    end
    
    if idAtEnd
        % We need to swap the eye(d) to the end
        AI = PermuteSystems(AI,[2,1],[d,size(A,1)]);
    end

end

