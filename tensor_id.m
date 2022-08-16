function AI = tensor_id(A,d,id_at_end)
%tensor_id Efficient tensoring of a matrix with the identity
%   AI = tensor_id(A,d) computes tensor_id(A,eye(d))
%   AI = tensor_id(A,d,id_at_end) computes kron(A,eye(d)) if id_at_end true, otherwise kron(eye(d),A)
%
%   This method is much faster when using overloaded kron to tensor sdp variables with Yalmip or cvx.
%   Here we make use of fact that kron(eye(d), A) is just a block-diagonal matrix with several
%   copies of A as the blocks.
%
%   Requires QETLAB for PermuteSystems

% Written by Alastair Abbott, last modified 19 April 2021

    % By default, we tensor eye(d) at end
    if nargin == 2
        id_at_end = true;
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
    
    if id_at_end
        % We need to swap the eye(d) to the end
        AI = PermuteSystems(AI,[2,1],[d,size(A,1)]);
    end

end

