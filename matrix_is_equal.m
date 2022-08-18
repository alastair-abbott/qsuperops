function is_equal = matrix_is_equal(A,B,tol)
%matrix_is_equal Determines matrix equality up to a given tolerance
%   is_equal = matrix_is_equal(A,B,tol) is true if max(|A-B|) < tol
%   If tol is omitted, then checks if A == B

% Written by Alastair Abbott, last modified 19 April 2021

    if nargin == 2
		tol = 0;
    end

    is_equal = all(all(abs(A-B) <= tol));
    
end

