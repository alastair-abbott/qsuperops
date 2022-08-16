function A_CJ = pure_CJ(A)
%pure_CJ Calculates the pure Choi-Jamialkowski isomorphism
%   A_CJ = pure_CJ(A) returns the pure Choi state |A>> = (id\oplus A)|id>>
	
	d = size(A,2);
	id_CJ = zeros(d^2,1);
	b = eye(d); % basis
	for i = 1:d
		id_CJ = id_CJ + kron(b(:,i),b(:,i));
	end
	
	A_CJ = kron(eye(d),A)*id_CJ;

end

