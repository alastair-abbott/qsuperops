function A_CJ = pureCJ(A)
%pureCJ Calculates the pure Choi-Jamialkowski isomorphism
%   A_CJ = pureCJ(A) returns the pure Choi state |A>> = (id\oplus A)|id>>
	
	d = size(A,2);
	id_CJ = zeros(d^2,1);
	b = eye(d); % basis
	for i = 1:d
		id_CJ = id_CJ + Tensor(b(:,i),b(:,i));
	end
	
	A_CJ = Tensor(eye(d),A)*id_CJ;

end

