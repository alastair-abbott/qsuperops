function Wr = random_superop(dims,parties,R,superop_class)
%random_superop Produces a random superoperator in the specified class
%   W = random_superop(dims,parties,superop_class,R)
%   dims, parties: specifies the desired structure of the superoperator
%   superop_class: a string, one of: QCPAR, QCFO, GEN (default)
%   R: the number of elements of the superinstrument (default: 1)
%
%   The sampling simply generates a random PSD matrix with the correct norm
%   and projects onto the space of valid superoperators. No claims are made about the uniformity
%   of this approach (it is certainly not uniform)
%   For superinstruments, we essentially measure an extra R-dimensional space in F and take the resulting
%   post-selected superoperators as the elements of the superinstrument
%
% Requires QETLAB for RandomDensityMatrix

% Written by Alastair Abbott 2022, last modified 20 September 2022

    if ~exist('superop_class','var')
        superop_class = 'GEN';
    end

    if ~exist('R','var')
        R = 1;
    end

    N = length(parties) - 2;

    % First we add an extra space to F that will split the superoperator into an R-element instrument
    dims_extended = [dims, R];
    d = prod(dims_extended);
    parties_extended = parties;
    parties_extended{N+2}{1} = [parties{N+2}{1}, length(dims_extended)];

    d_O = 1;
    d_O = d_O*prod(dims(parties{1}{1}));
    for n = 1:N
        d_O = d_O*prod(dims(parties{n+1}{2}));
    end

    % Generate a random PSD matrix
    W = d_O*RandomDensityMatrix(d);

    % Project onto valid subspace
    Wr = project_onto_superop_subspace(W,dims_extended,parties_extended,superop_class);

    % Measure register to make R-element superop if R>1
    if R > 1
        basis = eye(R);
        M = zeros(R,R,R);
        for r = 1:R
            M(:,:,r) = basis(:,r)*basis(r,:);
        end
        [Wr,~,~] = measure_superop_output(Wr,dims_extended,parties_extended,M,length(parties_extended{N+2}{1}));
    end

end