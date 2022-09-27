function Wr = random_superop(dims,parties,R,sampling_method)
%random_superop Produces a random superoperator in the specified class
%   W = random_superop(dims,parties,superop_class,R)
%   dims, parties: specifies the desired structure of the superoperator
%   R: the number of elements of the superinstrument (default: 1)
%   sampling_method: the method used to sample random processes (default: 'pure_projection')
%
%   Several sampling methods are possible:
%   'pure_projection': a random projector is generated and scaled to have the correct norm. This is then
%                      projected on the space of valid superoperators and then mixed with white noise
%                      so that the resulting operator is PSD. This is the default option.
%   'PSD_projection': as for 'pure_projection', except initially a random PSD matrix is generated.
%   'closest_valid': A random projector is generated and then an SDP is used to find the closest
%                    valid superoperator wrt XXX nrom
%   Note: no claims about the uniformity of these methods are made (none of them are uniform)
%   For superinstruments, we essentially measure an extra R-dimensional space in F and take the resulting
%   post-selected superoperators as the elements of the superinstrument
%
% Requires QETLAB for RandomDensityMatrix

% Written by Alastair Abbott 2022, last modified 27 September 2022

    if ~exist('R','var')
        R = 1;
    end

    if ~exist('sampling_method','var')
        sampling_method = 'pure_projection';
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
    d_I = d/d_O;

    % Generate a random operator
    switch lower(sampling_method)
        case 'pure_projection'
            W = d_O*pure_to_mixed(RandomStateVector(d));
        case 'psd_projection'
            W = d_O*RandomDensityMatrix(d);
        otherwise
            error('Error: unknown sampling method specified');
    end

    % Project onto valid subspace
    Wr = project_onto_valid_superops(W,dims_extended,parties_extended);

    % Resulting process may not be PSD; add noise to make it so
    eig_min = min(eig(Wr));
    q = eig_min*d_I/(eig_min*d_I-1);
    if eig_min < 0
        noisyW = eye(d)/d_I;
        Wr = q*noisyW + (1-q)*Wr;
    end

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