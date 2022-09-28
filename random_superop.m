function Wr = random_superop(dims,parties,R,sampling_method,yalmip_options)
%random_superop Produces a random superoperator in the specified class
%   W = random_superop(dims,parties,R,sampling_method,yalmip_options)
%   dims, parties: specifies the desired structure of the superoperator
%   R: the number of elements of the superinstrument (default: 1)
%   sampling_method: the method used to sample random processes (default: 'pure_projection')
%   yalmip_options: Yalmip options relevant for closest_norm methods (see below)
%
%   Several sampling methods are possible:
%   'pure_projection': a random projector is generated and scaled to have the correct norm. This is then
%                      projected on the space of valid superoperators and then mixed with white noise
%                      so that the resulting operator is PSD. This is the default option.
%   'PSD_projection': as for 'pure_projection', except initially a random PSD matrix is generated.
%   'closest_norm_X_pure': A random projector is generated and then an SDP is used to find the closest
%                        valid superoperator wrt the X norm, where X is one of ['1','2','inf','fro','nuclear','*']
%   'closest_norm_X_PSD': As above, but starting from a random PSD matrix
%   'closest_norm_pure','closest_norm_PSD': As above, but the default norm choice of '2' is used
%
%   Note: no claims about the uniformity of these methods are made (none of them are uniform)
%   For superinstruments, we essentially measure an extra R-dimensional space in F and take the resulting
%   post-selected superoperators as the elements of the superinstrument
%
% Requires QETLAB for RandomDensityMatrix

% Written by Alastair Abbott 2022, last modified 27 September 2022

    if ~exist('R','var')
        R = 1;
    end

    valid_sampling_methods = {'pure_projection','psd_projection','closest_norm_pure','closest_norm_PSD',...
        'closest_norm_1_pure','closest_norm_2_pure','closest_norm_inf_pure','closest_norm_fro_pure','closest_norm_nuclear_pure','closest_norm_*_pure',...
        'closest_norm_1_psd','closest_norm_2_psd','closest_norm_inf_psd','closest_norm_fro_psd','closest_norm_nuclear_psd','closest_norm_*_psd'};
    
    if ~exist('sampling_method','var')
        sampling_method = 'pure_projection';
    end
    sampling_method = lower(sampling_method);
    assert(contains(sampling_method,valid_sampling_methods),'Error: Unknown sampling method specified');

    % Setup sampling type information
    pure_seed = false;
    if contains(sampling_method,'pure')
        pure_seed = true;
    end

    projection_sampling = false;
    if contains(sampling_method,'projection')
        projection_sampling = true;
    else
        % extract the right norm
        norm_type = extractBetween(sampling_method,'norm_','_');
        if isempty(norm_type)
            norm_type = 2;
        else
            norm_type = norm_type{1};
            % if norm is 1,2,inf, need to convert to a numbert
            if ~isempty(str2num(norm_type))
                norm_type = str2num(norm_type);
            end
        end
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
    if pure_seed
        W = d_O*pure_to_mixed(RandomStateVector(d));
    else
        W = d_O*RandomDensityMatrix(d);
    end

    if projection_sampling
        % Project onto valid subspace
        Wr = project_onto_valid_superops(W,dims_extended,parties_extended);
    
        % Resulting process may not be PSD; add noise to make it so
        eig_min = min(eig(Wr));
        q = eig_min*d_I/(eig_min*d_I-1);
        if eig_min < 0
            noisyW = eye(d)/d_I;
            Wr = q*noisyW + (1-q)*Wr;
        end
    else
        % Otherwise we have an sdp to solve
        Wr = sdpvar(d,d,'hermitian','complex');
        constr = is_valid_superop(Wr,dims_extended,parties_extended);
        if ~exist('yalmip_options','var')
            % default options
            yalmip_options = sdpsettings();
        end
        yalmip_out = optimize(constr, norm(W-Wr,norm_type), yalmip_options);
        Wr = value(Wr);
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