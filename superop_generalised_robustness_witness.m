function [Sr_opt, yalmip_out, coeffs_opt] = superop_generalised_robustness_witness(Wr,dims_raw,parties_raw,superop_class,yalmip_options,unitary_ops,witness_basis)
%superop_generalised_robustness_witness Calculates the generalised robustness witness of a process or superinstrument
%   [Sr, yalmip_out] = superop_generalised_robustness_witness(Wr,dims,parties,superop_class,yalmip_options,unitary_ops)
%   [Sr, yalmip_out, coeffs_opt] = superop_generalised_robustness_witness(Wr,dims,parties,superop_class,yalmip_options,unitary_ops,witness_basis)
%
%   Compute the witness (wrt generalised robustness) of the superinstrument Wr with respect to given class of superoperators.  
%   superop_class: a string, on of: QCPAR, QCFO, convQCFO, QCCC, QCQC
%
%   yalmip_options: Provide settings to be passed to yalmip (e.g., choosing SDP solver)
%   unitary_ops: Calculate robustness corresponding to a witness with unitary operations for
%               operations A_1,...,A_N. False by default.
%   witness_basis: Optional input providing a basis for the witness to be expressed in.
%   In this case, the output coeffs will provide the optimal coeffs giving SrOpt{i} = \sum_b coeffs(i,b)*witnessBase(:,:,b)

% Written by Alastair Abbott (2022), last modified 23 August 2022

    %% Process the input
    % First put Wr in canonical ordering (this checks the input validity too)
    % The spaces P,AI,AO,...,F then correspond to dims 1,2,3,...,2*N+2
    if exist('parties_raw','var') && ~isempty(parties_raw)
        [Wr, dims, parties] = superop_to_canonical_ordering(Wr, dims_raw, parties_raw);
    else
        [Wr, dims, parties] = superop_to_canonical_ordering(Wr, dims_raw);
        parties_raw = parties;
    end
    
    if ~exist('superop_class','var')
        % By default, we do for QC-CCs
        disp('Calculating the random robustness witness wrt class QC-CC');
        superop_class = 'QCCC';
    end

    if ~exist('unitary_ops','var')
       unitary_ops = false; 
    end
    
    witness_from_basis = false;
    if exist('witness_basis','var')
        witness_from_basis = true;
    end

    % Record if input was just a process matrix so we can return a simple witness rather than a cell
    input_is_process_matrix = false;
    if ~iscell(Wr)
        input_is_process_matrix = true;
        Wr = {Wr};
    end
    
    R = length(Wr);
    N = length(parties) - 2;
    
    d = prod(dims);
    
    % Note that trivial spaces are explicitly listed as 1D after canonical ordering
    d_I = prod(dims(2:2:end));
    
    %% Setup the witness SDP
    % White noise process matrix
    noisyW = (1/d_I)*eye(d);

    if witness_from_basis
        % We decompose Sr{i} = \sum_b coeffs(i,b)*witnessBase(:,:,b)
        basisSize = size(witness_basis,3); % num elements in basis
        coeffs = sdpvar(R,basisSize,'full');
    end
    
    Sr = cell(1,R);
    for i = 1:R
        if witness_from_basis
            Sr{i} = zeros(d,d);
            for b = 1:basisSize
               Sr{i} = Sr{i} + coeffs(i,b)*witness_basis(:,:,b); 
            end
        else
            Sr{i} =  sdpvar(d,d,'hermitian','complex');
        end
        
    end
    
    switch upper(superop_class)
        case 'QCPAR'
            % disp('Calculating the random robustness witness wrt class QC-PAR (Parallel quantum circuits)');
            disp('TODO');
        case 'QCFO'
            % disp('Calculating the random robustness witness wrt class QC-FO');
            constr = superop_in_QCFO_witness_cone(Sr,dims,parties);
        case 'CONVQCFO'
            % disp('Calculating the random robustness witness wrt class conv(QC-FO)');
            disp('TODO');
        case 'QCCC'
            % disp('Calculating the random robustness witness wrt class QC-CC');
            constr = superop_in_QCCC_witness_cone(Sr,dims,parties);
        case 'QCQC'
            % disp('Calculating the random robustness witness wrt class QC-QC');
            disp('TODO');
        otherwise
            disp('Warning, invalid superoperator type specified. Calculating for QC-CCs')
            constr = superop_in_QCCC_witness_cone(Sr,dims,parties);
    end
    disp('');
    
    
    % Enforce witness to have unitary operations
    if unitary_ops
        if witness_from_basis
           disp('Warning: shouldn''t need to force unitary operations when a good witness basis is provided.'); 
        end
        for i = 1:R
            for n = 1:N
                constr = [constr, tr_replace(Sr{i},2*n,dims) == tr_replace(Sr{i},[2*n,2*n+1],dims)];
                constr = [constr, tr_replace(Sr{i},2*n+1,dims) == tr_replace(Sr{i},[2*n,2*n+1],dims)];
            end
        end
    end
    
    % normalisation condition
    norm = 0;
    for i = 1:R
        norm = norm + 1/R*trace(Sr{i}*noisyW);
    end
    constr = [constr, norm == 1];
    
    % objective: tr(S*W)
    objective = 0;
    for i = 1:R
        objective = objective + trace(Sr{i}*Wr{i});
    end
    
    %% Solve the SDP

    if ~exist('yalmip_options','var')
        % default options
        yalmip_options = sdpsettings();
    end
    yalmip_out = optimize(constr, real(objective), yalmip_options);
    
    if input_is_process_matrix
        Sr_opt = value(Sr{1});
    else
        Sr_opt = cell(1,R);
        for i = 1:R
           Sr_opt{i} = value(Sr{i}); 
        end
    end    
    % Need to put back in the non-canonical ordering of the input
    Sr_opt = superop_from_canonical_ordering(Sr_opt,dims_raw,parties_raw);
    
    if witness_from_basis
        coeffs_opt = value(coeffs);
    else
        coeffs_opt = [];
    end
end

