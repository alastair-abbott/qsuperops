function [r_opt, yalmip_out] = superop_generalised_robustness(Wr,dims,parties,superop_class,yalmip_options,unitary_ops)
%superop_generalised_robustness Calculates the generalised robustness of a process or superinstrument
%   [r, yalmip_out] = superop_generalised_robustness(Wr,dims,parties,superop_class,yalmip_options,unitary_ops)
%
%   Computes the random robustness of the superinstrument Wr with respect to the given class
%   of superoperators (by default, QC-CCs). 
%   superop_class: a string, on of: QCPAR, QCFO, convQCFO, QCCC, QCQC
%   yalmip_options: Provide settings to be passed to yalmip (e.g., choosing SDP solver)
%   unitary_ops: Calculate robustness corresponding to a witness with unitary operations for
%               operations A_1,...,A_N. False by default.

% Written by Alastair Abbott (2022), last modified 23 August 2022

    %% Process the input
    % First put Wr in canonical ordering (this checks the input validity too)
    % The spaces P,AI,AO,...,F then correspond to dims 1,2,3,...,2*N+2
    if exist('parties','var') && ~isempty(parties)
        [Wr, dims, parties] = superop_to_canonical_ordering(Wr, dims, parties);
    else
        [Wr, dims, parties] = superop_to_canonical_ordering(Wr, dims);
    end

    if ~exist('superop_class','var')
        % By default, we do for QC-CCs
        disp('Calculating the random robustness wrt class QC-CC');
        superop_class = 'QCCC';
    end

    if ~exist('unitary_ops','var')
       unitary_ops = false; 
    end
    
    if ~iscell(Wr)
        Wr = {Wr};
    end
    
    R = length(Wr);
    N = length(parties) - 2;
    
    d = prod(dims);
    
    % Note that trivial spaces are explicitly listed as 1D after canonical ordering
    d_O = prod(dims(1:2:(end-1)));
    
    %% Setup the robustness SDP
    % Noise process
    Omegar = cell(1,R);
    for r = 1:R
        Omegar{r} = sdpvar(d,d,'hermitian','complex');
    end
       
    % Omegar must be a valid superinstrument
    constr = superop_in_valid_cone(Omegar,dims,parties);

    % Robustness is given by norm of Omega
    Omega = Omegar{1};
    for r = 2:R
        Omega = Omega + Omegar{r};
    end
    r_g = real(trace(Omega))/d_O;
    
    % T is used to calculate robustness wrt unitary witnesses
    T = cell(1,R);
    Wr_admixed = cell(1,R);
    for r = 1:R
        if unitary_ops
            T{r} = sdpvar(d,d,'hermitian','complex');
        else
            T{r} = zeros(d,d);
        end
        % The elements of the superinstrment that should be in the QCCC cone
        Wr_admixed{r} = Wr{r} + Omegar{r} - T{r};
    end
    
    switch upper(superop_class)
        case 'QCPAR'
            % disp('Calculating the random robustness wrt class QC-PAR (Parallel quantum circuits)');
            constr = [constr, superop_in_QCPAR_cone(Wr_admixed,dims,parties)];
        case 'QCFO'
            % disp('Calculating the random robustness wrt class QC-FO');
            constr = [constr, superop_in_QCFO_cone(Wr_admixed,dims,parties)];
        case 'CONVQCFO'
            % disp('Calculating the random robustness wrt class conv(QC-FO)');
            constr = [constr, superop_in_convQCFO_cone(Wr_admixed,dims,parties)];
        case 'QCCC'
            % disp('Calculating the random robustness wrt class QC-CC');
            constr = [constr, superop_in_QCCC_cone(Wr_admixed,dims,parties)];
        case 'QCQC'
            % disp('Calculating the random robustness wrt class QC-QC');
            constr = [constr, superop_in_QCQC_cone(Wr_admixed,dims,parties)];
        otherwise
            disp('Warning, invalid superoperator type specified. Calculating for QC-CCs')
            constr = [constr, superop_in_QCCC_cone(Wr_admixed,dims,parties)];
    end
    disp('');
    
    % Enforce unitary witness in primal form
    if unitary_ops
        for i = 1:R
            T_proj = T{i};
            for n = 1:N
                % Constraint is to be in null space of (1-AI)AO and (1-AO)AI, for each party
                T_proj = T_proj - (tr_replace(T_proj,2*n,dims) - tr_replace(T_proj,[2*n,2*n+1],dims));
                T_proj = T_proj - (tr_replace(T_proj,2*n+1,dims) - tr_replace(T_proj,[2*n,2*n+1],dims));
            end
            constr = [constr, T_proj == 0];
        end
    end
    
    %% Solve the SDP
    if ~exist('yalmip_options','var')
        % default options
        yalmip_options = sdpsettings();
    end
    yalmip_out = optimize(constr, r_g, yalmip_options);
    if yalmip_out.problem == true
        disp(['WARNING: Yalmip encountered a problem ("', yalmip_out.info ,'"), result may be unreliable']);
    end
    
    r_opt = value(r_g);
    
end

