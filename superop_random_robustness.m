function [r_opt, yalmip_out] = superop_random_robustness(Wr,dims,parties,superop_class,yalmip_options,unitary_ops)
%superop_random_robustness Calculates the random robustness of a process or superinstrument
%   [r, yalmip_out] = superop_random_robustness(Wr,dims,parties,superop_class,yalmip_options,unitary_ops)
%
%   Computes the random robustness of the superinstrument Wr with respect to the given class
%   of superoperators (by default, QC-CCs). 
%   superop_class: 1 - QC-PAR; 2 - QC-FO; 3 (detault) - QC-CC; 4 - QC-QC
%   yalmip_options: Provide settings to be passed to yalmip (e.g., choosing SDP solver)
%   unitary_ops: Calculate robustness corresponding to a witness with unitary operations for
%               operations A_1,...,A_N. False by default.

% Written by Alastair Abbott (2021), last modified 16 August 2022

    %% Constants
    SUPEROP_CLASS_PAR = 1;
    SUPEROP_CLASS_QCFO = 2;
    SUPEROP_CLASS_QCCC = 3;
    SUPEROP_CLASS_QCQC = 4;

    %% Process the input
    % First put Wr in canonical ordering (this checks the input validity too)
    % The spaces P,AI,AO,...,F then correspond to dims 1,2,3,...,2*N+2
    [Wr, dims, parties] = superop_to_canonical_ordering(Wr, dims, parties);

    if ~exist('superop_class','var')
        % By default, we do for QC-CCs
        superop_class = SUPEROP_CLASS_QCCC;
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
    d_I = prod(dims(2:2:end));
    
    %% Setup the robustness SDP
    % White noise process matrix
    noisy_W = (1/d_I)*eye(d);
       
    r_rand = sdpvar(1);
    constr = [r_rand >= 0];
    
    % T is used to calculate robustness wrt unitary witnesses
    T = cell(1,R);
    Wr_admixed = cell(1,R);
    for i = 1:R
        if unitary_ops
            T{i} = sdpvar(d,d,'hermitian','complex');
        else
            T{i} = zeros(d,d);
        end
        % The elements of the superinstrment that should be in the QCCC cone
        Wr_admixed{i} = Wr{i} + r_rand*noisy_W/R - T{i};
    end
    
    switch superop_class
        case SUPEROP_CLASS_PAR
            disp('Calculating the random robustness wrt class QC-PAR (Parallel quantum circuits)');
            constr = [constr, superop_in_QCPAR_cone(Wr_admixed,dims,parties)];
        case SUPEROP_CLASS_QCFO
            disp('Calculating the random robustness wrt class QC-FO');
            disp('TODO');
        case SUPEROP_CLASS_QCCC
            disp('Calculating the random robustness wrt class QC-CC');
            constr = [constr, superop_in_QCCC_cone(Wr_admixed,dims,parties)];
        case SUPEROP_CLASS_QCQC
            disp('Calculating the random robustness wrt class QC-QC');
            disp('TODO');
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
    yalmip_out = optimize(constr, r_rand, yalmip_options);
    
    r_opt = value(r_rand);
    
end

