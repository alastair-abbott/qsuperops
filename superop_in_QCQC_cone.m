function cone_constraints = superop_in_QCQC_cone(Wr, dims, parties)
%superop_in_QCQC_cone Yalmip constraints for a superinstrument to be in the cone of QC-QCs
%   cone_constraints = superop_in_QCQC_cone(Wr, dims, parties)
%   Wr can be either a superinstrument or a superoperator/process matrix
%   Returns the yalmip SDP constraints for Wr to be in the cone of valid Quantum Circuits with Quantum Control of Causal order
%   Note: this doesn't enforce the normalisation of the superoperator
%   
%   Formulation based on: J. Wechs, H. Dourdent, A. A. Abbott, C. Branciard, PRX Quantum 2, 030335 (2021)
%   (See Propositions 7 and 15 therein.)
%
% Requires QETLAB for PermuteSystems, PartialTrace

% Written by Alastair Abbott 2022, last modified 22 August 2022

    % default tolerance
    if ~exist('tol','var')
        tol = 1e-12;
    end

    % First put Wr in canonical ordering (this checks the input validity too)
    % The spaces P,AI,AO,...,F then correspond to dims 1,2,3,...,2*N+2
    if exist('parties','var') && ~isempty(parties)
        [Wr, dims, parties] = superop_to_canonical_ordering(Wr, dims, parties);
    else
        [Wr, dims, parties] = superop_to_canonical_ordering(Wr, dims);
    end

    % Treat a process matrix as a 1-element superinstrument
    if ~iscell(Wr)
        Wr = {Wr};
    end
    
    R = length(Wr);
    N = length(parties) - 2;
    
    d = prod(dims);
    
    %% We treat each N separately
       
    % Need every Wr{i} >= 0
    cone_constraints = superop_in_PSD_cone(Wr);

    W = Wr{1};
    for r = 2:R
       W = W + Wr{r}; 
    end
    
    % We adopt the notation W_{K_{n-1}}_{k_n} to mimic the notation of the paper W_{(K_{n-1},k_n)}
    switch N
        case 1
            P = 1; d_P = dims(P);
            AI = 2;
            AO = 3; d_AO = dims(AO);
            F = 4;
            
            W_0_A = 1/d_AO*PartialTrace(W,[AO,F],dims);

            cone_constraints = [cone_constraints, PartialTrace(W,F,dims) == tensor_id(W_0_A,d_AO)];

            if d_P ~= 1 % Otherwise this is trvial
                % We only want it to be proportional to the identity, since we are only checking the cone
                % We need to also be careful because sometimes this constraint is trivially satisfied
                % but numerical error makes Yalmip think it's false
                diff = PartialTrace(W_0_A,2,dims([P,AI])) - trace(W_0_A)/d_P*eye(d_P);
                if isa(diff,'sdpvar')
                    cone_constraints = [cone_constraints, diff == 0];
                else
                    cone_constraints = [cone_constraints, matrix_is_equal(diff,zeros(d_P),tol)];
                end
            end 

        case 2
            % Recall we now have a conically ordered process and we've grouped the "sub"-spaces
            P = 1; d_P = dims(P);
            AI = 2; 
            AO = 3; d_AO = dims(AO);
            BI = 4; 
            BO = 5; d_BO = dims(BO);
            F = 6;
      
            d_W_A_B = prod(dims([P,AI,AO,BI]));
		    d_W_B_A = prod(dims([P,AI,BI,BO]));

            W_A_B = sdpvar(d_W_A_B,d_W_A_B,'hermitian','complex');
		    W_B_A = sdpvar(d_W_B_A,d_W_B_A,'hermitian','complex');
		    
		    W_0_A = 1/d_AO*PartialTrace(W_A_B,[3,4],dims([P,AI,AO,BI]));
		    W_0_B = 1/d_BO*PartialTrace(W_B_A,[2,4],dims([P,AI,BI,BO]));
		    
		    cone_constraints = [cone_constraints, W_A_B >= 0, W_B_A >= 0];
		    
		    cone_constraints = [cone_constraints, PartialTrace(W,F,dims) == tensor_id(W_A_B,d_BO) + PermuteSystems(tensor_id(W_B_A,d_AO),[P,AI,BI,BO,AO],dims([P,AI,BI,BO,AO]),0,1)];
		    cone_constraints = [cone_constraints, PartialTrace(W_A_B,4,dims([P,AI,AO,BI])) == tensor_id(W_0_A,d_AO)];
		    cone_constraints = [cone_constraints, PartialTrace(W_B_A,2,dims([P,AI,BI,BO])) == tensor_id(W_0_B,d_BO)];

            if d_P ~= 1 % Otherwise this is trvial
                % We only want it to be proportional to the identity, since we are only checking the cone
                % We need to also be careful because sometimes this constraint is trivially satisfied
                % but numerical error makes Yalmip think it's false
                diff = PartialTrace(W_0_A,2,dims([P,AI])) + PartialTrace(W_0_B,2,dims([P,BI])) - trace(W_0_A+W_0_B)/d_P*eye(d_P);
                if isa(diff,'sdpvar')
                    cone_constraints = [cone_constraints, diff == 0];
                else
                    cone_constraints = [cone_constraints, matrix_is_equal(diff,zeros(d_P),tol)];
                end
            end 

        case 3
            P = 1; d_P = dims(P);
            AI = 2; 
            AO = 3; d_AO = dims(AO);
            BI = 4; 
            BO = 5; d_BO = dims(BO);
            CI = 6;
            CO = 7; d_CO = dims(CO);
            F = 8;
            
            d_W_AB_C = prod(dims([P,AI,AO,BI,BO,CI]));
		    d_W_AC_B = prod(dims([P,AI,AO,BI,CI,CO]));
		    d_W_BC_A = prod(dims([P,AI,BI,BO,CI,CO]));

		    W_AB_C = sdpvar(d_W_AB_C,d_W_AB_C,'hermitian','complex');
            W_AC_B = sdpvar(d_W_AC_B,d_W_AC_B,'hermitian','complex');
            W_BC_A = sdpvar(d_W_BC_A,d_W_BC_A,'hermitian','complex');
    
		    d_W_A_B = prod(dims([P,AI,AO,BI]));
		    d_W_A_C = prod(dims([P,AI,AO,CI]));
		    d_W_B_A = prod(dims([P,AI,BI,BO]));
		    d_W_B_C = prod(dims([P,BI,BO,CI]));
		    d_W_C_A = prod(dims([P,AI,CI,CO]));
		    d_W_C_B = prod(dims([P,BI,CI,CO]));

		    W_A_B = sdpvar(d_W_A_B,d_W_A_B,'hermitian','complex');
            W_A_C = sdpvar(d_W_A_C,d_W_A_C,'hermitian','complex');
            W_B_A = sdpvar(d_W_B_A,d_W_B_A,'hermitian','complex');
            W_B_C = sdpvar(d_W_B_C,d_W_B_C,'hermitian','complex');
            W_C_A = sdpvar(d_W_C_A,d_W_C_A,'hermitian','complex');
            W_C_B = sdpvar(d_W_C_B,d_W_C_B,'hermitian','complex');
		    
		    W_0_A = 1/d_AO*(PartialTrace(W_A_B,[3,4],dims([P,AI,AO,BI])) + PartialTrace(W_A_C,[3,4],dims([P,AI,AO,CI])));
		    W_0_B = 1/d_BO*(PartialTrace(W_B_A,[2,4],dims([P,AI,BI,BO])) + PartialTrace(W_B_C,[3,4],dims([P,BI,BO,CI])));
		    W_0_C = 1/d_CO*(PartialTrace(W_C_A,[2,4],dims([P,AI,CI,CO])) + PartialTrace(W_C_B,[2,4],dims([P,BI,CI,CO])));
		    
		    cone_constraints = [cone_constraints, W_AB_C >= 0, W_AC_B >= 0, W_BC_A >= 0, W_A_B >= 0, W_A_C >= 0, W_B_A >= 0, W_B_C >= 0, W_C_A >= 0, W_C_B >= 0];		
		    
		    cone_constraints = [cone_constraints, PartialTrace(W,F,dims) == tensor_id(W_AB_C,d_CO) + PermuteSystems(tensor_id(W_AC_B,d_BO),[P,AI,AO,BI,CI,CO,BO],dims([P,AI,AO,BI,CI,CO,BO]),0,1) + PermuteSystems(tensor_id(W_BC_A,d_AO),[P,AI,BI,BO,CI,CO,AO],dims([P,AI,BI,BO,CI,CO,AO]),0,1)];
	    
		    cone_constraints = [cone_constraints, PartialTrace(W_AB_C,6,dims([P,AI,AO,BI,BO,CI])) == tensor_id(W_A_B,d_BO) + PermuteSystems(tensor_id(W_B_A,d_AO),[1,2,4,5,3],dims([P,AI,BI,BO,AO]),0,1)];
		    cone_constraints = [cone_constraints, PartialTrace(W_AC_B,4,dims([P,AI,AO,BI,CI,CO])) == tensor_id(W_A_C,d_CO) + PermuteSystems(tensor_id(W_C_A,d_AO),[1,2,4,5,3],dims([P,AI,CI,CO,AO]),0,1)];
		    cone_constraints = [cone_constraints, PartialTrace(W_BC_A,2,dims([P,AI,BI,BO,CI,CO])) == tensor_id(W_B_C,d_CO) + PermuteSystems(tensor_id(W_C_B,d_BO),[1,2,4,5,3],dims([P,BI,CI,CO,BO]),0,1)];
		    
		    cone_constraints = [cone_constraints, PartialTrace(W_A_B,4,dims([P,AI,AO,BI])) + PartialTrace(W_A_C,4,dims([P,AI,AO,CI])) == tensor_id(W_0_A,d_AO)];
		    cone_constraints = [cone_constraints, PartialTrace(W_B_A,2,dims([P,AI,BI,BO])) + PartialTrace(W_B_C,4,dims([P,BI,BO,CI])) == tensor_id(W_0_B,d_BO)];
		    cone_constraints = [cone_constraints, PartialTrace(W_C_A,2,dims([P,AI,CI,CO])) + PartialTrace(W_C_B,2,dims([P,BI,CI,CO])) == tensor_id(W_0_C,d_CO)];
            
            if d_P ~= 1 %Otherwise this is trvial
                % We only want it to be proportional to the identity, since we are only checking the cone
                % We need to also be careful because sometimes this constraint is trivially satisfied
                % but numerical error makes Yalmip think it's false
                diff = PartialTrace(W_0_A,2,dims([P,AI])) + PartialTrace(W_0_B,2,dims([P,BI])) + PartialTrace(W_0_C,2,dims([P,CI])) - trace(W_0_A+W_0_B+W_0_C)/d_P*eye(d_P);
                if isa(diff,'sdpvar')
                    cone_constraints = [cone_constraints, diff == 0];
                else
                    cone_constraints = [cone_constraints, matrix_is_equal(diff,zeros(d_P),tol)];
                end
            end     

        case 4
            P = 1; d_P = dims(P);
            AI = 2; 
            AO = 3; d_AO = dims(AO);
            BI = 4; 
            BO = 5; d_BO = dims(BO);
            CI = 6;
            CO = 7; d_CO = dims(CO);
            DI = 8;
            DO = 9; d_DO = dims(DO);
            F = 10;
            
            d_W_ABC_D = prod(dims([P,AI,AO,BI,BO,CI,CO,DI]));
		    d_W_ABD_C = prod(dims([P,AI,AO,BI,BO,CI,DI,DO]));
		    d_W_ACD_B = prod(dims([P,AI,AO,BI,CI,CO,DI,DO]));
		    d_W_BCD_A = prod(dims([P,AI,BI,BO,CI,CO,DI,DO]));

		    W_ABC_D = sdpvar(d_W_ABC_D,d_W_ABC_D,'hermitian','complex');
            W_ABD_C = sdpvar(d_W_ABD_C,d_W_ABD_C,'hermitian','complex');
            W_ACD_B = sdpvar(d_W_ACD_B,d_W_ACD_B,'hermitian','complex');
            W_BCD_A = sdpvar(d_W_BCD_A,d_W_BCD_A,'hermitian','complex');
		    
		    d_W_AB_C = prod(dims([P,AI,AO,BI,BO,CI]));
		    d_W_AB_D = prod(dims([P,AI,AO,BI,BO,DI]));
		    d_W_AC_B = prod(dims([P,AI,AO,BI,CI,CO]));
		    d_W_AC_D = prod(dims([P,AI,AO,CI,CO,DO]));    
		    d_W_AD_B = prod(dims([P,AI,AO,BI,DI,DO]));	    
		    d_W_AD_C = prod(dims([P,AI,AO,CI,DI,DO]));
		    d_W_BC_A = prod(dims([P,AI,BI,BO,CI,CO]));
		    d_W_BC_D = prod(dims([P,BI,BO,CI,CO,DI]));
		    d_W_BD_A = prod(dims([P,AI,BI,BO,DI,DO]));
		    d_W_BD_C = prod(dims([P,BI,BO,CI,DI,DO]));
		    d_W_CD_A = prod(dims([P,AI,CI,CO,DI,DO]));
		    d_W_CD_B = prod(dims([P,BI,CI,CO,DI,DO]));

		    W_AB_C = sdpvar(d_W_AB_C,d_W_AB_C,'hermitian','complex');
		    W_AB_D = sdpvar(d_W_AB_D,d_W_AB_D,'hermitian','complex');
		    W_AC_B = sdpvar(d_W_AC_B,d_W_AC_B,'hermitian','complex');
            W_AC_D = sdpvar(d_W_AC_D,d_W_AC_D,'hermitian','complex');
            W_AD_B = sdpvar(d_W_AD_B,d_W_AD_B,'hermitian','complex');
		    W_AD_C = sdpvar(d_W_AD_C,d_W_AD_C,'hermitian','complex');
		    W_BC_A = sdpvar(d_W_BC_A,d_W_BC_A,'hermitian','complex');
		    W_BC_D = sdpvar(d_W_BC_D,d_W_BC_D,'hermitian','complex');
            W_BD_A = sdpvar(d_W_BD_A,d_W_BD_A,'hermitian','complex');
            W_BD_C = sdpvar(d_W_BD_C,d_W_BD_C,'hermitian','complex');
            W_CD_A = sdpvar(d_W_CD_A,d_W_CD_A,'hermitian','complex');
		    W_CD_B = sdpvar(d_W_CD_B,d_W_CD_B,'hermitian','complex');
    
		    d_W_A_B = prod(dims([P,AI,AO,BI]));
		    d_W_A_C = prod(dims([P,AI,AO,CI]));
            d_W_A_D = prod(dims([P,AI,AO,DI]));
            d_W_B_A = prod(dims([P,AI,BI,BO]));
            d_W_B_C = prod(dims([P,BI,BO,CI]));
            d_W_B_D = prod(dims([P,BI,BO,DI]));
            d_W_C_A = prod(dims([P,AI,CI,CO]));
            d_W_C_B = prod(dims([P,BI,CI,CO]));
            d_W_C_D = prod(dims([P,CI,CO,DI]));
            d_W_D_A = prod(dims([P,AI,DI,DO]));
            d_W_D_B = prod(dims([P,BI,DI,DO]));
            d_W_D_C = prod(dims([P,CI,DI,DO]));

            W_A_B = sdpvar(d_W_A_B,d_W_A_B,'hermitian','complex');
		    W_A_C = sdpvar(d_W_A_C,d_W_A_C,'hermitian','complex');
		    W_A_D = sdpvar(d_W_A_D,d_W_A_D,'hermitian','complex');
		    W_B_A = sdpvar(d_W_B_A,d_W_B_A,'hermitian','complex');
		    W_B_C = sdpvar(d_W_B_C,d_W_B_C,'hermitian','complex');
		    W_B_D = sdpvar(d_W_B_D,d_W_B_D,'hermitian','complex');
		    W_C_A = sdpvar(d_W_C_A,d_W_C_A,'hermitian','complex');
		    W_C_B = sdpvar(d_W_C_B,d_W_C_B,'hermitian','complex');
		    W_C_D = sdpvar(d_W_C_D,d_W_C_D,'hermitian','complex');
		    W_D_A = sdpvar(d_W_D_A,d_W_D_A,'hermitian','complex');
		    W_D_B = sdpvar(d_W_D_B,d_W_D_B,'hermitian','complex');
		    W_D_C = sdpvar(d_W_D_C,d_W_D_C,'hermitian','complex');
		    
		    W_0_A = 1/d_AO*(PartialTrace(W_A_B,[3,4],dims([P,AI,AO,BI])) + PartialTrace(W_A_C,[3,4],dims([P,AI,AO,CI])) + PartialTrace(W_A_D,[3,4],dims([P,AI,AO,DI])));
		    W_0_B = 1/d_BO*(PartialTrace(W_B_A,[2,4],dims([P,AI,BI,BO])) + PartialTrace(W_B_C,[3,4],dims([P,BI,BO,CI])) + PartialTrace(W_B_D,[3,4],dims([P,BI,BO,DI])));
		    W_0_C = 1/d_CO*(PartialTrace(W_C_A,[2,4],dims([P,AI,CI,CO])) + PartialTrace(W_C_B,[2,4],dims([P,BI,CI,CO])) + PartialTrace(W_C_D,[3,4],dims([P,CI,CO,DI])));
		    W_0_D = 1/d_DO*(PartialTrace(W_D_A,[2,4],dims([P,AI,DI,DO])) + PartialTrace(W_D_B,[2,4],dims([P,BI,DI,DO])) + PartialTrace(W_D_C,[2,4],dims([P,CI,DI,DO])));
		    
		    cone_constraints = [cone_constraints, W_ABC_D >= 0, W_ABD_C >=0, W_ACD_B >=0, W_BCD_A >= 0, ...
			                    W_AB_C >= 0, W_AB_D >= 0, W_AC_B >= 0, W_AC_D >= 0, W_AD_B >= 0, W_AD_C >= 0, W_BC_A >= 0, W_BC_D >= 0, W_BD_A >= 0, W_BD_C >= 0, W_CD_A >= 0, W_CD_B >= 0 ...
			                    W_A_B >= 0, W_A_C >= 0, W_A_D >= 0, W_B_A >= 0, W_B_C >= 0, W_B_D >= 0, W_C_A >= 0, W_C_B >= 0, W_C_D >= 0, W_D_A >= 0, W_D_B >= 0, W_D_C >= 0];		
		    
		    cone_constraints = [cone_constraints, PartialTrace(W,F,dims) == tensor_id(W_ABC_D,d_DO) + PermuteSystems(tensor_id(W_ABD_C,d_CO),[P,AI,AO,BI,BO,CI,DI,DO,CO],dims([P,AI,AO,BI,BO,CI,DI,DO,CO]),0,1) + PermuteSystems(tensor_id(W_ACD_B,d_BO),[P,AI,AO,BI,CI,CO,DI,DO,BO],dims([P,AI,AO,BI,CI,CO,DI,DO,BO]),0,1) + PermuteSystems(tensor_id(W_BCD_A,d_AO),[P,AI,BI,BO,CI,CO,DI,DO,AO],dims([P,AI,BI,BO,CI,CO,DI,DO,AO]),0,1)];
	    
		    cone_constraints = [cone_constraints, PartialTrace(W_ABC_D,8,dims([P,AI,AO,BI,BO,CI,CO,DI])) == tensor_id(W_AB_C,d_CO) + PermuteSystems(tensor_id(W_AC_B,d_BO),[1,2,3,4,6,7,5],dims([P,AI,AO,BI,CI,CO,BO]),0,1) + PermuteSystems(tensor_id(W_BC_A,d_AO),[1,2,4,5,6,7,3],dims([P,AI,BI,BO,CI,CO,AO]),0,1)];
		    cone_constraints = [cone_constraints, PartialTrace(W_ABD_C,6,dims([P,AI,AO,BI,BO,CI,DI,DO])) == tensor_id(W_AB_D,d_DO) + PermuteSystems(tensor_id(W_AD_B,d_BO),[1,2,3,4,6,7,5],dims([P,AI,AO,BI,DI,DO,BO]),0,1) + PermuteSystems(tensor_id(W_BD_A,d_AO),[1,2,4,5,6,7,3],dims([P,AI,BI,BO,DI,DO,AO]),0,1)];
		    cone_constraints = [cone_constraints, PartialTrace(W_ACD_B,4,dims([P,AI,AO,BI,CI,CO,DI,DO])) == tensor_id(W_AC_D,d_DO) + PermuteSystems(tensor_id(W_AD_C,d_CO),[1,2,3,4,6,7,5],dims([P,AI,AO,CI,DI,DO,CO]),0,1) + PermuteSystems(tensor_id(W_CD_A,d_AO),[1,2,4,5,6,7,3],dims([P,AI,CI,CO,DI,DO,AO]),0,1)];
		    cone_constraints = [cone_constraints, PartialTrace(W_BCD_A,2,dims([P,AI,BI,BO,CI,CO,DI,DO])) == tensor_id(W_BC_D,d_DO) + PermuteSystems(tensor_id(W_BD_C,d_CO),[1,2,3,4,6,7,5],dims([P,BI,BO,CI,DI,DO,CO]),0,1) + PermuteSystems(tensor_id(W_CD_B,d_BO),[1,2,4,5,6,7,3],dims([P,BI,CI,CO,DI,DO,BO]),0,1)];
		    
		    cone_constraints = [cone_constraints, PartialTrace(W_AB_C,6,dims([P,AI,AO,BI,BO,CI])) + PartialTrace(W_AB_D,6,dims([P,AI,AO,BI,BO,DI])) == tensor_id(W_A_B,d_BO) + PermuteSystems(tensor_id(W_B_A,d_AO),[1,2,4,5,3],dims([P,AI,BI,BO,AO]),0,1)];
		    cone_constraints = [cone_constraints, PartialTrace(W_AC_B,4,dims([P,AI,AO,BI,CI,CO])) + PartialTrace(W_AC_D,6,dims([P,AI,AO,CI,CO,DI])) == tensor_id(W_A_C,d_CO) + PermuteSystems(tensor_id(W_C_A,d_AO),[1,2,4,5,3],dims([P,AI,CI,CO,AO]),0,1)];
		    cone_constraints = [cone_constraints, PartialTrace(W_AD_B,4,dims([P,AI,AO,BI,DI,DO])) + PartialTrace(W_AD_C,4,dims([P,AI,AO,CI,DI,DO])) == tensor_id(W_A_D,d_DO) + PermuteSystems(tensor_id(W_D_A,d_AO),[1,2,4,5,3],dims([P,AI,DI,DO,AO]),0,1)];
		    cone_constraints = [cone_constraints, PartialTrace(W_BC_A,2,dims([P,AI,BI,BO,CI,CO])) + PartialTrace(W_BC_D,6,dims([P,BI,BO,CI,CO,DI])) == tensor_id(W_B_C,d_CO) + PermuteSystems(tensor_id(W_C_B,d_BO),[1,2,4,5,3],dims([P,BI,CI,CO,BO]),0,1)];
		    cone_constraints = [cone_constraints, PartialTrace(W_BD_A,2,dims([P,AI,BI,BO,DI,DO])) + PartialTrace(W_BD_C,4,dims([P,BI,BO,CI,DI,DO])) == tensor_id(W_B_D,d_DO) + PermuteSystems(tensor_id(W_D_B,d_BO),[1,2,4,5,3],dims([P,BI,DI,DO,BO]),0,1)];
		    cone_constraints = [cone_constraints, PartialTrace(W_CD_A,2,dims([P,AI,CI,CO,DI,DO])) + PartialTrace(W_CD_B,2,dims([P,BI,CI,CO,DI,DO])) == tensor_id(W_C_D,d_DO) + PermuteSystems(tensor_id(W_D_C,d_CO),[1,2,4,5,3],dims([P,CI,DI,DO,CO]),0,1)];
		    
		    cone_constraints = [cone_constraints, PartialTrace(W_A_B,4,dims([P,AI,AO,BI])) + PartialTrace(W_A_C,4,dims([P,AI,AO,CI])) + PartialTrace(W_A_D,4,dims([P,AI,AO,DI])) == tensor_id(W_0_A,d_AO)];
		    cone_constraints = [cone_constraints, PartialTrace(W_B_A,2,dims([P,AI,BI,BO])) + PartialTrace(W_B_C,4,dims([P,BI,BO,CI])) + PartialTrace(W_B_D,4,dims([P,BI,BO,DI])) == tensor_id(W_0_B,d_BO)];
		    cone_constraints = [cone_constraints, PartialTrace(W_C_A,2,dims([P,AI,CI,CO])) + PartialTrace(W_C_B,2,dims([P,BI,CI,CO])) + PartialTrace(W_C_D,4,dims([P,CI,CO,DI])) == tensor_id(W_0_C,d_CO)];
		    cone_constraints = [cone_constraints, PartialTrace(W_D_A,2,dims([P,AI,DI,DO])) + PartialTrace(W_D_B,2,dims([P,BI,DI,DO])) + PartialTrace(W_D_C,2,dims([P,CI,DI,DO])) == tensor_id(W_0_D,d_DO)];
            
            if d_P ~= 1 %Otherwise this is trvial
                % We only want it to be proportional to the identity, since we are only checking the cone
                % We need to also be careful because sometimes this constraint is trivially satisfied
                % but numerical error makes Yalmip think it's false
                diff = PartialTrace(W_0_A,2,dims([P,AI])) + PartialTrace(W_0_B,2,dims([P,BI])) + PartialTrace(W_0_C,2,dims([P,CI])) + PartialTrace(W_0_D,2,dims([P,DI])) - trace(W_0_A+W_0_B+W_0_C+W_0_D)/d_P*eye(d_P);
                if isa(diff,'sdpvar')
                    cone_constraints = [cone_constraints, diff == 0];
                else
                    cone_constraints = [cone_constraints, matrix_is_equal(diff,zeros(d_P),tol)];
                end
            end        
        otherwise
            error('Currently only implemented up to N=4');
    end
end

