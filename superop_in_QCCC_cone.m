function cone_constraints = superop_in_QCCC_cone(Wr, dims, parties)
%superop_in_QCCC_cone Yalmip constraints for a superinstrument to be in the cone of QC-CCs
%   cone_constraints = superop_in_QCCC_cone(Wr, dims, parties)
%   Wr can be either a superinstrument or a superoperator/process matrix
%   Returns the yalmip SDP constraints for Wr to be in the cone of valid Quantum Circuits with Classical Control of Causal order
%   Note: this doesn't enforce the normalisation of the superoperator
%   
%   Formulation based on: J. Wechs, H. Dourdent, A. A. Abbott, C. Branciard, PRX Quantum 2, 030335 (2021)
%   (See Propositions 5 and 13 therein.)
%
% Requires QETLAB for PermuteSystems, PartialTrace

% Written by Alastair Abbott 2021, last modified 16 August 2022

    % default tolerance
    if ~exist('tol','var')
        tol = 1e-12;
    end

    % First put Wr in canonical ordering (this checks the input validity too)
    % The spaces P,AI,AO,...,F then correspond to dims 1,2,3,...,2*N+2
    if exist('parties','var')
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
       
    constr = [];
    
    switch N
        case 1
            P = 1; d_P = dims(P);
            AI = 2;
            AO = 3; d_AO = dims(AO);
            F = 4;
            
            W_AF = Wr{1};
            for i = 2:R
                W_AF = W_AF + Wr{i};
            end
            W_A = 1/d_AO*PartialTrace(W_AF,[AO,F],dims);

            constr = [constr, superop_in_PSD_cone(Wr)];
            constr = [constr, PartialTrace(W_AF,F,dims) == tensor_id(W_A,d_AO)];
            if d_P ~= 1 % Otherwise this is trvial
                % We only want it to be proportional to the identity, since we are only checking the cone
                % We need to also be careful because sometimes this constraint is trivially satisfied
                % but numerical error makes Yalmip think it's false
                diff = PartialTrace(W_A,2,dims([P,AI])) - trace(W_A)/d_P*eye(d_P);
                if isa(diff,'sdpvar')
                    constr = [constr, diff == 0];
                else
                    constr = [constr, matrix_is_equal(diff,zeros(d_P),tol)];
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
      
            Wr_ABF = cell(1,R);
            Wr_BAF = cell(1,R);
            for i = 1:R
                Wr_ABF{i} = sdpvar(d,d,'hermitian','complex');
                      
                % Define last element from others rather than enforcing an equality constraints
                Wr_BAF{i} = Wr{i} - Wr_ABF{i}; 
            end
            constr = [constr, superop_in_PSD_cone(Wr_ABF), superop_in_PSD_cone(Wr_BAF)];

            W_ABF = Wr_ABF{1};
            W_BAF = Wr_BAF{1};
            for i = 2:R
                W_ABF = W_ABF + Wr_ABF{i};
                W_BAF = W_BAF + Wr_BAF{i};
            end
            
            W_AB = 1/d_BO*PartialTrace(W_ABF,[BO,F],dims);
            W_BA = 1/d_AO*PartialTrace(W_BAF,[AO,F],dims);
            W_A = 1/d_AO*PartialTrace(W_AB,[AO,BI],dims(P:BI));
            W_B = 1/d_BO*PartialTrace(W_BA,[2,4],dims([P,AI,BI,BO]));
            
            constr = [constr, PartialTrace(W_ABF,F,dims) == tensor_id(W_AB,d_BO)];
            constr = [constr, PartialTrace(W_BAF,F,dims) == PermuteSystems(tensor_id(W_BA,d_AO),[P,AI,BI,BO,AO],dims([P,AI,BI,BO,AO]),0,1)];
            constr = [constr, PartialTrace(W_AB,BI,dims(P:BI)) == tensor_id(W_A,d_AO)];
            constr = [constr, PartialTrace(W_BA,AI,dims([P,AI,BI,BO])) == tensor_id(W_B,d_BO)];
            if d_P ~= 1 % Otherwise this is trvial
                % We only want it to be proportional to the identity, since we are only checking the cone
                % We need to also be careful because sometimes this constraint is trivially satisfied
                % but numerical error makes Yalmip think it's false
                diff = PartialTrace(W_A,2,dims([P,AI])) + PartialTrace(W_B,2,dims([P,BI])) - trace(W_A+W_B)/d_P*eye(d_P);
                if isa(diff,'sdpvar')
                    constr = [constr, diff == 0];
                else
                    constr = [constr, matrix_is_equal(diff,zeros(d_P),tol)];
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
            
            Wr_ABCF = cell(1,R);
            Wr_ACBF = cell(1,R);
            Wr_BACF = cell(1,R);
            Wr_BCAF = cell(1,R);
            Wr_CABF = cell(1,R);
            Wr_CBAF = cell(1,R);
            for i = 1:R
                Wr_ABCF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_ACBF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_BACF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_BCAF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_CABF{i} = sdpvar(d,d,'hermitian','complex');
                
                % Define last element from others rather than enforcing an equality constraints
                Wr_CBAF{i} = Wr{i} - (Wr_ABCF{i} + Wr_ACBF{i} + Wr_BACF{i} + Wr_BCAF{i} + Wr_CABF{i}); 
            end
            constr = [constr, superop_in_PSD_cone(Wr_ABCF), superop_in_PSD_cone(Wr_ACBF), superop_in_PSD_cone(Wr_BACF), ...
                              superop_in_PSD_cone(Wr_BCAF), superop_in_PSD_cone(Wr_CABF), superop_in_PSD_cone(Wr_CBAF)];
            
            W_ABCF = Wr_ABCF{1};
            W_ACBF = Wr_ACBF{1};
            W_BACF = Wr_BACF{1};
            W_BCAF = Wr_BCAF{1};
            W_CABF = Wr_CABF{1};
            W_CBAF = Wr_CBAF{1};
            for i = 2:R
                W_ABCF = W_ABCF + Wr_ABCF{i};
                W_ACBF = W_ACBF + Wr_ACBF{i};
                W_BACF = W_BACF + Wr_BACF{i};
                W_BCAF = W_BCAF + Wr_BCAF{i};
                W_CABF = W_CABF + Wr_CABF{i};
                W_CBAF = W_CBAF + Wr_CBAF{i};
            end
            
            W_ABC = 1/d_CO*PartialTrace(W_ABCF,[CO,F],dims);
            W_ACB = 1/d_BO*PartialTrace(W_ACBF,[BO,F],dims);
            W_BAC = 1/d_CO*PartialTrace(W_BACF,[CO,F],dims);
            W_BCA = 1/d_AO*PartialTrace(W_BCAF,[AO,F],dims);
            W_CAB = 1/d_BO*PartialTrace(W_CABF,[BO,F],dims);
            W_CBA = 1/d_AO*PartialTrace(W_CBAF,[AO,F],dims);
            
            W_AB = 1/d_BO*PartialTrace(W_ABC,[5,6],dims([P,AI,AO,BI,BO,CI]));
            W_AC = 1/d_CO*PartialTrace(W_ACB,[4,6],dims([P,AI,AO,BI,CI,CO]));
            W_BA = 1/d_AO*PartialTrace(W_BAC,[3,6],dims([P,AI,AO,BI,BO,CI]));
            W_BC = 1/d_CO*PartialTrace(W_BCA,[2,6],dims([P,AI,BI,BO,CI,CO]));
            W_CA = 1/d_AO*PartialTrace(W_CAB,[3,4],dims([P,AI,AO,BI,CI,CO]));
            W_CB = 1/d_BO*PartialTrace(W_CBA,[2,4],dims([P,AI,BI,BO,CI,CO]));
            
            W_A = 1/d_AO*( PartialTrace(W_AB,[3,4],dims([P,AI,AO,BI])) + PartialTrace(W_AC,[3,4],dims([P,AI,AO,CI])) );
            W_B = 1/d_BO*( PartialTrace(W_BA,[2,4],dims([P,AI,BI,BO])) + PartialTrace(W_BC,[3,4],dims([P,BI,BO,CI])) );
            W_C = 1/d_CO*( PartialTrace(W_CA,[2,4],dims([P,AI,CI,CO])) + PartialTrace(W_CB,[2,4],dims([P,BI,CI,CO])) );
            
            constr = [constr, PartialTrace(W_ABCF,F,dims) == tensor_id(W_ABC,d_CO)];
            constr = [constr, PartialTrace(W_ACBF,F,dims) == PermuteSystems(tensor_id(W_ACB,d_BO),[P,AI,AO,BI,CI,CO,BO],dims([P,AI,AO,BI,CI,CO,BO]),0,1)];
            constr = [constr, PartialTrace(W_BACF,F,dims) == tensor_id(W_BAC,d_CO)];
            constr = [constr, PartialTrace(W_BCAF,F,dims) == PermuteSystems(tensor_id(W_BCA,d_AO),[P,AI,BI,BO,CI,CO,AO],dims([P,AI,BI,BO,CI,CO,AO]),0,1)];
            constr = [constr, PartialTrace(W_CABF,F,dims) == PermuteSystems(tensor_id(W_CAB,d_BO),[P,AI,AO,BI,CI,CO,BO],dims([P,AI,AO,BI,CI,CO,BO]),0,1)];
            constr = [constr, PartialTrace(W_CBAF,F,dims) == PermuteSystems(tensor_id(W_CBA,d_AO),[P,AI,BI,BO,CI,CO,AO],dims([P,AI,BI,BO,CI,CO,AO]),0,1)];
            
            constr = [constr, PartialTrace(W_ABC,6,dims([P,AI,AO,BI,BO,CI])) == tensor_id(W_AB,d_BO)];
            constr = [constr, PartialTrace(W_ACB,4,dims([P,AI,AO,BI,CI,CO])) == tensor_id(W_AC,d_CO)];
            constr = [constr, PartialTrace(W_BAC,6,dims([P,AI,AO,BI,BO,CI])) == PermuteSystems(tensor_id(W_BA,d_AO),[1,2,4,5,3],dims([P,AI,BI,BO,AO]),0,1)];
            constr = [constr, PartialTrace(W_BCA,2,dims([P,AI,BI,BO,CI,CO])) == tensor_id(W_BC,d_CO)];
            constr = [constr, PartialTrace(W_CAB,4,dims([P,AI,AO,BI,CI,CO])) == PermuteSystems(tensor_id(W_CA,d_AO),[1,2,4,5,3],dims([P,AI,CI,CO,AO]),0,1)];
            constr = [constr, PartialTrace(W_CBA,2,dims([P,AI,BI,BO,CI,CO])) == PermuteSystems(tensor_id(W_CB,d_BO),[1,2,4,5,3],dims([P,BI,CI,CO,BO]),0,1)];
            
            constr = [constr, PartialTrace(W_AB,4,dims([P,AI,AO,BI])) + PartialTrace(W_AC,4,dims([P,AI,AO,CI])) == tensor_id(W_A,d_AO)];
            constr = [constr, PartialTrace(W_BA,2,dims([P,AI,BI,BO])) + PartialTrace(W_BC,4,dims([P,BI,BO,CI])) == tensor_id(W_B,d_BO)];
            constr = [constr, PartialTrace(W_CA,2,dims([P,AI,CI,CO])) + PartialTrace(W_CB,2,dims([P,BI,CI,CO])) == tensor_id(W_C,d_CO)];
            
            if d_P ~= 1 %Otherwise this is trvial
                % We only want it to be proportional to the identity, since we are only checking the cone
                % We need to also be careful because sometimes this constraint is trivially satisfied
                % but numerical error makes Yalmip think it's false
                diff = PartialTrace(W_A,2,dims([P,AI])) + PartialTrace(W_B,2,dims([P,BI])) + PartialTrace(W_C,2,dims([P,CI])) - trace(W_A+W_B+W_C)/d_P*eye(d_P);
                if isa(diff,'sdpvar')
                    constr = [constr, diff == 0];
                else
                    constr = [constr, matrix_is_equal(diff,zeros(d_P),tol)];
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
            
            Wr_ABCDF = cell(1,R);
            Wr_ABDCF = cell(1,R);
            Wr_ACBDF = cell(1,R);
            Wr_ACDBF = cell(1,R);
            Wr_ADBCF = cell(1,R);
            Wr_ADCBF = cell(1,R);
            Wr_BACDF = cell(1,R);
            Wr_BADCF = cell(1,R);
            Wr_BCADF = cell(1,R);
            Wr_BCDAF = cell(1,R);
            Wr_BDACF = cell(1,R);
            Wr_BDCAF = cell(1,R);
            Wr_CABDF = cell(1,R);
            Wr_CADBF = cell(1,R);
            Wr_CBADF = cell(1,R);
            Wr_CBDAF = cell(1,R);
            Wr_CDABF = cell(1,R);
            Wr_CDBAF = cell(1,R);
            Wr_DABCF = cell(1,R);
            Wr_DACBF = cell(1,R);
            Wr_DBACF = cell(1,R);
            Wr_DBCAF = cell(1,R);
            Wr_DCABF = cell(1,R);
            Wr_DCBAF = cell(1,R);
            for i = 1:R
                Wr_ABCDF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_ABDCF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_ACBDF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_ACDBF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_ADBCF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_ADCBF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_BACDF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_BADCF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_BCADF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_BCDAF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_BDACF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_BDCAF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_CABDF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_CADBF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_CBADF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_CBDAF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_CDABF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_CDBAF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_DABCF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_DACBF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_DBACF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_DBCAF{i} = sdpvar(d,d,'hermitian','complex');
                Wr_DCABF{i} = sdpvar(d,d,'hermitian','complex');
                
                % Define last element from others rather than enforcing an equality constraints
                % Equal noise on each element of superinstrument
                Wr_DCBAF{i} = Wr{i} ...
                               - (Wr_ABCDF{i} + Wr_ABDCF{i} + Wr_ACBDF{i} + Wr_ACDBF{i} + Wr_ADBCF{i} + Wr_ADCBF{i} ...
                                  + Wr_BACDF{i} + Wr_BADCF{i} + Wr_BCADF{i} + Wr_BCDAF{i} + Wr_BDACF{i} + Wr_BDCAF{i} ...
                                  + Wr_CABDF{i} + Wr_CADBF{i} + Wr_CBADF{i} + Wr_CBDAF{i} + Wr_CDABF{i} +  Wr_CDBAF{i} ...
                                  + Wr_DABCF{i} + Wr_DACBF{i} + Wr_DBACF{i} + Wr_DBCAF{i} + Wr_DCABF{i});
            end
            constr = [constr, superop_in_PSD_cone(Wr_ABCDF), superop_in_PSD_cone(Wr_ABDCF), superop_in_PSD_cone(Wr_ACBDF), superop_in_PSD_cone(Wr_ACDBF), ...
                              superop_in_PSD_cone(Wr_ADBCF), superop_in_PSD_cone(Wr_ADCBF), superop_in_PSD_cone(Wr_BACDF), superop_in_PSD_cone(Wr_BADCF), ...
                              superop_in_PSD_cone(Wr_BCADF), superop_in_PSD_cone(Wr_BCDAF), superop_in_PSD_cone(Wr_BDACF), superop_in_PSD_cone(Wr_BDCAF), ...
                              superop_in_PSD_cone(Wr_CABDF), superop_in_PSD_cone(Wr_CADBF), superop_in_PSD_cone(Wr_CBADF), superop_in_PSD_cone(Wr_CBDAF), ...
                              superop_in_PSD_cone(Wr_CDABF), superop_in_PSD_cone(Wr_CDBAF), superop_in_PSD_cone(Wr_DABCF), superop_in_PSD_cone(Wr_DACBF), ...
                              superop_in_PSD_cone(Wr_DBACF), superop_in_PSD_cone(Wr_DBCAF), superop_in_PSD_cone(Wr_DCABF), superop_in_PSD_cone(Wr_DCBAF)];
            
            W_ABCDF = Wr_ABCDF{1};
            W_ABDCF = Wr_ABDCF{1};
            W_ACBDF = Wr_ACBDF{1};
            W_ACDBF = Wr_ACDBF{1};
            W_ADBCF = Wr_ADBCF{1};
            W_ADCBF = Wr_ADCBF{1};
            W_BACDF = Wr_BACDF{1};
            W_BADCF = Wr_BADCF{1};
            W_BCADF = Wr_BCADF{1};
            W_BCDAF = Wr_BCDAF{1};
            W_BDACF = Wr_BDACF{1};
            W_BDCAF = Wr_BDCAF{1};
            W_CABDF = Wr_CABDF{1};
            W_CADBF = Wr_CADBF{1};
            W_CBADF = Wr_CBADF{1};
            W_CBDAF = Wr_CBDAF{1};
            W_CDABF = Wr_CDABF{1};
            W_CDBAF = Wr_CDBAF{1};
            W_DABCF = Wr_DABCF{1};
            W_DACBF = Wr_DACBF{1};
            W_DBACF = Wr_DBACF{1};
            W_DBCAF = Wr_DBCAF{1};
            W_DCABF = Wr_DCABF{1};
            W_DCBAF = Wr_DCBAF{1};
            for i = 2:n
                W_ABCDF = W_ABCDF + Wr_ABCDF{i};
                W_ABDCF = W_ABDCF + Wr_ABDCF{i};
                W_ACBDF = W_ACBDF + Wr_ACBDF{i};
                W_ACDBF = W_ACDBF + Wr_ACDBF{i};
                W_ADBCF = W_ADBCF + Wr_ADBCF{i};
                W_ADCBF = W_ADCBF + Wr_ADCBF{i};
                W_BACDF = W_BACDF + Wr_BACDF{i};
                W_BADCF = W_BADCF + Wr_BADCF{i};
                W_BCADF = W_BCADF + Wr_BCADF{i};
                W_BCDAF = W_BCDAF + Wr_BCDAF{i};
                W_BDACF = W_BDACF + Wr_BDACF{i};
                W_BDCAF = W_BDCAF + Wr_BDCAF{i};
                W_CABDF = W_CABDF + Wr_CABDF{i};
                W_CADBF = W_CADBF + Wr_CADBF{i};
                W_CBADF = W_CBADF + Wr_CBADF{i};
                W_CBDAF = W_CBDAF + Wr_CBDAF{i};
                W_CDABF = W_CDABF + Wr_CDABF{i};
                W_CDBAF = W_CDBAF + Wr_CDBAF{i};
                W_DABCF = W_DABCF + Wr_DABCF{i};
                W_DACBF = W_DACBF + Wr_DACBF{i};
                W_DBACF = W_DBACF + Wr_DBACF{i};
                W_DBCAF = W_DBCAF + Wr_DBCAF{i};
                W_DCABF = W_DCABF + Wr_DCABF{i};
                W_DCBAF = W_DCBAF + Wr_DCBAF{i};
            end
            
            W_ABCD = 1/d_DO*PartialTrace(W_ABCDF,[DO,F],dims);
            W_ABDC = 1/d_CO*PartialTrace(W_ABDCF,[CO,F],dims);
            W_ACBD = 1/d_DO*PartialTrace(W_ACBDF,[DO,F],dims);
            W_ACDB = 1/d_BO*PartialTrace(W_ACDBF,[BO,F],dims);
            W_ADBC = 1/d_CO*PartialTrace(W_ADBCF,[CO,F],dims);
            W_ADCB = 1/d_BO*PartialTrace(W_ADCBF,[BO,F],dims);
            W_BACD = 1/d_DO*PartialTrace(W_BACDF,[DO,F],dims);
            W_BADC = 1/d_CO*PartialTrace(W_BADCF,[CO,F],dims);
            W_BCAD = 1/d_DO*PartialTrace(W_BCADF,[DO,F],dims);
            W_BCDA = 1/d_AO*PartialTrace(W_BCDAF,[AO,F],dims);
            W_BDAC = 1/d_CO*PartialTrace(W_BDACF,[CO,F],dims);
            W_BDCA = 1/d_AO*PartialTrace(W_BDCAF,[AO,F],dims);
            W_CABD = 1/d_DO*PartialTrace(W_CABDF,[DO,F],dims);
            W_CADB = 1/d_BO*PartialTrace(W_CADBF,[BO,F],dims);
            W_CBAD = 1/d_DO*PartialTrace(W_CBADF,[DO,F],dims);
            W_CBDA = 1/d_AO*PartialTrace(W_CBDAF,[AO,F],dims);
            W_CDAB = 1/d_BO*PartialTrace(W_CDABF,[BO,F],dims);
            W_CDBA = 1/d_AO*PartialTrace(W_CDBAF,[AO,F],dims);
            W_DABC = 1/d_CO*PartialTrace(W_DABCF,[CO,F],dims);
            W_DACB = 1/d_BO*PartialTrace(W_DACBF,[BO,F],dims);
            W_DBAC = 1/d_CO*PartialTrace(W_DBACF,[CO,F],dims);
            W_DBCA = 1/d_AO*PartialTrace(W_DBCAF,[AO,F],dims);
            W_DCAB = 1/d_BO*PartialTrace(W_DCABF,[BO,F],dims);
            W_DCBA = 1/d_AO*PartialTrace(W_DCBAF,[AO,F],dims);
            
            W_ABC = 1/d_CO*PartialTrace(W_ABCD,[7,8],dims([P,AI,AO,BI,BO,CI,CO,DI]));
            W_ABD = 1/d_DO*PartialTrace(W_ABDC,[6,8],dims([P,AI,AO,BI,BO,CI,DI,DO]));
            W_ACB = 1/d_BO*PartialTrace(W_ACBD,[5,8],dims([P,AI,AO,BI,BO,CI,CO,DI]));
            W_ACD = 1/d_DO*PartialTrace(W_ACDB,[4,8],dims([P,AI,AO,BI,CI,CO,DI,DO]));
            W_ADB = 1/d_BO*PartialTrace(W_ADBC,[5,6],dims([P,AI,AO,BI,BO,CI,DI,DO]));
            W_ADC = 1/d_CO*PartialTrace(W_ADCB,[4,6],dims([P,AI,AO,BI,CI,CO,DI,DO]));
            W_BAC = 1/d_CO*PartialTrace(W_BACD,[7,8],dims([P,AI,AO,BI,BO,CI,CO,DI]));
            W_BAD = 1/d_DO*PartialTrace(W_BADC,[6,8],dims([P,AI,AO,BI,BO,CI,DI,DO]));
            W_BCA = 1/d_AO*PartialTrace(W_BCAD,[3,8],dims([P,AI,AO,BI,BO,CI,CO,DI]));
            W_BCD = 1/d_DO*PartialTrace(W_BCDA,[2,8],dims([P,AI,BI,BO,CI,CO,DI,DO]));
            W_BDA = 1/d_AO*PartialTrace(W_BDAC,[3,6],dims([P,AI,AO,BI,BO,CI,DI,DO]));
            W_BDC = 1/d_CO*PartialTrace(W_BDCA,[2,6],dims([P,AI,BI,BO,CI,CO,DI,DO]));
            W_CAB = 1/d_BO*PartialTrace(W_CABD,[5,8],dims([P,AI,AO,BI,BO,CI,CO,DI]));
            W_CAD = 1/d_DO*PartialTrace(W_CADB,[4,8],dims([P,AI,AO,BI,CI,CO,DI,DO]));
            W_CBA = 1/d_AO*PartialTrace(W_CBAD,[3,8],dims([P,AI,AO,BI,BO,CI,CO,DI]));
            W_CBD = 1/d_DO*PartialTrace(W_CBDA,[2,8],dims([P,AI,BI,BO,CI,CO,DI,DO]));
            W_CDA = 1/d_AO*PartialTrace(W_CDAB,[3,4],dims([P,AI,AO,BI,CI,CO,DI,DO]));
            W_CDB = 1/d_BO*PartialTrace(W_CDBA,[2,4],dims([P,AI,BI,BO,CI,CO,DI,DO]));
            W_DAB = 1/d_BO*PartialTrace(W_DABC,[5,6],dims([P,AI,AO,BI,BO,CI,DI,DO]));
            W_DAC = 1/d_CO*PartialTrace(W_DACB,[4,6],dims([P,AI,AO,BI,CI,CO,DI,DO]));
            W_DBA = 1/d_AO*PartialTrace(W_DBAC,[3,6],dims([P,AI,AO,BI,BO,CI,DI,DO]));
            W_DBC = 1/d_CO*PartialTrace(W_DBCA,[2,6],dims([P,AI,BI,BO,CI,CO,DI,DO]));
            W_DCA = 1/d_AO*PartialTrace(W_DCAB,[3,4],dims([P,AI,AO,BI,CI,CO,DI,DO]));
            W_DCB = 1/d_BO*PartialTrace(W_DCBA,[2,4],dims([P,AI,BI,BO,CI,CO,DI,DO]));
            
            W_AB = 1/d_BO*( PartialTrace(W_ABC,[5,6],dims([P,AI,AO,BI,BO,CI])) + PartialTrace(W_ABD,[5,6],dims([P,AI,AO,BI,BO,DI])) );
            W_AC = 1/d_CO*( PartialTrace(W_ACB,[4,6],dims([P,AI,AO,BI,CI,CO])) + PartialTrace(W_ACD,[5,6],dims([P,AI,AO,CI,CO,DI])) );
            W_AD = 1/d_DO*( PartialTrace(W_ADB,[4,6],dims([P,AI,AO,BI,DI,DO])) + PartialTrace(W_ADC,[4,6],dims([P,AI,AO,CI,DI,DO])) );
            W_BA = 1/d_AO*( PartialTrace(W_BAC,[3,6],dims([P,AI,AO,BI,BO,CI])) + PartialTrace(W_BAD,[3,6],dims([P,AI,AO,BI,BO,DI])) );
            W_BC = 1/d_CO*( PartialTrace(W_BCA,[2,6],dims([P,AI,BI,BO,CI,CO])) + PartialTrace(W_BCD,[5,6],dims([P,BI,BO,CI,CO,DI])) );
            W_BD = 1/d_DO*( PartialTrace(W_BDA,[2,6],dims([P,AI,BI,BO,DI,DO])) + PartialTrace(W_BDC,[4,6],dims([P,BI,BO,CI,DI,DO])) );
            W_CA = 1/d_AO*( PartialTrace(W_CAB,[3,4],dims([P,AI,AO,BI,CI,CO])) + PartialTrace(W_CAD,[3,6],dims([P,AI,AO,CI,CO,DI])) );
            W_CB = 1/d_BO*( PartialTrace(W_CBA,[2,4],dims([P,AI,BI,BO,CI,CO])) + PartialTrace(W_CBD,[3,6],dims([P,BI,BO,CI,CO,DI])) );
            W_CD = 1/d_DO*( PartialTrace(W_CDA,[2,6],dims([P,AI,CI,CO,DI,DO])) + PartialTrace(W_CDB,[2,6],dims([P,BI,CI,CO,DI,DO])) );
            W_DA = 1/d_AO*( PartialTrace(W_DAB,[3,4],dims([P,AI,AO,BI,DI,DO])) + PartialTrace(W_DAC,[3,4],dims([P,AI,AO,CI,DI,DO])) );
            W_DB = 1/d_BO*( PartialTrace(W_DBA,[2,4],dims([P,AI,BI,BO,DI,DO])) + PartialTrace(W_DBC,[3,4],dims([P,BI,BO,CI,DI,DO])) );
            W_DC = 1/d_CO*( PartialTrace(W_DCA,[2,4],dims([P,AI,CI,CO,DI,DO])) + PartialTrace(W_DCB,[2,4],dims([P,BI,CI,CO,DI,DO])) );
            
            W_A = 1/d_AO*( PartialTrace(W_AB,[3,4],dims([P,AI,AO,BI])) + PartialTrace(W_AC,[3,4],dims([P,AI,AO,CI])) + PartialTrace(W_AD,[3,4],dims([P,AI,AO,DI])) );
            W_B = 1/d_BO*( PartialTrace(W_BA,[2,4],dims([P,AI,BI,BO])) + PartialTrace(W_BC,[3,4],dims([P,BI,BO,CI])) + PartialTrace(W_BD,[3,4],dims([P,BI,BO,DI])) );
            W_C = 1/d_CO*( PartialTrace(W_CA,[2,4],dims([P,AI,CI,CO])) + PartialTrace(W_CB,[2,4],dims([P,BI,CI,CO])) + PartialTrace(W_CD,[3,4],dims([P,CI,CO,DI])) );
            W_D = 1/d_DO*( PartialTrace(W_DA,[2,4],dims([P,AI,DI,DO])) + PartialTrace(W_DB,[2,4],dims([P,BI,DI,DO])) + PartialTrace(W_DC,[2,4],dims([P,CI,DI,DO])) );
            
            constr = [constr, PartialTrace(W_ABCDF,F,dims) == tensor_id(W_ABCD,d_DO)];
            constr = [constr, PartialTrace(W_ABDCF,F,dims) == PermuteSystems(tensor_id(W_ABDC,d_CO),[P,AI,AO,BI,BO,CI,DI,DO,CO],dims([P,AI,AO,BI,BO,CI,DI,DO,CO]),0,1)];
            constr = [constr, PartialTrace(W_ACBDF,F,dims) == tensor_id(W_ACBD,d_DO)];
            constr = [constr, PartialTrace(W_ACDBF,F,dims) == PermuteSystems(tensor_id(W_ACDB,d_BO),[P,AI,AO,BI,CI,CO,DI,DO,BO],dims([P,AI,AO,BI,CI,CO,DI,DO,BO]),0,1)];
            constr = [constr, PartialTrace(W_ADBCF,F,dims) == PermuteSystems(tensor_id(W_ADBC,d_CO),[P,AI,AO,BI,BO,CI,DI,DO,CO],dims([P,AI,AO,BI,BO,CI,DI,DO,CO]),0,1)];
            constr = [constr, PartialTrace(W_ADCBF,F,dims) == PermuteSystems(tensor_id(W_ADCB,d_BO),[P,AI,AO,BI,CI,CO,DI,DO,BO],dims([P,AI,AO,BI,CI,CO,DI,DO,BO]),0,1)];
            constr = [constr, PartialTrace(W_BACDF,F,dims) == tensor_id(W_BACD,d_DO)];
            constr = [constr, PartialTrace(W_BADCF,F,dims) == PermuteSystems(tensor_id(W_BADC,d_CO),[P,AI,AO,BI,BO,CI,DI,DO,CO],dims([P,AI,AO,BI,BO,CI,DI,DO,CO]),0,1)];
            constr = [constr, PartialTrace(W_BCADF,F,dims) == tensor_id(W_BCAD,d_DO)];
            constr = [constr, PartialTrace(W_BCDAF,F,dims) == PermuteSystems(tensor_id(W_BCDA,d_AO),[P,AI,BI,BO,CI,CO,DI,DO,AO],dims([P,AI,BI,BO,CI,CO,DI,DO,AO]),0,1)];
            constr = [constr, PartialTrace(W_BDACF,F,dims) == PermuteSystems(tensor_id(W_BDAC,d_CO),[P,AI,AO,BI,BO,CI,DI,DO,CO],dims([P,AI,AO,BI,BO,CI,DI,DO,CO]),0,1)];
            constr = [constr, PartialTrace(W_BDCAF,F,dims) == PermuteSystems(tensor_id(W_BDCA,d_AO),[P,AI,BI,BO,CI,CO,DI,DO,AO],dims([P,AI,BI,BO,CI,CO,DI,DO,AO]),0,1)];
            constr = [constr, PartialTrace(W_CABDF,F,dims) == tensor_id(W_CABD,d_DO)];
            constr = [constr, PartialTrace(W_CADBF,F,dims) == PermuteSystems(tensor_id(W_CADB,d_BO),[P,AI,AO,BI,CI,CO,DI,DO,BO],dims([P,AI,AO,BI,CI,CO,DI,DO,BO]),0,1)];
            constr = [constr, PartialTrace(W_CBADF,F,dims) == tensor_id(W_CBAD,d_DO)];
            constr = [constr, PartialTrace(W_CBDAF,F,dims) == PermuteSystems(tensor_id(W_CBDA,d_AO),[P,AI,BI,BO,CI,CO,DI,DO,AO],dims([P,AI,BI,BO,CI,CO,DI,DO,AO]),0,1)];
            constr = [constr, PartialTrace(W_CDABF,F,dims) == PermuteSystems(tensor_id(W_CDAB,d_BO),[P,AI,AO,BI,CI,CO,DI,DO,BO],dims([P,AI,AO,BI,CI,CO,DI,DO,BO]),0,1)];
            constr = [constr, PartialTrace(W_CDBAF,F,dims) == PermuteSystems(tensor_id(W_CDBA,d_AO),[P,AI,BI,BO,CI,CO,DI,DO,AO],dims([P,AI,BI,BO,CI,CO,DI,DO,AO]),0,1)];
            constr = [constr, PartialTrace(W_DABCF,F,dims) == PermuteSystems(tensor_id(W_DABC,d_CO),[P,AI,AO,BI,BO,CI,DI,DO,CO],dims([P,AI,AO,BI,BO,CI,DI,DO,CO]),0,1)];
            constr = [constr, PartialTrace(W_DACBF,F,dims) == PermuteSystems(tensor_id(W_DACB,d_BO),[P,AI,AO,BI,CI,CO,DI,DO,BO],dims([P,AI,AO,BI,CI,CO,DI,DO,BO]),0,1)];
            constr = [constr, PartialTrace(W_DBACF,F,dims) == PermuteSystems(tensor_id(W_DBAC,d_CO),[P,AI,AO,BI,BO,CI,DI,DO,CO],dims([P,AI,AO,BI,BO,CI,DI,DO,CO]),0,1)];
            constr = [constr, PartialTrace(W_DBCAF,F,dims) == PermuteSystems(tensor_id(W_DBCA,d_AO),[P,AI,BI,BO,CI,CO,DI,DO,AO],dims([P,AI,BI,BO,CI,CO,DI,DO,AO]),0,1)];
            constr = [constr, PartialTrace(W_DCABF,F,dims) == PermuteSystems(tensor_id(W_DCAB,d_BO),[P,AI,AO,BI,CI,CO,DI,DO,BO],dims([P,AI,AO,BI,CI,CO,DI,DO,BO]),0,1)];
            constr = [constr, PartialTrace(W_DCBAF,F,dims) == PermuteSystems(tensor_id(W_DCBA,d_AO),[P,AI,BI,BO,CI,CO,DI,DO,AO],dims([P,AI,BI,BO,CI,CO,DI,DO,AO]),0,1)];
            
            constr = [constr, PartialTrace(W_ABCD,8,dims([P,AI,AO,BI,BO,CI,CO,DI])) == tensor_id(W_ABC,d_CO)];
            constr = [constr, PartialTrace(W_ABDC,6,dims([P,AI,AO,BI,BO,CI,DI,DO])) == tensor_id(W_ABD,d_DO)];
            constr = [constr, PartialTrace(W_ACBD,8,dims([P,AI,AO,BI,BO,CI,CO,DI])) == PermuteSystems(tensor_id(W_ACB,d_BO),[1,2,3,4,6,7,5],dims([P,AI,AO,BI,CI,CO,BO]),0,1)];
            constr = [constr, PartialTrace(W_ACDB,4,dims([P,AI,AO,BI,CI,CO,DI,DO])) == tensor_id(W_ACD,d_DO)];
            constr = [constr, PartialTrace(W_ADBC,6,dims([P,AI,AO,BI,BO,CI,DI,DO])) == PermuteSystems(tensor_id(W_ADB,d_BO),[1,2,3,4,6,7,5],dims([P,AI,AO,BI,DI,DO,BO]),0,1)];
            constr = [constr, PartialTrace(W_ADCB,4,dims([P,AI,AO,BI,CI,CO,DI,DO])) == PermuteSystems(tensor_id(W_ADC,d_CO),[1,2,3,4,6,7,5],dims([P,AI,AO,CI,DI,DO,CO]),0,1)];
            constr = [constr, PartialTrace(W_BACD,8,dims([P,AI,AO,BI,BO,CI,CO,DI])) == tensor_id(W_BAC,d_CO)];
            constr = [constr, PartialTrace(W_BADC,6,dims([P,AI,AO,BI,BO,CI,DI,DO])) == tensor_id(W_BAD,d_DO)];
            constr = [constr, PartialTrace(W_BCAD,8,dims([P,AI,AO,BI,BO,CI,CO,DI])) == PermuteSystems(tensor_id(W_BCA,d_AO),[1,2,4,5,6,7,3],dims([P,AI,BI,BO,CI,CO,AO]),0,1)];
            constr = [constr, PartialTrace(W_BCDA,2,dims([P,AI,BI,BO,CI,CO,DI,DO])) == tensor_id(W_BCD,d_DO)];
            constr = [constr, PartialTrace(W_BDAC,6,dims([P,AI,AO,BI,BO,CI,DI,DO])) == PermuteSystems(tensor_id(W_BDA,d_AO),[1,2,4,5,6,7,3],dims([P,AI,BI,BO,DI,DO,AO]),0,1)];
            constr = [constr, PartialTrace(W_BDCA,2,dims([P,AI,BI,BO,CI,CO,DI,DO])) == PermuteSystems(tensor_id(W_BDC,d_CO),[1,2,3,4,6,7,5],dims([P,BI,BO,CI,DI,DO,CO]),0,1)];
            constr = [constr, PartialTrace(W_CABD,8,dims([P,AI,AO,BI,BO,CI,CO,DI])) == PermuteSystems(tensor_id(W_CAB,d_BO),[1,2,3,4,6,7,5],dims([P,AI,AO,BI,CI,CO,BO]),0,1)];
            constr = [constr, PartialTrace(W_CADB,4,dims([P,AI,AO,BI,CI,CO,DI,DO])) == tensor_id(W_CAD,d_DO)];
            constr = [constr, PartialTrace(W_CBAD,8,dims([P,AI,AO,BI,BO,CI,CO,DI])) == PermuteSystems(tensor_id(W_CBA,d_AO),[1,2,4,5,6,7,3],dims([P,AI,BI,BO,CI,CO,AO]),0,1)];
            constr = [constr, PartialTrace(W_CBDA,2,dims([P,AI,BI,BO,CI,CO,DI,DO])) == tensor_id(W_CBD,d_DO)];
            constr = [constr, PartialTrace(W_CDAB,4,dims([P,AI,AO,BI,CI,CO,DI,DO])) == PermuteSystems(tensor_id(W_CDA,d_AO),[1,2,4,5,6,7,3],dims([P,AI,CI,CO,DI,DO,AO]),0,1)];
            constr = [constr, PartialTrace(W_CDBA,2,dims([P,AI,BI,BO,CI,CO,DI,DO])) == PermuteSystems(tensor_id(W_CDB,d_BO),[1,2,4,5,6,7,3],dims([P,BI,CI,CO,DI,DO,BO]),0,1)];
            constr = [constr, PartialTrace(W_DABC,6,dims([P,AI,AO,BI,BO,CI,DI,DO])) == PermuteSystems(tensor_id(W_DAB,d_BO),[1,2,3,4,6,7,5],dims([P,AI,AO,BI,DI,DO,BO]),0,1)];
            constr = [constr, PartialTrace(W_DACB,4,dims([P,AI,AO,BI,CI,CO,DI,DO])) == PermuteSystems(tensor_id(W_DAC,d_CO),[1,2,3,4,6,7,5],dims([P,AI,AO,CI,DI,DO,CO]),0,1)];
            constr = [constr, PartialTrace(W_DBAC,6,dims([P,AI,AO,BI,BO,CI,DI,DO])) == PermuteSystems(tensor_id(W_DBA,d_AO),[1,2,4,5,6,7,3],dims([P,AI,BI,BO,DI,DO,AO]),0,1)];
            constr = [constr, PartialTrace(W_DBCA,2,dims([P,AI,BI,BO,CI,CO,DI,DO])) == PermuteSystems(tensor_id(W_DBC,d_CO),[1,2,3,4,6,7,5],dims([P,BI,BO,CI,DI,DO,CO]),0,1)];
            constr = [constr, PartialTrace(W_DCAB,4,dims([P,AI,AO,BI,CI,CO,DI,DO])) == PermuteSystems(tensor_id(W_DCA,d_AO),[1,2,4,5,6,7,3],dims([P,AI,CI,CO,DI,DO,AO]),0,1)];
            constr = [constr, PartialTrace(W_DCBA,2,dims([P,AI,BI,BO,CI,CO,DI,DO])) == PermuteSystems(tensor_id(W_DCB,d_BO),[1,2,4,5,6,7,3],dims([P,BI,CI,CO,DI,DO,BO]),0,1)];
            
            constr = [constr, PartialTrace(W_ABC,6,dims([P,AI,AO,BI,BO,CI])) + PartialTrace(W_ABD,6,dims([P,AI,AO,BI,BO,DI])) == tensor_id(W_AB,d_BO)];
            constr = [constr, PartialTrace(W_ACB,4,dims([P,AI,AO,BI,CI,CO])) + PartialTrace(W_ACD,6,dims([P,AI,AO,CI,CO,DI])) == tensor_id(W_AC,d_CO)];
            constr = [constr, PartialTrace(W_ADB,4,dims([P,AI,AO,BI,DI,DO])) + PartialTrace(W_ADC,4,dims([P,AI,AO,CI,DI,DO])) == tensor_id(W_AD,d_DO)];
            constr = [constr, PartialTrace(W_BAC,6,dims([P,AI,AO,BI,BO,CI])) + PartialTrace(W_BAD,6,dims([P,AI,AO,BI,BO,DI])) == PermuteSystems(tensor_id(W_BA,d_AO),[1,2,4,5,3],dims([P,AI,BI,BO,AO]),0,1)];
            constr = [constr, PartialTrace(W_BCA,2,dims([P,AI,BI,BO,CI,CO])) + PartialTrace(W_BCD,6,dims([P,BI,BO,CI,CO,DI])) == tensor_id(W_BC,d_CO)];
            constr = [constr, PartialTrace(W_BDA,2,dims([P,AI,BI,BO,DI,DO])) + PartialTrace(W_BDC,4,dims([P,BI,BO,CI,DI,DO])) == tensor_id(W_BD,d_DO)];
            constr = [constr, PartialTrace(W_CAB,4,dims([P,AI,AO,BI,CI,CO])) + PartialTrace(W_CAD,6,dims([P,AI,AO,CI,CO,DI])) == PermuteSystems(tensor_id(W_CA,d_AO),[1,2,4,5,3],dims([P,AI,CI,CO,AO]),0,1)];
            constr = [constr, PartialTrace(W_CBA,2,dims([P,AI,BI,BO,CI,CO])) + PartialTrace(W_CBD,6,dims([P,BI,BO,CI,CO,DI])) == PermuteSystems(tensor_id(W_CB,d_BO),[1,2,4,5,3],d([P,BI,CI,CO,BO]),0,1)];
            constr = [constr, PartialTrace(W_CDA,2,dims([P,AI,CI,CO,DI,DO])) + PartialTrace(W_CDB,2,dims([P,BI,CI,CO,DI,DO])) == tensor_id(W_CD,d_DO)];
            constr = [constr, PartialTrace(W_DAB,4,dims([P,AI,AO,BI,DI,DO])) + PartialTrace(W_DAC,4,dims([P,AI,AO,CI,DI,DO])) == PermuteSystems(tensor_id(W_DA,d_AO),[1,2,4,5,3],dims([P,AI,DI,DO,AO]),0,1)];
            constr = [constr, PartialTrace(W_DBA,2,dims([P,AI,BI,BO,DI,DO])) + PartialTrace(W_DBC,4,dims([P,BI,BO,CI,DI,DO])) == PermuteSystems(tensor_id(W_DB,d_BO),[1,2,4,5,3],dims([P,BI,DI,DO,BO]),0,1)];
            constr = [constr, PartialTrace(W_DCA,2,dims([P,AI,CI,CO,DI,DO])) + PartialTrace(W_DCB,2,dims([P,BI,CI,CO,DI,DO])) == PermuteSystems(tensor_id(W_DC,d_CO),[1,2,4,5,3],dims([P,CI,DI,DO,CO]),0,1)];
            
            constr = [constr, PartialTrace(W_AB,4,dims([P,AI,AO,BI])) + PartialTrace(W_AC,4,dims([P,AI,AO,CI])) + PartialTrace(W_AD,4,dims([P,AI,AO,DI])) == tensor_id(W_A,d_AO)];
            constr = [constr, PartialTrace(W_BA,2,dims([P,AI,BI,BO])) + PartialTrace(W_BC,4,dims([P,BI,BO,CI])) + PartialTrace(W_BD,4,dims([P,BI,BO,DI])) == tensor_id(W_B,d_BO)];
            constr = [constr, PartialTrace(W_CA,2,dims([P,AI,CI,CO])) + PartialTrace(W_CB,2,dims([P,BI,CI,CO])) + PartialTrace(W_CD,4,dims([P,CI,CO,DI])) == tensor_id(W_C,d_CO)];
            constr = [constr, PartialTrace(W_DA,2,dims([P,AI,DI,DO])) + PartialTrace(W_DB,2,dims([P,BI,DI,DO])) + PartialTrace(W_DC,2,dims([P,CI,DI,DO])) == tensor_id(W_D,d_DO)];
            
            if d_P ~= 1 %Otherwise this is trvial
                % We only want it to be proportional to the identity, since we are only checking the cone
                % We need to also be careful because sometimes this constraint is trivially satisfied
                % but numerical error makes Yalmip think it's false
                diff = PartialTrace(W_A,2,dims([P,AI])) + PartialTrace(W_B,2,dims([P,BI])) + PartialTrace(W_C,2,dims([P,CI])) + PartialTrace(W_D,2,dims([P,DI])) - trace(W_A+W_B+W_C+W_D)/d_P*eye(d_P);
                if isa(diff,'sdpvar')
                    constr = [constr, diff == 0];
                else
                    constr = [constr, matrix_is_equal(diff,zeros(d_P),tol)];
                end
            end        
        otherwise
            error('Currently only implemented up to N=4');
    end
    cone_constraints = constr;
end

