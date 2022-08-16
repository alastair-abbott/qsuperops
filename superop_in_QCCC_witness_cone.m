function cone_constraints = superop_in_QCCC_witness_cone(Sr, dims, parties)
%superop_in_QCCC_witness_cone Yalmip constraints for a set of matrices to be in the cone of witnesses for QC-CCs
%   cone_constraints = superop_in_QCCC_witness_cone(Sr, dims, parties) 
%   Sr can be either a single witness S, or a superinstrument witness Sr
%   Returns the yalmip constraints, i.e. for the dual cone
%
% Requires QETLAB for PermuteSystems, PartialTrace

% Written by Alastair Abbott 2021, last modified 16 August 2022

    % First put Sr in canonical ordering (this checks the input validity too)
    % The spaces P,AI,AO,...,F then correspond to dims 1,2,3,...,2*N+2
    [Sr, dims, parties] = superop_to_canonical_ordering(Sr, dims, parties);

    % treat a process matrix witness as a single element superinstrument witness
    if ~iscell(Sr)
        Sr = {Sr};
    end
    
    R = length(Sr);
    N = length(parties) - 2;
    
    d = prod(dims);
    
    %% We treat each N separately
    
    constr = [];
    
    switch N
        case 1
            %%
            P = 1; d_P = dims(P);
            AI = 2; 
            AO = 3; 
            F = 4;
            
            S_PAF = sdpvar(d,d,'hermitian','complex');
            
            if d_P ~= 1
                S_P = sdpvar(d,d,'hermitian','complex');
            else
                S_P = zeros(d,d);
            end
            
            T_PAF = cell(1,R);
            T_P = cell(1,R);
            for i = 1:R
                % %% Full versions of constraints to aid readability, before eliminating equality constraints
                % T_PAF{i} = sdpvar(d,d,'hermitian','complex');
                %
                % T_P{i} = T_PAF{i} + S_PAF;
                %
                % constr = [constr, Sr{i} == T_P{i} + S_P];
                
                T_P{i} = Sr{i} - S_P;
                T_PAF{i} = T_P{i} - S_PAF;
                
                constr = [constr, T_PAF{i} >= 0];
            end
            
            S_PAF_proj = S_PAF - (trRep(S_PAF,F,dims) - trRep(S_PAF,[AO,F],dims));
            
            constr = [constr, S_PAF_proj == 0];
            
            if d_P ~= 1
                constr = [constr, S_P - (trRep(S_P,[AI,AO,F],dims) - trRep(S_P,[P,AI,AO,F],dims)) == 0];
            end
            
        case 2
            %%
            % Recall we now have a canically ordered process and we've grouped the "sub"-spaces
            P = 1; d_P = dims(P);
            AI = 2; 
            AO = 3; 
            BI = 4; 
            BO = 5; 
            F = 6;
            
            S_PABF = sdpvar(d,d,'hermitian','complex');
            S_PBAF = sdpvar(d,d,'hermitian','complex');
            
            if d_P ~= 1
                S_P = sdpvar(d,d,'hermitian','complex');
            else
                S_P = zeros(d,d);
            end
            
            T_PABF = cell(1,R);
            T_PBAF = cell(1,R);
            T_P = cell(1,R);
            for i = 1:R
                % %% Full versions of constraints to aid readability, before eliminating equality constraints
                % T_PABF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PBAF{i} = sdpvar(d,d,'hermitian','complex');
                %
                % T_P{i} = T_PABF{i} + S_PABF;
                % constr = [constr, T_P{i} == T_PBAF{i} + S_PBAF];
                %
                % constr = [constr, Sr{i} == T_P{i} + S_P];
                
                T_P{i} = Sr{i} - S_P;
                T_PABF{i} = T_P{i} - S_PABF;
                T_PBAF{i} = T_P{i} - S_PBAF;
                
                constr = [constr, T_PABF{i} >= 0, T_PBAF{i} >= 0];
            end
            
            S_PABF_proj = S_PABF - (trRep(S_PABF,F,dims) - trRep(S_PABF,[BO,F],dims));
            S_PABF_proj = S_PABF_proj - (trRep(S_PABF_proj,[BI,BO,F],dims) - trRep(S_PABF_proj,[AO,BI,BO,F],dims));
            S_PBAF_proj = S_PBAF - (trRep(S_PBAF,F,dims) - trRep(S_PBAF,[AO,F],dims));
            S_PBAF_proj = S_PBAF_proj - (trRep(S_PBAF_proj,[AI,AO,F],dims) - trRep(S_PBAF_proj,[AI,AO,BO,F],dims));
            
            constr = [constr, S_PABF_proj == 0, S_PBAF_proj == 0];
            
            if d_P ~= 1
                constr = [constr, S_P - (trRep(S_P,[AI,AO,BI,BO,F],dims) - trRep(S_P,[P,AI,AO,BI,BO,F],dims)) == 0];
            end
            
        case 3
            %%
            P = 1; d_P = dims(P);
            AI = 2; 
            AO = 3; 
            BI = 4; 
            BO = 5;
            CI = 6;
            CO = 7; 
            F = 8;
            
            S_PABCF = sdpvar(d,d,'hermitian','complex');
            S_PACBF = sdpvar(d,d,'hermitian','complex');
            S_PBACF = sdpvar(d,d,'hermitian','complex');
            S_PBCAF = sdpvar(d,d,'hermitian','complex');
            S_PCABF = sdpvar(d,d,'hermitian','complex');
            S_PCBAF = sdpvar(d,d,'hermitian','complex');
            
            S_PA = sdpvar(d,d,'hermitian','complex');
            S_PB = sdpvar(d,d,'hermitian','complex');
            S_PC = sdpvar(d,d,'hermitian','complex');
            
            if d_P ~= 1
                S_P = sdpvar(d,d,'hermitian','complex');
            else
                S_P = zeros(d,d);
            end
            
            T_PABCF = cell(1,R);
            T_PACBF = cell(1,R);
            T_PBACF = cell(1,R);
            T_PBCAF = cell(1,R);
            T_PCABF = cell(1,R);
            T_PCBAF = cell(1,R);
            T_PA = cell(1,R);
            T_PB = cell(1,R);
            T_PC = cell(1,R);
            T_P = cell(1,R);
            for i = 1:R
                % %% Full versions of constraints to aid readability, before eliminating equality constraints
                % T_PABCF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PACBF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PBACF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PBCAF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PCABF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PCBAF{i} = sdpvar(d,d,'hermitian','complex');
                %
                % T_PA{i} = T_PABCF{i} + S_PABCF;
                % constr = [constr, T_PA{i} == T_PACBF{i} + S_PACBF];
                % T_PB{i} = T_PBACF{i} + S_PBACF;
                % constr = [constr, T_PB{i} == T_PBCAF{i} + S_PBCAF];
                % T_PC{i} = T_PCABF{i} + S_PCABF;
                % constr = [constr, T_PC{i} == T_PCBAF{i} + S_PCBAF];
                %
                % T_P{i} = T_PA{i} + S_PA;
                % constr = [constr, T_P{i} == T_PB{i} + S_PB, T_P{i} == T_PC{i} + S_PC];
                %
                % constr = [constr, Sr{i} == T_P{i} + S_P];
                
                T_P{i} = Sr{i} - S_P;
                
                T_PA{i} = T_P{i} - S_PA;
                T_PB{i} = T_P{i} - S_PB;
                T_PC{i} = T_P{i} - S_PC;
                
                T_PABCF{i} = T_PA{i} - S_PABCF;
                T_PACBF{i} = T_PA{i} - S_PACBF;

                T_PBACF{i} = T_PB{i} - S_PBACF;
                T_PBCAF{i} = T_PB{i} - S_PBCAF;
                
                T_PCABF{i} = T_PC{i} - S_PCABF;
                T_PCBAF{i} = T_PC{i} - S_PCBAF; 
                
                constr = [constr, T_PABCF{i} >= 0, T_PACBF{i} >= 0, T_PBACF{i} >= 0, T_PBCAF{i} >= 0, T_PCABF{i} >= 0, T_PCBAF{i} >= 0];
            end
            
            S_PABCF_proj = S_PABCF - ( trRep(S_PABCF,F,dims) - trRep(S_PABCF,[CO,F],dims) );
            S_PABCF_proj = S_PABCF_proj - ( trRep(S_PABCF_proj,[CI,CO,F],dims) - trRep(S_PABCF_proj,[BO,CI,CO,F],dims) );
            S_PACBF_proj = S_PACBF - ( trRep(S_PACBF,F,dims) - trRep(S_PACBF,[BO,F],dims) );
            S_PACBF_proj = S_PACBF_proj - ( trRep(S_PACBF_proj,[BI,BO,F],dims) - trRep(S_PACBF_proj,[BI,BO,CO,F],dims) );
            S_PBACF_proj = S_PBACF - ( trRep(S_PBACF,F,dims) - trRep(S_PBACF,[CO,F],dims) );
            S_PBACF_proj = S_PBACF_proj - ( trRep(S_PBACF_proj,[CI,CO,F],dims) - trRep(S_PBACF_proj,[AO,CI,CO,F],dims) );
            S_PBCAF_proj = S_PBCAF - ( trRep(S_PBCAF,F,dims) - trRep(S_PBCAF,[AO,F],dims) );
            S_PBCAF_proj = S_PBCAF_proj - ( trRep(S_PBCAF_proj,[AI,AO,F],dims) - trRep(S_PBCAF_proj,[AI,AO,CO,F],dims) );
            S_PCABF_proj = S_PCABF - ( trRep(S_PCABF,F,dims) - trRep(S_PCABF,[BO,F],dims) );
            S_PCABF_proj = S_PCABF_proj - ( trRep(S_PCABF_proj,[BI,BO,F],dims) - trRep(S_PCABF_proj,[AO,BI,BO,F],dims) );
            S_PCBAF_proj = S_PCBAF - ( trRep(S_PCBAF,F,dims) - trRep(S_PCBAF,[AO,F],dims) );
            S_PCBAF_proj = S_PCBAF_proj - ( trRep(S_PCBAF_proj,[AI,AO,F],dims) - trRep(S_PCBAF_proj,[AI,AO,BO,F],dims) );
              
            constr = [constr, S_PABCF_proj == 0, S_PACBF_proj == 0, S_PBACF_proj == 0, S_PBCAF_proj == 0, S_PCABF_proj == 0, S_PCBAF_proj == 0];
       
            constr = [constr, S_PA - ( trRep(S_PA,[BI,BO,CI,CO,F],dims) - trRep(S_PA,[AO,BI,BO,CI,CO,F],dims) ) == 0];
            constr = [constr, S_PB - ( trRep(S_PB,[AI,AO,CI,CO,F],dims) - trRep(S_PB,[AI,AO,BO,CI,CO,F],dims) ) == 0];
            constr = [constr, S_PC - ( trRep(S_PC,[AI,AO,BI,BO,F],dims) - trRep(S_PC,[AI,AO,BI,BO,CO,F],dims) ) == 0];
            
            if d_P ~= 1
                constr = [constr, S_P - ( trRep(S_P,[AI,AO,BI,BO,CI,CO,F],dims) - trRep(S_P,[P,AI,AO,BI,BO,CI,CO,F],dims) ) == 0];
            end
            
        case 4
            %%
            P = 1; d_P = dims(P);
            AI = 2; 
            AO = 3; 
            BI = 4; 
            BO = 5;
            CI = 6;
            CO = 7;
            DI = 8;
            DO = 9;
            F = 10;
            
            S_PABCDF = sdpvar(d,d,'hermitian','complex');
            S_PABDCF = sdpvar(d,d,'hermitian','complex');
            S_PACBDF = sdpvar(d,d,'hermitian','complex');
            S_PACDBF = sdpvar(d,d,'hermitian','complex');
            S_PADBCF = sdpvar(d,d,'hermitian','complex');
            S_PADCBF = sdpvar(d,d,'hermitian','complex');
            S_PBACDF = sdpvar(d,d,'hermitian','complex');
            S_PBADCF = sdpvar(d,d,'hermitian','complex');
            S_PBCADF = sdpvar(d,d,'hermitian','complex');
            S_PBCDAF = sdpvar(d,d,'hermitian','complex');
            S_PBDACF = sdpvar(d,d,'hermitian','complex');
            S_PBDCAF = sdpvar(d,d,'hermitian','complex');
            S_PCABDF = sdpvar(d,d,'hermitian','complex');
            S_PCADBF = sdpvar(d,d,'hermitian','complex');
            S_PCBADF = sdpvar(d,d,'hermitian','complex');
            S_PCBDAF = sdpvar(d,d,'hermitian','complex');
            S_PCDABF = sdpvar(d,d,'hermitian','complex');
            S_PCDBAF = sdpvar(d,d,'hermitian','complex');
            S_PDABCF = sdpvar(d,d,'hermitian','complex');
            S_PDACBF = sdpvar(d,d,'hermitian','complex');
            S_PDBACF = sdpvar(d,d,'hermitian','complex');
            S_PDBCAF = sdpvar(d,d,'hermitian','complex');
            S_PDCABF = sdpvar(d,d,'hermitian','complex');
            S_PDCBAF = sdpvar(d,d,'hermitian','complex');
            
            S_PAB = sdpvar(d,d,'hermitian','complex');
            S_PAC = sdpvar(d,d,'hermitian','complex');
            S_PAD = sdpvar(d,d,'hermitian','complex');
            S_PBA = sdpvar(d,d,'hermitian','complex');
            S_PBC = sdpvar(d,d,'hermitian','complex');
            S_PBD = sdpvar(d,d,'hermitian','complex');
            S_PCA = sdpvar(d,d,'hermitian','complex');
            S_PCB = sdpvar(d,d,'hermitian','complex');
            S_PCD = sdpvar(d,d,'hermitian','complex');
            S_PDA = sdpvar(d,d,'hermitian','complex');
            S_PDB = sdpvar(d,d,'hermitian','complex');
            S_PDC = sdpvar(d,d,'hermitian','complex');
            
            S_PA = sdpvar(d,d,'hermitian','complex');
            S_PB = sdpvar(d,d,'hermitian','complex');
            S_PC = sdpvar(d,d,'hermitian','complex');
            S_PD = sdpvar(d,d,'hermitian','complex');
            
            if d_P ~= 1
                S_P = sdpvar(d,d,'hermitian','complex');
            else
                S_P = zeros(d,d);
            end
            
            T_PABCDF = cell(1,R);
            T_PABDCF = cell(1,R);
            T_PACBDF = cell(1,R);
            T_PACDBF = cell(1,R);
            T_PADBCF = cell(1,R);
            T_PADCBF = cell(1,R);
            T_PBACDF = cell(1,R);
            T_PBADCF = cell(1,R);
            T_PBCADF = cell(1,R);
            T_PBCDAF = cell(1,R);
            T_PBDACF = cell(1,R);
            T_PBDCAF = cell(1,R);
            T_PCABDF = cell(1,R);
            T_PCADBF = cell(1,R);
            T_PCBADF = cell(1,R);
            T_PCBDAF = cell(1,R);
            T_PCDABF = cell(1,R);
            T_PCDBAF = cell(1,R);
            T_PDABCF = cell(1,R);
            T_PDACBF = cell(1,R);
            T_PDBACF = cell(1,R);
            T_PDBCAF = cell(1,R);
            T_PDCABF = cell(1,R);
            T_PDCBAF = cell(1,R);
            T_PAB = cell(1,R);
            T_PAC = cell(1,R);
            T_PAD = cell(1,R);
            T_PBA = cell(1,R);
            T_PBC = cell(1,R);
            T_PBD = cell(1,R);
            T_PCA = cell(1,R);
            T_PCB = cell(1,R);
            T_PCD = cell(1,R);
            T_PDA = cell(1,R);
            T_PDB = cell(1,R);
            T_PDC = cell(1,R);
            T_PA = cell(1,R);
            T_PB = cell(1,R);
            T_PC = cell(1,R);
            T_PD = cell(1,R);
            T_P = cell(1,R);
            for i = 1:R
                % %% Full versions of constraints to aid readability, before eliminating equality constraints
                % T_PABCDF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PABDCF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PACBDF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PACDBF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PADBCF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PADCBF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PBACDF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PBADCF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PBCADF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PBCDAF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PBDACF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PBDCAF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PCABDF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PCADBF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PCBADF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PCBDAF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PCDABF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PCDBAF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PDABCF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PDACBF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PDBACF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PDBCAF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PDCABF{i} = sdpvar(d,d,'hermitian','complex');
                % T_PDCBAF{i} = sdpvar(d,d,'hermitian','complex');
                %
                % T_PAB{i} = T_PABCDF{i} + S_PABCDF;
                % constr = [constr, T_PAB{i} == T_PABDCF{i} + S_PABDCF];
                % T_PAC{i} = T_PACBDF{i} + S_PACBDF;
                % constr = [constr, T_PAC{i} == T_PACDBF{i} + S_PACDBF];
                % T_PAD{i} = T_PADBCF{i} + S_PADBCF;
                % constr = [constr, T_PAD{i} == T_PADCBF{i} + S_PADCBF];
                % T_PBA{i} = T_PBACDF{i} + S_PBACDF;
                % constr = [constr, T_PBA{i} == T_PBADCF{i} + S_PBADCF];
                % T_PBC{i} = T_PBCADF{i} + S_PBCADF;
                % constr = [constr, T_PBC{i} == T_PBCDAF{i} + S_PBCDAF];
                % T_PBD{i} = T_PBDACF{i} + S_PBDACF;
                % constr = [constr, T_PBD{i} == T_PBDCAF{i} + S_PBDCAF];
                % T_PCA{i} = T_PCABDF{i} + S_PCABDF;
                % constr = [constr, T_PCA{i} == T_PCADBF{i} + S_PCADBF];
                % T_PCB{i} = T_PCBADF{i} + S_PCBADF;
                % constr = [constr, T_PCB{i} == T_PCBDAF{i} + S_PCBDAF];
                % T_PCD{i} = T_PCDABF{i} + S_PCDABF;
                % constr = [constr, T_PCD{i} == T_PCDBAF{i} + S_PCDBAF];
                % T_PDA{i} = T_PDABCF{i} + S_PDABCF;
                % constr = [constr, T_PDA{i} == T_PDACBF{i} + S_PDACBF];
                % T_PDB{i} = T_PDBACF{i} + S_PDBACF;
                % constr = [constr, T_PDB{i} == T_PDBCAF{i} + S_PDBCAF];
                % T_PDC{i} = T_PDCABF{i} + S_PDCABF;
                % constr = [constr, T_PDC{i} == T_PDCBAF{i} + S_PDCBAF];
                %
                % T_PA{i} = T_PAB{i} + S_PAB;
                % constr = [constr, T_PA{i} == T_PAC{i} + S_PAC, T_PA{i} == T_PAD{i} + S_PAD];
                % T_PB{i} = T_PBA{i} + S_PBA;
                % constr = [constr, T_PB{i} == T_PBC{i} + S_PBC, T_PB{i} == T_PBD{i} + S_PBD];
                % T_PC{i} = T_PCA{i} + S_PCA;
                % constr = [constr, T_PC{i} == T_PCB{i} + S_PCB, T_PC{i} == T_PCD{i} + S_PCD];
                % T_PD{i} = T_PDA{i} + S_PDA;
                % constr = [constr, T_PD{i} == T_PDB{i} + S_PDB, T_PD{i} == T_PDC{i} + S_PDC];
                %
                % T_P{i} = T_PA{i} + S_PA;
                % constr = [constr, T_P{i} == T_PB{i} + S_PB, T_P{i} == T_PC{i} + S_PC, T_P{i} == T_PD{i} + S_PD];
                %
                % constr = [constr, Sr{i} == T_P{i} + S_P];
                
                T_P{i} = Sr{i} - S_P;
                
                T_PA{i} = T_P{i} - S_PA;
                T_PB{i} = T_P{i} - S_PB;
                T_PC{i} = T_P{i} - S_PC;
                T_PD{i} = T_P{i} - S_PD;
                
                T_PAB{i} = T_PA{i} - S_PAB;
                T_PAC{i} = T_PA{i} - S_PAC;
                T_PAD{i} = T_PA{i} - S_PAD;
                
                T_PBA{i} = T_PB{i} - S_PBA;
                T_PBC{i} = T_PB{i} - S_PBC;
                T_PBD{i} = T_PB{i} - S_PBD;
                
                T_PCA{i} = T_PC{i} - S_PCA;
                T_PCB{i} = T_PC{i} - S_PCB;
                T_PCD{i} = T_PC{i} - S_PCD;
                
                T_PDA{i} = T_PD{i} - S_PDA;
                T_PDB{i} = T_PD{i} - S_PDB;
                T_PDC{i} = T_PD{i} - S_PDC;
                
                T_PABCDF{i} = T_PAB{i} - S_PABCDF;
                T_PABDCF{i} = T_PAB{i} - S_PABDCF;
                T_PACBDF{i} = T_PAC{i} - S_PACBDF;
                T_PACDBF{i} = T_PAC{i} - S_PACDBF;
                T_PADBCF{i} = T_PAD{i} - S_PADBCF;
                T_PADCBF{i} = T_PAD{i} - S_PADCBF;
                
                T_PBACDF{i} = T_PBA{i} - S_PBACDF;
                T_PBADCF{i} = T_PBA{i} - S_PBADCF;
                T_PBCADF{i} = T_PBC{i} - S_PBCADF;
                T_PBCDAF{i} = T_PBC{i} - S_PBCDAF;
                T_PBDACF{i} = T_PBD{i} - S_PBDACF;
                T_PBDCAF{i} = T_PBD{i} - S_PBDCAF;
                
                T_PCABDF{i} = T_PCA{i} - S_PCABDF;
                T_PCADBF{i} = T_PCA{i} - S_PCADBF;
                T_PCBADF{i} = T_PCB{i} - S_PCBADF;
                T_PCBDAF{i} = T_PCB{i} - S_PCBDAF; 
                T_PCDABF{i} = T_PCD{i} - S_PCDABF;
                T_PCDBAF{i} = T_PCD{i} - S_PCDBAF;

                T_PDABCF{i} = T_PDA{i} - S_PDABCF;
                T_PDACBF{i} = T_PDA{i} - S_PDACBF;
                T_PDBACF{i} = T_PDB{i} - S_PDBACF;
                T_PDBCAF{i} = T_PDB{i} - S_PDBCAF;
                T_PDCABF{i} = T_PDC{i} - S_PDCABF;
                T_PDCBAF{i} = T_PDC{i} - S_PDCBAF;
                
                constr = [constr, T_PABCDF{i} >= 0, T_PABDCF{i} >= 0, T_PACBDF{i} >= 0, T_PACDBF{i} >= 0, T_PADBCF{i} >= 0, T_PADCBF{i} >= 0, ...
                                  T_PBACDF{i} >= 0, T_PBADCF{i} >= 0, T_PBCADF{i} >= 0, T_PBCDAF{i} >= 0, T_PBDACF{i} >= 0, T_PBDCAF{i} >= 0, ...
                                  T_PCABDF{i} >= 0, T_PCADBF{i} >= 0, T_PCBADF{i} >= 0, T_PCBDAF{i} >= 0, T_PCDABF{i} >= 0, T_PCDBAF{i} >= 0, ...
                                  T_PDABCF{i} >= 0, T_PDACBF{i} >= 0, T_PDBACF{i} >= 0, T_PDBCAF{i} >= 0, T_PDCABF{i} >= 0, T_PDCBAF{i} >= 0, ];
            end
            
            S_PABCDF_proj = S_PABCDF - ( trRep(S_PABCDF,F,dims) - trRep(S_PABCDF,[DO,F],dims) );
            S_PABCDF_proj = S_PABCDF_proj - ( trRep(S_PABCDF_proj,[DI,DO,F],dims) - trRep(S_PABCDF_proj,[CO,DI,DO,F],dims) );
            S_PABDCF_proj = S_PABDCF - ( trRep(S_PABDCF,F,dims) - trRep(S_PABDCF,[CO,F],dims) );
            S_PABDCF_proj = S_PABDCF_proj - ( trRep(S_PABDCF_proj,[CI,CO,F],dims) - trRep(S_PABDCF_proj,[CI,CO,DO,F],dims));
            S_PACBDF_proj = S_PACBDF - ( trRep(S_PACBDF,F,dims) - trRep(S_PACBDF,[DO,F],dims) );
            S_PACBDF_proj = S_PACBDF_proj - ( trRep(S_PACBDF_proj,[DI,DO,F],dims) - trRep(S_PACBDF_proj,[BO,DI,DO,F],dims) );
            S_PACDBF_proj = S_PACDBF - ( trRep(S_PACDBF,F,dims) - trRep(S_PACDBF,[BO,F],dims) );
            S_PACDBF_proj = S_PACDBF_proj - ( trRep(S_PACDBF_proj,[BI,BO,F],dims) - trRep(S_PACDBF_proj,[BI,BO,DO,F],dims) );
            S_PADBCF_proj = S_PADBCF - ( trRep(S_PADBCF,F,dims) - trRep(S_PADBCF,[CO,F],dims) );
            S_PADBCF_proj = S_PADBCF_proj - ( trRep(S_PADBCF_proj,[CI,CO,F],dims) - trRep(S_PADBCF_proj,[BO,CI,CO,F],dims) );
            S_PADCBF_proj = S_PADCBF - ( trRep(S_PADCBF,F,dims) - trRep(S_PADCBF,[BO,F],dims) );
            S_PADCBF_proj = S_PADCBF_proj - ( trRep(S_PADCBF_proj,[BI,BO,F],dims) - trRep(S_PADCBF_proj,[BI,BO,CO,F],dims) );
            S_PBACDF_proj = S_PBACDF - ( trRep(S_PBACDF,F,dims) - trRep(S_PBACDF,[DO,F],dims) );
            S_PBACDF_proj = S_PBACDF_proj - ( trRep(S_PBACDF_proj,[DI,DO,F],dims) - trRep(S_PBACDF_proj,[CO,DI,DO,F],dims) );
            S_PBADCF_proj = S_PBADCF - ( trRep(S_PBADCF,F,dims) - trRep(S_PBADCF,[CO,F],dims) );
            S_PBADCF_proj = S_PBADCF_proj - ( trRep(S_PBADCF_proj,[CI,CO,F],dims) - trRep(S_PBADCF_proj,[CI,CO,DO,F],dims) );
            S_PBCADF_proj = S_PBCADF - ( trRep(S_PBCADF,F,dims) - trRep(S_PBCADF,[DO,F],dims) );
            S_PBCADF_proj = S_PBCADF_proj - ( trRep(S_PBCADF_proj,[DI,DO,F],dims) - trRep(S_PBCADF_proj,[AO,DI,DO,F],dims) );
            S_PBCDAF_proj = S_PBCDAF - ( trRep(S_PBCDAF,F,dims) - trRep(S_PBCDAF,[AO,F],dims) );
            S_PBCDAF_proj = S_PBCDAF_proj - ( trRep(S_PBCDAF_proj,[AI,AO,F],dims) - trRep(S_PBCDAF_proj,[AI,AO,DO,F],dims) );
            S_PBDACF_proj = S_PBDACF - ( trRep(S_PBDACF,F,dims) - trRep(S_PBDACF,[CO,F],dims) );
            S_PBDACF_proj = S_PBDACF_proj - ( trRep(S_PBDACF_proj,[CI,CO,F],dims) - trRep(S_PBDACF_proj,[AO,CI,CO,F],dims) );
            S_PBDCAF_proj = S_PBDCAF - ( trRep(S_PBDCAF,F,dims) - trRep(S_PBDCAF,[AO,F],dims) );
            S_PBDCAF_proj = S_PBDCAF_proj - ( trRep(S_PBDCAF_proj,[AI,AO,F],dims) - trRep(S_PBDCAF_proj,[AI,AO,CO,F],dims) );
            S_PCABDF_proj = S_PCABDF - ( trRep(S_PCABDF,F,dims) - trRep(S_PCABDF,[DO,F],dims) );
            S_PCABDF_proj = S_PCABDF_proj - ( trRep(S_PCABDF_proj,[DI,DO,F],dims) - trRep(S_PCABDF_proj,[BO,DI,DO,F],dims) );
            S_PCADBF_proj = S_PCADBF - ( trRep(S_PCADBF,F,dims) - trRep(S_PCADBF,[BO,F],dims) );
            S_PCADBF_proj = S_PCADBF_proj - ( trRep(S_PCADBF_proj,[BI,BO,F],dims) - trRep(S_PCADBF_proj,[BI,BO,DO,F],dims) );
            S_PCBADF_proj = S_PCBADF - ( trRep(S_PCBADF,F,dims) - trRep(S_PCBADF,[DO,F],dims) );
            S_PCBADF_proj = S_PCBADF_proj - ( trRep(S_PCBADF_proj,[DI,DO,F],dims) - trRep(S_PCBADF_proj,[AO,DI,DO,F],dims) );
            S_PCBDAF_proj = S_PCBDAF - ( trRep(S_PCBDAF,F,dims) - trRep(S_PCBDAF,[AO,F],dims) );
            S_PCBDAF_proj = S_PCBDAF_proj - ( trRep(S_PCBDAF_proj,[AI,AO,F],dims) - trRep(S_PCBDAF_proj,[AI,AO,DO,F],dims) );
            S_PCDABF_proj = S_PCDABF - ( trRep(S_PCDABF,F,dims) - trRep(S_PCDABF,[BO,F],dims) );
            S_PCDABF_proj = S_PCDABF_proj - ( trRep(S_PCDABF_proj,[BI,BO,F],dims) - trRep(S_PCDABF_proj,[AO,BI,BO,F],dims) );
            S_PCDBAF_proj = S_PCDBAF - ( trRep(S_PCDBAF,F,dims) - trRep(S_PCDBAF,[AO,F],dims) );
            S_PCDBAF_proj = S_PCDBAF_proj - ( trRep(S_PCDBAF_proj,[AI,AO,F],dims) - trRep(S_PCDBAF_proj,[AI,AO,BO,F],dims) );
            S_PDABCF_proj = S_PDABCF - ( trRep(S_PDABCF,F,dims) - trRep(S_PDABCF,[CO,F],dims) );
            S_PDABCF_proj = S_PDABCF_proj - ( trRep(S_PDABCF_proj,[CI,CO,F],dims) - trRep(S_PDABCF_proj,[BO,CI,CO,F],dims) );
            S_PDACBF_proj = S_PDACBF - ( trRep(S_PDACBF,F,dims) - trRep(S_PDACBF,[BO,F],dims) );
            S_PDACBF_proj = S_PDACBF_proj - ( trRep(S_PDACBF_proj,[BI,BO,F],dims) - trRep(S_PDACBF_proj,[BI,BO,CO,F],dims) );
            S_PDBACF_proj = S_PDBACF - ( trRep(S_PDBACF,F,dims) - trRep(S_PDBACF,[CO,F],dims) );
            S_PDBACF_proj = S_PDBACF_proj - ( trRep(S_PDBACF_proj,[CI,CO,F],dims) - trRep(S_PDBACF_proj,[AO,CI,CO,F],dims) );
            S_PDBCAF_proj = S_PDBCAF - ( trRep(S_PDBCAF,F,dims) - trRep(S_PDBCAF,[AO,F],dims) );
            S_PDBCAF_proj = S_PDBCAF_proj - ( trRep(S_PDBCAF_proj,[AI,AO,F],dims) - trRep(S_PDBCAF_proj,[AI,AO,CO,F],dims) );
            S_PDCABF_proj = S_PDCABF - ( trRep(S_PDCABF,F,dims) - trRep(S_PDCABF,[BO,F],dims) );
            S_PDCABF_proj = S_PDCABF_proj - ( trRep(S_PDCABF_proj,[BI,BO,F],dims) - trRep(S_PDCABF_proj,[AO,BI,BO,F],dims) );
            S_PDCBAF_proj = S_PDCBAF - ( trRep(S_PDCBAF,F,dims) - trRep(S_PDCBAF,[AO,F],dims) );
            S_PDCBAF_proj = S_PDCBAF_proj - ( trRep(S_PDCBAF_proj,[AI,AO,F],dims) - trRep(S_PDCBAF_proj,[AI,AO,BO,F],dims) );
            
            constr = [constr, S_PABCDF_proj == 0, S_PABDCF_proj == 0, S_PACBDF_proj == 0, S_PACDBF_proj == 0, S_PADBCF_proj == 0, S_PADCBF_proj == 0, ...
                              S_PBACDF_proj == 0, S_PBADCF_proj == 0, S_PBCADF_proj == 0, S_PBCDAF_proj == 0, S_PBDACF_proj == 0, S_PBDCAF_proj == 0, ...
                              S_PCABDF_proj == 0, S_PCADBF_proj == 0, S_PCBADF_proj == 0, S_PCBDAF_proj == 0, S_PCDABF_proj == 0, S_PCDBAF_proj == 0, ...
                              S_PDABCF_proj == 0, S_PDACBF_proj == 0, S_PDBACF_proj == 0, S_PDBCAF_proj == 0, S_PDCABF_proj == 0, S_PDCBAF_proj == 0];
            
            constr = [constr, S_PAB - ( trRep(S_PAB,[CI,CO,DI,DO,F],dims) - trRep(S_PAB,[BO,CI,CO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PAC - ( trRep(S_PAC,[BI,BO,DI,DO,F],dims) - trRep(S_PAC,[BI,BO,CO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PAD - ( trRep(S_PAD,[BI,BO,CI,CO,F],dims) - trRep(S_PAD,[BI,BO,CI,CO,DO,F],dims) ) == 0];
            constr = [constr, S_PBA - ( trRep(S_PBA,[CI,CO,DI,DO,F],dims) - trRep(S_PBA,[AO,CI,CO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PBC - ( trRep(S_PBC,[AI,AO,DI,DO,F],dims) - trRep(S_PBC,[AI,AO,CO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PBD - ( trRep(S_PBD,[AI,AO,CI,CO,F],dims) - trRep(S_PBD,[AI,AO,CI,CO,DO,F],dims) ) == 0];
            constr = [constr, S_PCA - ( trRep(S_PCA,[BI,BO,DI,DO,F],dims) - trRep(S_PCA,[AO,BI,BO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PCB - ( trRep(S_PCB,[AI,AO,DI,DO,F],dims) - trRep(S_PCB,[AI,AO,BO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PCD - ( trRep(S_PCD,[AI,AO,BI,BO,F],dims) - trRep(S_PCD,[AI,AO,BI,BO,DO,F],dims) ) == 0];
            constr = [constr, S_PDA - ( trRep(S_PDA,[BI,BO,CI,CO,F],dims) - trRep(S_PDA,[AO,BI,BO,CI,CO,F],dims) ) == 0];
            constr = [constr, S_PDB - ( trRep(S_PDB,[AI,AO,CI,CO,F],dims) - trRep(S_PDB,[AI,AO,BO,CI,CO,F],dims) ) == 0];
            constr = [constr, S_PDC - ( trRep(S_PDC,[AI,AO,BI,BO,F],dims) - trRep(S_PDC,[AI,AO,BI,BO,CO,F],dims) ) == 0];
            
            constr = [constr, S_PA - ( trRep(S_PA,[BI,BO,CI,CO,DI,DO,F],dims) - trRep(S_PA,[AO,BI,BO,CI,CO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PB - ( trRep(S_PB,[AI,AO,CI,CO,DI,DO,F],dims) - trRep(S_PB,[AI,AO,BO,CI,CO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PC - ( trRep(S_PC,[AI,AO,BI,BO,DI,DO,F],dims) - trRep(S_PC,[AI,AO,BI,BO,CO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PD - ( trRep(S_PD,[AI,AO,BI,BO,CI,CO,F],dims) - trRep(S_PD,[AI,AO,BI,BO,CI,CO,DO,F],dims) ) == 0];
            
            if d_P ~= 1
                constr = [constr, S_P - ( trRep(S_P,[AI,AO,BI,BO,CI,CO,DI,DO,F],dims) - trRep(S_P,[P,AI,AO,BI,BO,CI,CO,DI,DO,F],dims) ) == 0];
            end
            
        otherwise
            error('Currently only implemented up to N=4');
    end
    
    cone_constraints = constr;
end

