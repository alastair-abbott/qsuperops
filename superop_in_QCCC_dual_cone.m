function cone_constraints = superop_in_QCCC_dual_cone(Sr, dims, parties)
%superop_in_QCCC_dual_cone Yalmip constraints for a set of matrices to be in the cone of witnesses for QC-CCs
%   cone_constraints = superop_in_QCCC_dual_cone(Sr, dims, parties) 
%   Sr can be either a single witness S, or a superinstrument witness Sr
%   Returns the yalmip constraints, i.e. for the dual cone
%
% Requires QETLAB for PermuteSystems, PartialTrace

% Written by Alastair Abbott 2021, last modified 16 August 2022

    % First put Sr in canonical ordering (this checks the input validity too)
    % The spaces P,AI,AO,...,F then correspond to dims 1,2,3,...,2*N+2
    if exist('parties','var') && ~isempty(parties)
        [Sr, dims, parties] = superop_to_canonical_ordering(Sr, dims, parties);
    else
        [Sr, dims, parties] = superop_to_canonical_ordering(Sr, dims);
    end

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
            for r = 1:R
                % %% Full versions of constraints to aid readability, before eliminating equality constraints
                % T_PAF{r} = sdpvar(d,d,'hermitian','complex');
                %
                % T_P{r} = T_PAF{r} + S_PAF;
                %
                % constr = [constr, Sr{r} == T_P{r} + S_P];
                
                T_P{r} = Sr{r} - S_P;
                T_PAF{r} = T_P{r} - S_PAF;
            end
            constr = [constr, superop_in_PSD_cone(T_PAF)];
            
            S_PAF_proj = S_PAF - (tr_replace(S_PAF,F,dims) - tr_replace(S_PAF,[AO,F],dims));
            
            constr = [constr, S_PAF_proj == 0];
            
            if d_P ~= 1
                constr = [constr, S_P - (tr_replace(S_P,[AI,AO,F],dims) - tr_replace(S_P,[P,AI,AO,F],dims)) == 0];
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
            for r = 1:R
                % %% Full versions of constraints to aid readability, before eliminating equality constraints
                % T_PABF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PBAF{r} = sdpvar(d,d,'hermitian','complex');
                %
                % T_P{r} = T_PABF{r} + S_PABF;
                % constr = [constr, T_P{r} == T_PBAF{r} + S_PBAF];
                %
                % constr = [constr, Sr{r} == T_P{r} + S_P];
                
                T_P{r} = Sr{r} - S_P;
                T_PABF{r} = T_P{r} - S_PABF;
                T_PBAF{r} = T_P{r} - S_PBAF;
            end
            constr = [constr, superop_in_PSD_cone(T_PABF), superop_in_PSD_cone(T_PBAF)];
            
            S_PABF_proj = S_PABF - (tr_replace(S_PABF,F,dims) - tr_replace(S_PABF,[BO,F],dims));
            S_PABF_proj = S_PABF_proj - (tr_replace(S_PABF_proj,[BI,BO,F],dims) - tr_replace(S_PABF_proj,[AO,BI,BO,F],dims));
            S_PBAF_proj = S_PBAF - (tr_replace(S_PBAF,F,dims) - tr_replace(S_PBAF,[AO,F],dims));
            S_PBAF_proj = S_PBAF_proj - (tr_replace(S_PBAF_proj,[AI,AO,F],dims) - tr_replace(S_PBAF_proj,[AI,AO,BO,F],dims));
            
            constr = [constr, S_PABF_proj == 0, S_PBAF_proj == 0];
            
            if d_P ~= 1
                constr = [constr, S_P - (tr_replace(S_P,[AI,AO,BI,BO,F],dims) - tr_replace(S_P,[P,AI,AO,BI,BO,F],dims)) == 0];
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
            for r = 1:R
                % %% Full versions of constraints to aid readability, before eliminating equality constraints
                % T_PABCF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PACBF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PBACF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PBCAF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PCABF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PCBAF{r} = sdpvar(d,d,'hermitian','complex');
                %
                % T_PA{r} = T_PABCF{r} + S_PABCF;
                % constr = [constr, T_PA{r} == T_PACBF{r} + S_PACBF];
                % T_PB{r} = T_PBACF{r} + S_PBACF;
                % constr = [constr, T_PB{r} == T_PBCAF{r} + S_PBCAF];
                % T_PC{r} = T_PCABF{r} + S_PCABF;
                % constr = [constr, T_PC{r} == T_PCBAF{r} + S_PCBAF];
                %
                % T_P{r} = T_PA{r} + S_PA;
                % constr = [constr, T_P{r} == T_PB{r} + S_PB, T_P{r} == T_PC{r} + S_PC];
                %
                % constr = [constr, Sr{r} == T_P{r} + S_P];
                
                T_P{r} = Sr{r} - S_P;
                
                T_PA{r} = T_P{r} - S_PA;
                T_PB{r} = T_P{r} - S_PB;
                T_PC{r} = T_P{r} - S_PC;
                
                T_PABCF{r} = T_PA{r} - S_PABCF;
                T_PACBF{r} = T_PA{r} - S_PACBF;

                T_PBACF{r} = T_PB{r} - S_PBACF;
                T_PBCAF{r} = T_PB{r} - S_PBCAF;
                
                T_PCABF{r} = T_PC{r} - S_PCABF;
                T_PCBAF{r} = T_PC{r} - S_PCBAF; 
            end
            constr = [constr, superop_in_PSD_cone(T_PABCF), superop_in_PSD_cone(T_PACBF),...
                              superop_in_PSD_cone(T_PBACF), superop_in_PSD_cone(T_PBCAF),...
                              superop_in_PSD_cone(T_PCABF), superop_in_PSD_cone(T_PCBAF)];
            
            S_PABCF_proj = S_PABCF - ( tr_replace(S_PABCF,F,dims) - tr_replace(S_PABCF,[CO,F],dims) );
            S_PABCF_proj = S_PABCF_proj - ( tr_replace(S_PABCF_proj,[CI,CO,F],dims) - tr_replace(S_PABCF_proj,[BO,CI,CO,F],dims) );
            S_PACBF_proj = S_PACBF - ( tr_replace(S_PACBF,F,dims) - tr_replace(S_PACBF,[BO,F],dims) );
            S_PACBF_proj = S_PACBF_proj - ( tr_replace(S_PACBF_proj,[BI,BO,F],dims) - tr_replace(S_PACBF_proj,[BI,BO,CO,F],dims) );
            S_PBACF_proj = S_PBACF - ( tr_replace(S_PBACF,F,dims) - tr_replace(S_PBACF,[CO,F],dims) );
            S_PBACF_proj = S_PBACF_proj - ( tr_replace(S_PBACF_proj,[CI,CO,F],dims) - tr_replace(S_PBACF_proj,[AO,CI,CO,F],dims) );
            S_PBCAF_proj = S_PBCAF - ( tr_replace(S_PBCAF,F,dims) - tr_replace(S_PBCAF,[AO,F],dims) );
            S_PBCAF_proj = S_PBCAF_proj - ( tr_replace(S_PBCAF_proj,[AI,AO,F],dims) - tr_replace(S_PBCAF_proj,[AI,AO,CO,F],dims) );
            S_PCABF_proj = S_PCABF - ( tr_replace(S_PCABF,F,dims) - tr_replace(S_PCABF,[BO,F],dims) );
            S_PCABF_proj = S_PCABF_proj - ( tr_replace(S_PCABF_proj,[BI,BO,F],dims) - tr_replace(S_PCABF_proj,[AO,BI,BO,F],dims) );
            S_PCBAF_proj = S_PCBAF - ( tr_replace(S_PCBAF,F,dims) - tr_replace(S_PCBAF,[AO,F],dims) );
            S_PCBAF_proj = S_PCBAF_proj - ( tr_replace(S_PCBAF_proj,[AI,AO,F],dims) - tr_replace(S_PCBAF_proj,[AI,AO,BO,F],dims) );
              
            constr = [constr, S_PABCF_proj == 0, S_PACBF_proj == 0, S_PBACF_proj == 0, S_PBCAF_proj == 0, S_PCABF_proj == 0, S_PCBAF_proj == 0];
       
            constr = [constr, S_PA - ( tr_replace(S_PA,[BI,BO,CI,CO,F],dims) - tr_replace(S_PA,[AO,BI,BO,CI,CO,F],dims) ) == 0];
            constr = [constr, S_PB - ( tr_replace(S_PB,[AI,AO,CI,CO,F],dims) - tr_replace(S_PB,[AI,AO,BO,CI,CO,F],dims) ) == 0];
            constr = [constr, S_PC - ( tr_replace(S_PC,[AI,AO,BI,BO,F],dims) - tr_replace(S_PC,[AI,AO,BI,BO,CO,F],dims) ) == 0];
            
            if d_P ~= 1
                constr = [constr, S_P - ( tr_replace(S_P,[AI,AO,BI,BO,CI,CO,F],dims) - tr_replace(S_P,[P,AI,AO,BI,BO,CI,CO,F],dims) ) == 0];
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
            for r = 1:R
                % %% Full versions of constraints to aid readability, before eliminating equality constraints
                % T_PABCDF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PABDCF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PACBDF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PACDBF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PADBCF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PADCBF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PBACDF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PBADCF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PBCADF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PBCDAF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PBDACF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PBDCAF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PCABDF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PCADBF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PCBADF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PCBDAF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PCDABF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PCDBAF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PDABCF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PDACBF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PDBACF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PDBCAF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PDCABF{r} = sdpvar(d,d,'hermitian','complex');
                % T_PDCBAF{r} = sdpvar(d,d,'hermitian','complex');
                %
                % T_PAB{r} = T_PABCDF{r} + S_PABCDF;
                % constr = [constr, T_PAB{r} == T_PABDCF{r} + S_PABDCF];
                % T_PAC{r} = T_PACBDF{r} + S_PACBDF;
                % constr = [constr, T_PAC{r} == T_PACDBF{r} + S_PACDBF];
                % T_PAD{r} = T_PADBCF{r} + S_PADBCF;
                % constr = [constr, T_PAD{r} == T_PADCBF{r} + S_PADCBF];
                % T_PBA{r} = T_PBACDF{r} + S_PBACDF;
                % constr = [constr, T_PBA{r} == T_PBADCF{r} + S_PBADCF];
                % T_PBC{r} = T_PBCADF{r} + S_PBCADF;
                % constr = [constr, T_PBC{r} == T_PBCDAF{r} + S_PBCDAF];
                % T_PBD{r} = T_PBDACF{r} + S_PBDACF;
                % constr = [constr, T_PBD{r} == T_PBDCAF{r} + S_PBDCAF];
                % T_PCA{r} = T_PCABDF{r} + S_PCABDF;
                % constr = [constr, T_PCA{r} == T_PCADBF{r} + S_PCADBF];
                % T_PCB{r} = T_PCBADF{r} + S_PCBADF;
                % constr = [constr, T_PCB{r} == T_PCBDAF{r} + S_PCBDAF];
                % T_PCD{r} = T_PCDABF{r} + S_PCDABF;
                % constr = [constr, T_PCD{r} == T_PCDBAF{r} + S_PCDBAF];
                % T_PDA{r} = T_PDABCF{r} + S_PDABCF;
                % constr = [constr, T_PDA{r} == T_PDACBF{r} + S_PDACBF];
                % T_PDB{r} = T_PDBACF{r} + S_PDBACF;
                % constr = [constr, T_PDB{r} == T_PDBCAF{r} + S_PDBCAF];
                % T_PDC{r} = T_PDCABF{r} + S_PDCABF;
                % constr = [constr, T_PDC{r} == T_PDCBAF{r} + S_PDCBAF];
                %
                % T_PA{r} = T_PAB{r} + S_PAB;
                % constr = [constr, T_PA{r} == T_PAC{r} + S_PAC, T_PA{r} == T_PAD{r} + S_PAD];
                % T_PB{r} = T_PBA{r} + S_PBA;
                % constr = [constr, T_PB{r} == T_PBC{r} + S_PBC, T_PB{r} == T_PBD{r} + S_PBD];
                % T_PC{r} = T_PCA{r} + S_PCA;
                % constr = [constr, T_PC{r} == T_PCB{r} + S_PCB, T_PC{r} == T_PCD{r} + S_PCD];
                % T_PD{r} = T_PDA{r} + S_PDA;
                % constr = [constr, T_PD{r} == T_PDB{r} + S_PDB, T_PD{r} == T_PDC{r} + S_PDC];
                %
                % T_P{r} = T_PA{r} + S_PA;
                % constr = [constr, T_P{r} == T_PB{r} + S_PB, T_P{r} == T_PC{r} + S_PC, T_P{r} == T_PD{r} + S_PD];
                %
                % constr = [constr, Sr{r} == T_P{r} + S_P];
                
                T_P{r} = Sr{r} - S_P;
                
                T_PA{r} = T_P{r} - S_PA;
                T_PB{r} = T_P{r} - S_PB;
                T_PC{r} = T_P{r} - S_PC;
                T_PD{r} = T_P{r} - S_PD;
                
                T_PAB{r} = T_PA{r} - S_PAB;
                T_PAC{r} = T_PA{r} - S_PAC;
                T_PAD{r} = T_PA{r} - S_PAD;
                
                T_PBA{r} = T_PB{r} - S_PBA;
                T_PBC{r} = T_PB{r} - S_PBC;
                T_PBD{r} = T_PB{r} - S_PBD;
                
                T_PCA{r} = T_PC{r} - S_PCA;
                T_PCB{r} = T_PC{r} - S_PCB;
                T_PCD{r} = T_PC{r} - S_PCD;
                
                T_PDA{r} = T_PD{r} - S_PDA;
                T_PDB{r} = T_PD{r} - S_PDB;
                T_PDC{r} = T_PD{r} - S_PDC;
                
                T_PABCDF{r} = T_PAB{r} - S_PABCDF;
                T_PABDCF{r} = T_PAB{r} - S_PABDCF;
                T_PACBDF{r} = T_PAC{r} - S_PACBDF;
                T_PACDBF{r} = T_PAC{r} - S_PACDBF;
                T_PADBCF{r} = T_PAD{r} - S_PADBCF;
                T_PADCBF{r} = T_PAD{r} - S_PADCBF;
                
                T_PBACDF{r} = T_PBA{r} - S_PBACDF;
                T_PBADCF{r} = T_PBA{r} - S_PBADCF;
                T_PBCADF{r} = T_PBC{r} - S_PBCADF;
                T_PBCDAF{r} = T_PBC{r} - S_PBCDAF;
                T_PBDACF{r} = T_PBD{r} - S_PBDACF;
                T_PBDCAF{r} = T_PBD{r} - S_PBDCAF;
                
                T_PCABDF{r} = T_PCA{r} - S_PCABDF;
                T_PCADBF{r} = T_PCA{r} - S_PCADBF;
                T_PCBADF{r} = T_PCB{r} - S_PCBADF;
                T_PCBDAF{r} = T_PCB{r} - S_PCBDAF; 
                T_PCDABF{r} = T_PCD{r} - S_PCDABF;
                T_PCDBAF{r} = T_PCD{r} - S_PCDBAF;

                T_PDABCF{r} = T_PDA{r} - S_PDABCF;
                T_PDACBF{r} = T_PDA{r} - S_PDACBF;
                T_PDBACF{r} = T_PDB{r} - S_PDBACF;
                T_PDBCAF{r} = T_PDB{r} - S_PDBCAF;
                T_PDCABF{r} = T_PDC{r} - S_PDCABF;
                T_PDCBAF{r} = T_PDC{r} - S_PDCBAF;
            end
            constr = [constr, superop_in_PSD_cone(T_PABCDF), superop_in_PSD_cone(T_PABDCF), superop_in_PSD_cone(T_PACBDF), superop_in_PSD_cone(T_PACDBF),...
                              superop_in_PSD_cone(T_PADBCF), superop_in_PSD_cone(T_PADCBF), superop_in_PSD_cone(T_PBACDF), superop_in_PSD_cone(T_PBADCF),...
                              superop_in_PSD_cone(T_PBADCF), superop_in_PSD_cone(T_PBCDAF), superop_in_PSD_cone(T_PBDACF), superop_in_PSD_cone(T_PBDCAF),...
                              superop_in_PSD_cone(T_PCABDF), superop_in_PSD_cone(T_PCADBF), superop_in_PSD_cone(T_PCBADF), superop_in_PSD_cone(T_PCBDAF),...
                              superop_in_PSD_cone(T_PCDABF), superop_in_PSD_cone(T_PCDBAF), superop_in_PSD_cone(T_PDABCF), superop_in_PSD_cone(T_PDACBF),...
                              superop_in_PSD_cone(T_PDBACF), superop_in_PSD_cone(T_PDBCAF), superop_in_PSD_cone(T_PDCABF), superop_in_PSD_cone(T_PDCBAF)];
            
            S_PABCDF_proj = S_PABCDF - ( tr_replace(S_PABCDF,F,dims) - tr_replace(S_PABCDF,[DO,F],dims) );
            S_PABCDF_proj = S_PABCDF_proj - ( tr_replace(S_PABCDF_proj,[DI,DO,F],dims) - tr_replace(S_PABCDF_proj,[CO,DI,DO,F],dims) );
            S_PABDCF_proj = S_PABDCF - ( tr_replace(S_PABDCF,F,dims) - tr_replace(S_PABDCF,[CO,F],dims) );
            S_PABDCF_proj = S_PABDCF_proj - ( tr_replace(S_PABDCF_proj,[CI,CO,F],dims) - tr_replace(S_PABDCF_proj,[CI,CO,DO,F],dims));
            S_PACBDF_proj = S_PACBDF - ( tr_replace(S_PACBDF,F,dims) - tr_replace(S_PACBDF,[DO,F],dims) );
            S_PACBDF_proj = S_PACBDF_proj - ( tr_replace(S_PACBDF_proj,[DI,DO,F],dims) - tr_replace(S_PACBDF_proj,[BO,DI,DO,F],dims) );
            S_PACDBF_proj = S_PACDBF - ( tr_replace(S_PACDBF,F,dims) - tr_replace(S_PACDBF,[BO,F],dims) );
            S_PACDBF_proj = S_PACDBF_proj - ( tr_replace(S_PACDBF_proj,[BI,BO,F],dims) - tr_replace(S_PACDBF_proj,[BI,BO,DO,F],dims) );
            S_PADBCF_proj = S_PADBCF - ( tr_replace(S_PADBCF,F,dims) - tr_replace(S_PADBCF,[CO,F],dims) );
            S_PADBCF_proj = S_PADBCF_proj - ( tr_replace(S_PADBCF_proj,[CI,CO,F],dims) - tr_replace(S_PADBCF_proj,[BO,CI,CO,F],dims) );
            S_PADCBF_proj = S_PADCBF - ( tr_replace(S_PADCBF,F,dims) - tr_replace(S_PADCBF,[BO,F],dims) );
            S_PADCBF_proj = S_PADCBF_proj - ( tr_replace(S_PADCBF_proj,[BI,BO,F],dims) - tr_replace(S_PADCBF_proj,[BI,BO,CO,F],dims) );
            S_PBACDF_proj = S_PBACDF - ( tr_replace(S_PBACDF,F,dims) - tr_replace(S_PBACDF,[DO,F],dims) );
            S_PBACDF_proj = S_PBACDF_proj - ( tr_replace(S_PBACDF_proj,[DI,DO,F],dims) - tr_replace(S_PBACDF_proj,[CO,DI,DO,F],dims) );
            S_PBADCF_proj = S_PBADCF - ( tr_replace(S_PBADCF,F,dims) - tr_replace(S_PBADCF,[CO,F],dims) );
            S_PBADCF_proj = S_PBADCF_proj - ( tr_replace(S_PBADCF_proj,[CI,CO,F],dims) - tr_replace(S_PBADCF_proj,[CI,CO,DO,F],dims) );
            S_PBCADF_proj = S_PBCADF - ( tr_replace(S_PBCADF,F,dims) - tr_replace(S_PBCADF,[DO,F],dims) );
            S_PBCADF_proj = S_PBCADF_proj - ( tr_replace(S_PBCADF_proj,[DI,DO,F],dims) - tr_replace(S_PBCADF_proj,[AO,DI,DO,F],dims) );
            S_PBCDAF_proj = S_PBCDAF - ( tr_replace(S_PBCDAF,F,dims) - tr_replace(S_PBCDAF,[AO,F],dims) );
            S_PBCDAF_proj = S_PBCDAF_proj - ( tr_replace(S_PBCDAF_proj,[AI,AO,F],dims) - tr_replace(S_PBCDAF_proj,[AI,AO,DO,F],dims) );
            S_PBDACF_proj = S_PBDACF - ( tr_replace(S_PBDACF,F,dims) - tr_replace(S_PBDACF,[CO,F],dims) );
            S_PBDACF_proj = S_PBDACF_proj - ( tr_replace(S_PBDACF_proj,[CI,CO,F],dims) - tr_replace(S_PBDACF_proj,[AO,CI,CO,F],dims) );
            S_PBDCAF_proj = S_PBDCAF - ( tr_replace(S_PBDCAF,F,dims) - tr_replace(S_PBDCAF,[AO,F],dims) );
            S_PBDCAF_proj = S_PBDCAF_proj - ( tr_replace(S_PBDCAF_proj,[AI,AO,F],dims) - tr_replace(S_PBDCAF_proj,[AI,AO,CO,F],dims) );
            S_PCABDF_proj = S_PCABDF - ( tr_replace(S_PCABDF,F,dims) - tr_replace(S_PCABDF,[DO,F],dims) );
            S_PCABDF_proj = S_PCABDF_proj - ( tr_replace(S_PCABDF_proj,[DI,DO,F],dims) - tr_replace(S_PCABDF_proj,[BO,DI,DO,F],dims) );
            S_PCADBF_proj = S_PCADBF - ( tr_replace(S_PCADBF,F,dims) - tr_replace(S_PCADBF,[BO,F],dims) );
            S_PCADBF_proj = S_PCADBF_proj - ( tr_replace(S_PCADBF_proj,[BI,BO,F],dims) - tr_replace(S_PCADBF_proj,[BI,BO,DO,F],dims) );
            S_PCBADF_proj = S_PCBADF - ( tr_replace(S_PCBADF,F,dims) - tr_replace(S_PCBADF,[DO,F],dims) );
            S_PCBADF_proj = S_PCBADF_proj - ( tr_replace(S_PCBADF_proj,[DI,DO,F],dims) - tr_replace(S_PCBADF_proj,[AO,DI,DO,F],dims) );
            S_PCBDAF_proj = S_PCBDAF - ( tr_replace(S_PCBDAF,F,dims) - tr_replace(S_PCBDAF,[AO,F],dims) );
            S_PCBDAF_proj = S_PCBDAF_proj - ( tr_replace(S_PCBDAF_proj,[AI,AO,F],dims) - tr_replace(S_PCBDAF_proj,[AI,AO,DO,F],dims) );
            S_PCDABF_proj = S_PCDABF - ( tr_replace(S_PCDABF,F,dims) - tr_replace(S_PCDABF,[BO,F],dims) );
            S_PCDABF_proj = S_PCDABF_proj - ( tr_replace(S_PCDABF_proj,[BI,BO,F],dims) - tr_replace(S_PCDABF_proj,[AO,BI,BO,F],dims) );
            S_PCDBAF_proj = S_PCDBAF - ( tr_replace(S_PCDBAF,F,dims) - tr_replace(S_PCDBAF,[AO,F],dims) );
            S_PCDBAF_proj = S_PCDBAF_proj - ( tr_replace(S_PCDBAF_proj,[AI,AO,F],dims) - tr_replace(S_PCDBAF_proj,[AI,AO,BO,F],dims) );
            S_PDABCF_proj = S_PDABCF - ( tr_replace(S_PDABCF,F,dims) - tr_replace(S_PDABCF,[CO,F],dims) );
            S_PDABCF_proj = S_PDABCF_proj - ( tr_replace(S_PDABCF_proj,[CI,CO,F],dims) - tr_replace(S_PDABCF_proj,[BO,CI,CO,F],dims) );
            S_PDACBF_proj = S_PDACBF - ( tr_replace(S_PDACBF,F,dims) - tr_replace(S_PDACBF,[BO,F],dims) );
            S_PDACBF_proj = S_PDACBF_proj - ( tr_replace(S_PDACBF_proj,[BI,BO,F],dims) - tr_replace(S_PDACBF_proj,[BI,BO,CO,F],dims) );
            S_PDBACF_proj = S_PDBACF - ( tr_replace(S_PDBACF,F,dims) - tr_replace(S_PDBACF,[CO,F],dims) );
            S_PDBACF_proj = S_PDBACF_proj - ( tr_replace(S_PDBACF_proj,[CI,CO,F],dims) - tr_replace(S_PDBACF_proj,[AO,CI,CO,F],dims) );
            S_PDBCAF_proj = S_PDBCAF - ( tr_replace(S_PDBCAF,F,dims) - tr_replace(S_PDBCAF,[AO,F],dims) );
            S_PDBCAF_proj = S_PDBCAF_proj - ( tr_replace(S_PDBCAF_proj,[AI,AO,F],dims) - tr_replace(S_PDBCAF_proj,[AI,AO,CO,F],dims) );
            S_PDCABF_proj = S_PDCABF - ( tr_replace(S_PDCABF,F,dims) - tr_replace(S_PDCABF,[BO,F],dims) );
            S_PDCABF_proj = S_PDCABF_proj - ( tr_replace(S_PDCABF_proj,[BI,BO,F],dims) - tr_replace(S_PDCABF_proj,[AO,BI,BO,F],dims) );
            S_PDCBAF_proj = S_PDCBAF - ( tr_replace(S_PDCBAF,F,dims) - tr_replace(S_PDCBAF,[AO,F],dims) );
            S_PDCBAF_proj = S_PDCBAF_proj - ( tr_replace(S_PDCBAF_proj,[AI,AO,F],dims) - tr_replace(S_PDCBAF_proj,[AI,AO,BO,F],dims) );
            
            constr = [constr, S_PABCDF_proj == 0, S_PABDCF_proj == 0, S_PACBDF_proj == 0, S_PACDBF_proj == 0, S_PADBCF_proj == 0, S_PADCBF_proj == 0, ...
                              S_PBACDF_proj == 0, S_PBADCF_proj == 0, S_PBCADF_proj == 0, S_PBCDAF_proj == 0, S_PBDACF_proj == 0, S_PBDCAF_proj == 0, ...
                              S_PCABDF_proj == 0, S_PCADBF_proj == 0, S_PCBADF_proj == 0, S_PCBDAF_proj == 0, S_PCDABF_proj == 0, S_PCDBAF_proj == 0, ...
                              S_PDABCF_proj == 0, S_PDACBF_proj == 0, S_PDBACF_proj == 0, S_PDBCAF_proj == 0, S_PDCABF_proj == 0, S_PDCBAF_proj == 0];
            
            constr = [constr, S_PAB - ( tr_replace(S_PAB,[CI,CO,DI,DO,F],dims) - tr_replace(S_PAB,[BO,CI,CO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PAC - ( tr_replace(S_PAC,[BI,BO,DI,DO,F],dims) - tr_replace(S_PAC,[BI,BO,CO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PAD - ( tr_replace(S_PAD,[BI,BO,CI,CO,F],dims) - tr_replace(S_PAD,[BI,BO,CI,CO,DO,F],dims) ) == 0];
            constr = [constr, S_PBA - ( tr_replace(S_PBA,[CI,CO,DI,DO,F],dims) - tr_replace(S_PBA,[AO,CI,CO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PBC - ( tr_replace(S_PBC,[AI,AO,DI,DO,F],dims) - tr_replace(S_PBC,[AI,AO,CO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PBD - ( tr_replace(S_PBD,[AI,AO,CI,CO,F],dims) - tr_replace(S_PBD,[AI,AO,CI,CO,DO,F],dims) ) == 0];
            constr = [constr, S_PCA - ( tr_replace(S_PCA,[BI,BO,DI,DO,F],dims) - tr_replace(S_PCA,[AO,BI,BO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PCB - ( tr_replace(S_PCB,[AI,AO,DI,DO,F],dims) - tr_replace(S_PCB,[AI,AO,BO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PCD - ( tr_replace(S_PCD,[AI,AO,BI,BO,F],dims) - tr_replace(S_PCD,[AI,AO,BI,BO,DO,F],dims) ) == 0];
            constr = [constr, S_PDA - ( tr_replace(S_PDA,[BI,BO,CI,CO,F],dims) - tr_replace(S_PDA,[AO,BI,BO,CI,CO,F],dims) ) == 0];
            constr = [constr, S_PDB - ( tr_replace(S_PDB,[AI,AO,CI,CO,F],dims) - tr_replace(S_PDB,[AI,AO,BO,CI,CO,F],dims) ) == 0];
            constr = [constr, S_PDC - ( tr_replace(S_PDC,[AI,AO,BI,BO,F],dims) - tr_replace(S_PDC,[AI,AO,BI,BO,CO,F],dims) ) == 0];
            
            constr = [constr, S_PA - ( tr_replace(S_PA,[BI,BO,CI,CO,DI,DO,F],dims) - tr_replace(S_PA,[AO,BI,BO,CI,CO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PB - ( tr_replace(S_PB,[AI,AO,CI,CO,DI,DO,F],dims) - tr_replace(S_PB,[AI,AO,BO,CI,CO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PC - ( tr_replace(S_PC,[AI,AO,BI,BO,DI,DO,F],dims) - tr_replace(S_PC,[AI,AO,BI,BO,CO,DI,DO,F],dims) ) == 0];
            constr = [constr, S_PD - ( tr_replace(S_PD,[AI,AO,BI,BO,CI,CO,F],dims) - tr_replace(S_PD,[AI,AO,BI,BO,CI,CO,DO,F],dims) ) == 0];
            
            if d_P ~= 1
                constr = [constr, S_P - ( tr_replace(S_P,[AI,AO,BI,BO,CI,CO,DI,DO,F],dims) - tr_replace(S_P,[P,AI,AO,BI,BO,CI,CO,DI,DO,F],dims) ) == 0];
            end
            
        otherwise
            error('Currently only implemented up to N=4');
    end
    
    cone_constraints = constr;
end

