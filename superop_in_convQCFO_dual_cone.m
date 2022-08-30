function cone_constraints = superop_in_convQCFO_witness_cone(Sr, dims, parties)
%superop_in_convQCFO_witness_cone Yalmip constraints for a set of matrices to be in the cone of witnesses for convex mixtures of QC-FOs
%   cone_constraints = superop_in_convQCFO_witness_cone(Sr, dims, parties) 
%   Sr can be either a single witness S, or a superinstrument witness Sr
%   Returns the yalmip constraints, i.e. for the dual cone
%
% Requires QETLAB for PermuteSystems, PartialTrace

% Written by Alastair Abbott 2022, last modified 30 August 2022

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
            %% Only one order possible, reuse QCFO code
            
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties)];
        case 2
            %% Recall we now have a conically ordered process and we've grouped the "sub"-spaces
            
            % Sr has to be in witness cone for both possible orders
            % We can just reinterpret the parties and check like that
            parties_BAF_permuted = parties([1,3,2,4]);

            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_BAF_permuted)];

        case 3
            %% Recall we now have a conically ordered process and we've grouped the "sub"-spaces
            
            % Sr has to be in witness cone for all possible orders
            % We can just reinterpret the parties and check like that
            parties_ACBF_permuted = parties([1,2,4,3,5]);
            parties_BACF_permuted = parties([1,3,2,4,5]);
            parties_BCAF_permuted = parties([1,3,4,2,5]);
            parties_CABF_permuted = parties([1,4,3,2,5]);
            parties_CBAF_permuted = parties([1,4,2,3,5]);

            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_ACBF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_BACF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_BCAF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_CABF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_CBAF_permuted)];
            
        case 4
            %% Recall we now have a conically ordered process and we've grouped the "sub"-spaces
            
            % Sr has to be in witness cone for all possible orders
            % We can just reinterpret the parties and check like that
            parties_ABDCF_permuted = parties([1,2,3,5,4,6]);
            parties_ACBDF_permuted = parties([1,2,4,3,5,6]);
            parties_ACDBF_permuted = parties([1,2,4,5,3,6]);
            parties_ADBCF_permuted = parties([1,2,5,3,4,6]);
            parties_ADCBF_permuted = parties([1,2,5,4,3,6]);
            parties_BACDF_permuted = parties([1,3,2,4,5,6]);
            parties_BADCF_permuted = parties([1,3,2,5,4,6]);
            parties_BCADF_permuted = parties([1,3,4,2,5,6]);
            parties_BCDAF_permuted = parties([1,3,4,5,2,6]);
            parties_BDACF_permuted = parties([1,3,5,2,4,6]);
            parties_BDCAF_permuted = parties([1,3,5,4,2,6]);
            parties_CABDF_permuted = parties([1,4,2,3,5,6]);
            parties_CADBF_permuted = parties([1,4,2,5,3,6]);
            parties_CBADF_permuted = parties([1,4,3,2,5,6]);
            parties_CBDAF_permuted = parties([1,4,3,5,2,6]);
            parties_CDABF_permuted = parties([1,4,5,2,3,6]);
            parties_CDBAF_permuted = parties([1,4,5,3,2,6]);
            parties_DABCF_permuted = parties([1,5,2,3,4,6]);
            parties_DACBF_permuted = parties([1,5,2,4,3,6]);
            parties_DBACF_permuted = parties([1,5,3,2,4,6]);
            parties_DBCAF_permuted = parties([1,5,3,4,2,6]);
            parties_DCABF_permuted = parties([1,5,4,2,3,6]);
            parties_DCBAF_permuted = parties([1,5,4,3,2,6]);

            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_ABDCF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_ACBDF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_ACDBF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_ADBCF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_ADCBF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_BACDF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_BADCF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_BCADF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_BCDAF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_BDACF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_BDCAF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_CABDF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_CADBF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_CBADF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_CBDAF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_CDABF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_CDBAF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_DABCF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_DACBF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_DBACF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_DBCAF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_DCABF_permuted)];
            constr = [constr, superop_in_QCFO_witness_cone(Sr,dims,parties_DCBAF_permuted)];
            
        otherwise
            error('Currently only implemented up to N=4');
    end
    
    cone_constraints = constr;
end

