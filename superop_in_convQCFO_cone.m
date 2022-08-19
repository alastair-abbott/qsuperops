function cone_constraints = superop_in_convQCFO_cone(Wr, dims, parties)
%superop_in_convQCFO_cone Yalmip constraints for a superinstrument to be in the cone of convex mixtures of QC-FOs
%   cone_constraints = superop_in_convQCFO_cone(Wr, dims, parties)
%   Wr can be either a superinstrument or a superoperator/process matrix
%   Returns the yalmip SDP constraints for Wr to be in the convex hull of conees of valid Quantum Circuits with Fixed-Order
%   This is similar to QC-CC, but without the possibility for dynamical causal structures
%   Note: this doesn't enforce the normalisation of the superoperator

% Written by Alastair Abbott 2022, last modified 19 August 2022

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
            % In this case only one order is possible, so just reuse QC-FO code
            constr = superop_in_QCFO_cone(Wr,dims,parties);
        case 2
            % Recall we now have a conically ordered process and we've grouped the "sub"-spaces
      
            Wr_ABF = cell(1,R);
            Wr_BAF = cell(1,R);
            for r = 1:R
                Wr_ABF{r} = sdpvar(d,d,'hermitian','complex');
                      
                % Define last element from others rather than enforcing an equality constraints
                Wr_BAF{r} = Wr{r} - Wr_ABF{r}; 
            end
            
            % Wr_ABF anf Wr_BAF should be in QC-FO cone for respective order
            % We can just reinterpret the parties and check like that
            parties_BAF_permuted = parties([1,3,2,4]);

            constr = [constr, superop_in_QCFO_cone(Wr_ABF,dims,parties)];
            constr = [constr, superop_in_QCFO_cone(Wr_BAF,dims,parties_BAF_permuted)];

        case 3   
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
            
            % Each Wr_perm should be be in QC-FO cone for respective order
            parties_ACBF_permuted = parties([1,2,4,3,5]);
            parties_BACF_permuted = parties([1,3,2,4,5]);
            parties_BCAF_permuted = parties([1,3,4,2,5]);
            parties_CABF_permuted = parties([1,4,3,2,5]);
            parties_CBAF_permuted = parties([1,4,2,3,5]);

            constr = [constr, superop_in_QCFO_cone(Wr_ABCF,dims,parties)];
            constr = [constr, superop_in_QCFO_cone(Wr_ACBF,dims,parties_ACBF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_BACF,dims,parties_BACF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_BCAF,dims,parties_BCAF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_CABF,dims,parties_CABF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_CBAF,dims,parties_CBAF_permuted)];

        case 4
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
            
            % Each Wr_perm should be be in QC-FO cone for respective order
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

            constr = [constr, superop_in_QCFO_cone(Wr_ABCDF,dims,parties)];
            constr = [constr, superop_in_QCFO_cone(Wr_ABDCF,dims,parties_ABDCF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_ACBDF,dims,parties_ACBDF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_ACDBF,dims,parties_ACDBF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_ADBCF,dims,parties_ADBCF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_ADCBF,dims,parties_ADCBF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_BACDF,dims,parties_BACDF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_BADCF,dims,parties_BADCF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_BCADF,dims,parties_BCADF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_BCDAF,dims,parties_BCDAF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_BDACF,dims,parties_BDACF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_BDCAF,dims,parties_BDCAF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_CABDF,dims,parties_CABDF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_CADBF,dims,parties_CADBF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_CBADF,dims,parties_CBADF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_CBDAF,dims,parties_CBDAF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_CDABF,dims,parties_CDABF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_CDBAF,dims,parties_CDBAF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_DABCF,dims,parties_DABCF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_DACBF,dims,parties_DACBF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_DBACF,dims,parties_DBACF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_DBCAF,dims,parties_DBCAF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_DCABF,dims,parties_DCABF_permuted)];
            constr = [constr, superop_in_QCFO_cone(Wr_DCBAF,dims,parties_DCBAF_permuted)];

        otherwise
            error('Currently only implemented up to N=4');
    end
    cone_constraints = constr;
end

