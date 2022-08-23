function [W,dims,parties] = OCB_process_matrix(eta1,eta2)
%OCB_process_matrix Generates the bipartite OCB process matrix and a generalised family of similar processes
%   [W,dims,parties] = OCB_process_matrix()
%   [W,dims,parties] = OCB_process_matrix(eta1,eta2)
%   
%   Ordering of systems is taken to be AI AO BI BO
%   The OCB matrix is described in O. Oreshkov, F. Costa, C. Brukner, Nature Commun 3, 1092 (2012)
%   and defined as W = 1/4*[ID + 1/sqrt(2)*(1.sz.sz.1 + sz.1.sx.sz)]
%
%   A generalised family of similar processes was described in C. Branciard, Sci. Rep. 6, 26018 (2016)
%   as W = 1/4*[ID + eta1*1.sz.sz.1 + eta2*sz.1.sx.sz]. W is valid if eta1^2 + eta2^2 <= 1.
%   The OCB process is recovered for eta1=eta2=1/sqrt(2), the default parameters here

% Written by Alastair Abbott, last modified 16 August 2022

    %% 

    if nargin < 2
        eta1 = 1/sqrt(2);
        eta2 = 1/sqrt(2);
    end
    assert(eta1^2 + eta2^2 <= 1, 'Error: inputs must satisfy eta1^2 + eta2^2 <= 1')

    dims = [2,2,2,2];
    parties = {{[]},{1,2},{3,4},{[]}};

    id = eye(2);
    sx = [0,1;1,0];
    sz = [1,0;0,-1];

    W = 1/4*(eye(16) + eta1*Tensor(id,sz,sz,id) + eta2*Tensor(sz,id,sx,sz));
end