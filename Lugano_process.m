function [W, dims, parties] = Lugano_process(purification)
%Lugano_process Generates the tripartite Lugano process 
%   [W, dims, parties] = Lugano_process()
%   [W, dims, parties] = Lugano_process(purification)
%
%   The Lugano process is described in A. Baumeler, S. Wolf, New J. Phys 18, 013036 (2016)
%   If the optional argument is set to true, the unitary purification of this classical process
%   as described in M. Araújo, A. Feix, M. Navascués, C. Brukner, Quantum 1, 10 (2017) is returned instead
%
% Requires QETLAB for PermuteSystems

% Written by Alastair Abbott, last modified 17 August 2022
 
    if ~exist('purification','var')
        purification = false;
    end

    % Define the Lugano "function"
    f = @(a,b,c)[~b & c, ~c & a, ~a & b];
    % easy way to reference basis kets
    ket = @(x)[1-x;x];

    if ~purification
        dims = 2*ones(1,6);
        d = prod(dims);
        % Note that it's easiest to express the process if we don't use canonical ordering
        AO = 1;
        BO = 2;
        CO = 3;
        AI = 4;
        BI = 5;
        CI = 6;

        parties = {{[]}, {AI,AO}, {BI,BO}, {CI,CO}, {[]}};

        W = zeros(d,d);
        for a = 0:1
            for b = 0:1
                for c = 0:1
                    f_abc = f(a,b,c);
                    w_term = Tensor(ket(a),ket(b),ket(c),ket(f_abc(1)),ket(f_abc(2)),ket(f_abc(3)));
                    W = W + pure_to_mixed(w_term);
                end
            end
        end

    else
        dims = 2*ones(1,12);
        d = prod(dims);
        % Note that it's easiest to express the process if we don't use canonical ordering
        P = [4,5,6];
        AO = 1;
        BO = 2;
        CO = 3;
        AI = 10;
        BI = 11;
        CI = 12;
        F = [7,8,9];

        parties = {{P}, {AI,AO}, {BI,BO}, {CI,CO}, {F}};

        w = zeros(d,1);
        for a = 0:1
            for b = 0:1
                for c = 0:1
                    for i = 0:1
                        for j = 0:1
                            for k = 0:1
                                f_abc = f(a,b,c);
                                w = w + Tensor(ket(a),ket(b),ket(c),ket(i),ket(j),ket(k),...
                                                ket(a),ket(b),ket(c),ket(xor(i,f_abc(1))),ket(xor(j,f_abc(2))),ket(xor(k,f_abc(3))));
                            end
                        end
                    end
                end
            end
        end
        W = pure_to_mixed(w);

    end
       
end
