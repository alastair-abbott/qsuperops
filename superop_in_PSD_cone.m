function cone_constraints = superop_in_PSD_cone(Wr,tol)
%superop_in_PSD_cone Determines whether Wr is a superinstrument of PSD operators, up to tolerance
%   cone_constraints = superop_in_PSD_cone(Wr, tol)
%   Wr can be either a superinstrument or a superoperator/process matrix
%   Depending on whether Wr is an sdpvar or not, returns yalmip constraints or Boolean
%   tol is an optional argument, with default value 1e-6, which only impacts numerical checks (not sdpvars)

% Written by Alastair Abbott 2022, last modified 18 August 2022

    % default tolerance
    if nargin < 2
        tol = 1e-6;
    end

    % Treat a process matrix as a 1-element superinstrument
    if ~iscell(Wr)
        Wr = {Wr};
    end
    
    R = length(Wr);

    % Keep the logical and yalmip constraints separate until the end
    constraints_logical = [true];
    constraints_yalmip = [true];

    for r = 1:R
        % Each element of superinstrument should be PSD
        if isa(Wr{r},'sdpvar')
            constraints_yalmip = [constraints_yalmip, Wr{r} >= 0];
        else
            constraints_logical = [constraints_logical, min(eig(Wr{r})) >= -tol];
        end
    end

    % Combine the two types of constraints
    constraints_logical = all(constraints_logical);
    if constraints_logical == false
        cone_constraints = false;
    else
        cone_constraints = constraints_yalmip;
    end

end

