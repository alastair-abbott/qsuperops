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

    % Initialise like this to ensure it stays a logical vector for numeric input
    cone_constraints = [true];

    for r = 1:R
        % Each element of superinstrument should be PSD
        if isa(Wr{r},'sdpvar')
            cone_constraints = [cone_constraints, Wr{r} >= 0];
        else
            cone_constraints = [cone_constraints, all(eig(Wr{r}) > -tol)];
        end
    end

    % If all the elements of Wr were numeric, we take logical AND.
    % If we had at least one sdpvar, then cone_constraints are yalmip constraints
    if isa(cone_constraints,'logical')
        cone_constraints = all(cone_constraints);
    end

end

