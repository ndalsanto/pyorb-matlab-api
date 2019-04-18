function [sol] = solve_parameter( param, fem_specifics )
% Assemble fom matrix for elliptic scalar problems
% input=
%           param: vector of parameters
%           fem_specifics: struct containing the information to build the
%           mesh and the fespace
% output=
%           sol: struct containing the solution

    fom_problem = load_fom_problem( fem_specifics );    
    sol = fom_problem.solve_parameter( param );
    
end