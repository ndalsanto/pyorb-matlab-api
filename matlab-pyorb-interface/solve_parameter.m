function [sol] = solve_parameter( param, fem_specifics )
% Assemble fom matrix for elliptic scalar problems
% input=
%           param: vector of parameters
%           fem_specifics: struct containing the information to build the
%           mesh and the fespace
% output=
%           sol: struct containing the solution

    disp('Solving FOM problem from MATLAB pyorb interface')

    fom_problem = initialize_fom_simulation( fem_specifics );
    
    sol = fom_problem.solve_parameter( param );
    
end