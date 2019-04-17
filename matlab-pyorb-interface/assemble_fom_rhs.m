function [array] = assemble_fom_rhs( param, fem_specifics, varargin )
% Assemble fom matrix for elliptic scalar problems
% input=
%           param: vector of parameters
%           fem_specifics: struct containing the information to build the
%           mesh and the fespace
% output=
%           array: struct containing the stiffness matrix in COO format

    fom_problem = initialize_fom_simulation( fem_specifics );
    
    if nargin < 4
        array = fom_problem.assemble_fom_rhs( param );
    end
    if nargin >= 4
        array = fom_problem.assemble_fom_rhs( param, varargin{1}, varargin{2} );
    end

end

