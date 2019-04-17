function [array] = assemble_fom_rhs( param, fem_specifics, varargin )
% Assemble fom matrix for elliptic scalar problems
% input=
%           param: vector of parameters
%           fem_specifics: struct containing the information to build the
%           mesh and the fespace
% output=
%           array: struct containing the stiffness matrix in COO format

%     current_model = fem_specifics.model;
% 
%     if ( ~strcmp(current_model, 'navier_stokes') )
%         
%         if nargin == 2
%             array = assemble_elliptic_fom_rhs( param, fem_specifics );
%         end
%         if nargin == 4
%             array = assemble_elliptic_fom_rhs( param, fem_specifics, varargin{1}, varargin{2} );
%         end
%     end

    disp('Assemble FOM rhs from MATLAB pyorb interface')

    fom_problem = initialize_fom_simulation( fem_specifics );
    
    if nargin < 4
        array = fom_problem.assemble_fom_rhs( param );
    end
    if nargin >= 4
        array = fom_problem.assemble_fom_rhs( param, varargin{1}, varargin{2} );
    end
    
    
    
%     if ( strcmp(current_model, 'navier_stokes') )
%         if  nargin == 2
%             disp('Computing FULL FEM NS rhs')
%             array = assemble_ns_fom_rhs( param, fem_specifics );
%         end
%         if nargin == 4
%             disp('Calling assemble_ns_fom_rhs with varargin')
%             array = assemble_ns_fom_rhs( param, fem_specifics, varargin{1}, varargin{2} );
%         end
%     end

end

