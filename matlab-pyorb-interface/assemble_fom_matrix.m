function [array] = assemble_fom_matrix( param, fem_specifics, varargin )
% Assemble fom matrix for elliptic scalar problems
% input=
%           param: vector of parameters
%           fem_specifics: struct containing the information to build the
%           mesh and the fespace
% output=
%           array: struct containing the stiffness matrix in COO format

    disp('Assemble FOM matrix from MATLAB pyorb interface')

    fom_problem = initialize_fom_simulation( fem_specifics );
    
    if nargin < 4
        array = fom_problem.assemble_fom_matrix( param );
    end
    
    if nargin >= 4
        array = fom_problem.assemble_fom_matrix( param, varargin{1}, varargin{2} );
    end
    
%     if (strcmp( 'navier_stokes', fem_specifics.model ) )
%         
%         if nargin == 2
%             array = assemble_ns_non_affine_matrix( param, fem_specifics );
%         end
%         if nargin > 2
%             disp('Assembling non affine matrix on reduced mesh')
%             tic
%             array = assemble_ns_non_affine_matrix( param, fem_specifics, varargin{1}, varargin{2} );
%             toc
%         end
%     end
    
end

