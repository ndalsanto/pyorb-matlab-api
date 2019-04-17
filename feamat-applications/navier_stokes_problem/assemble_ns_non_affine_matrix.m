function [array] = assemble_ns_non_affine_matrix( navier_stokes_fom_problem, param, varargin )
% Assemble fom nonaffine part of Navier-Stokes matrix 
% input=
%           param: vector of parameters
%           fem_specifics: struct containing the information to build the
%           mesh and the fespace
%           varargin can contain the element_list of the subelements where 
%           we want to integrate as first argument and the indices of the
%           matrix we want to extract
% output=
%           array: struct containing the stiffness matrix in COO format

    n_nodes_u = size(navier_stokes_fom_problem.fespace_u.nodes,1);
    ccc = build_obstruction_function( );
    cc = @(x) ccc(x, param);

    disp('assemble_ns_non_affine_matrix --- assemble_ns_non_affine_matrix')
    disp(nargin)
    
    if nargin == 2
        M_scalar = assemble_cu( cc, navier_stokes_fom_problem.fespace_u );
        M_scalar = apply_dirichlet_bc_matrix(M_scalar, navier_stokes_fom_problem.fespace_u, 0);
        M = [M_scalar, sparse( n_nodes_u, n_nodes_u ); ...
             sparse( n_nodes_u, n_nodes_u ), M_scalar ];

        [ i, j, val ] = find( M );
        array.A = [ i, j, val ];
        
        return
    end

%     if nargin > 2
    disp( 'Accessing subelements ' )
%         element_list = varargin{1};
    M_scalar = assemble_cu_elementlist( cc, navier_stokes_fom_problem.fespace_u, varargin{1} );
    M_scalar = apply_dirichlet_bc_matrix(M_scalar, navier_stokes_fom_problem.fespace_u, 0);
    M = [M_scalar, sparse( n_nodes_u, n_nodes_u ); ...
         sparse( n_nodes_u, n_nodes_u ), M_scalar ];

    if nargin > 3
        indeces_list = varargin{2}; % supposed to be a matrix of size nb_of_indices x 2 
%         elements_A = zeros( size( indeces_list, 1 ), 1 );
%         for ii = 1:size( indeces_list , 1 )
%             elements_A(ii) = M(indeces_list(ii, 1), indeces_list(ii, 2));
%         end
%         array.A = [indeces_list, elements_A];

        disp( 'Accessing in new new way ' );
        array.A = [ indeces_list, full( ( M( sub2ind(size(M), indeces_list(:,1), indeces_list(:,2) ) ) ) ) ];

        return
    else
        [ i, j, val ] = find( M );
        array.A = [ i, j, val ];
    end
%     end
    
end







