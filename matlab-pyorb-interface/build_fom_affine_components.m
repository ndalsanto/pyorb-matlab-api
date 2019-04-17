function [array] = build_fom_affine_components( operator, fem_specifics )
% Assemble fom affine matrix for elliptic scalar problems
% input=
%           operator: operator corresponding to stiffness matrix or RHS 
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the chosen model
% output=
%           array: struct containing the affine stiffness matrices in COO format
    

    disp('Building affine FOM components from MATLAB pyorb interface')

    fom_problem = initialize_fom_simulation( fem_specifics );
    
    array = fom_problem.build_fom_affine_components( operator );

%     if ~strcmp( 'navier_stokes', considered_model )
% 
%     else
% 
%         [fespace_u, fespace_p] = set_ns_simulation( fem_specifics );
%         n_nodes_u = size( fespace_u.nodes, 1 );
%         n_nodes_p = size( fespace_p.nodes, 1 );
%         dirichlet_stokes_functions    = @(x) [0 0;0 0;0 0;0 0;0 0;0 0]';
%         neumann_functions             = @(x) [0 0;0 0;0 0;0 0;0 0;0 0]';
%         f = [0; 0];
% 
%         if operator == 'A'
% 
%             disp( strcat( 'Retrieving operator A affine components for Navier Stokes - ', num2str(n_nodes_u) ) );
%             
%             u_lifting_stokes = build_stokes_lifting( fespace_u, fespace_p, fem_specifics );
%             u_lifting_stokes(1+2*n_nodes_u:end) = 0;
% 
%             C_1 = sparse( 2*n_nodes_u+n_nodes_p, 2*n_nodes_u+n_nodes_p );
%             C_1 = add_convective_term( C_1, u_lifting_stokes, fespace_u );
%             C_1 = add_flipped_convective_term( C_1, u_lifting_stokes, fespace_u );
% 
%             [i,j,val] = find( C_1(1:2*n_nodes_u, 1:2*n_nodes_u) );
%             array.A0 = [i,j,val];
% 
%             [A_stokes, ~] = assembler_steady_stokes( fespace_u, fespace_p, f, 1., dirichlet_stokes_functions, ...
%                                                      neumann_functions);
% 
%             [i,j,val] = find( A_stokes(1:2*n_nodes_u, 1:2*n_nodes_u) );
%             array.A1 = [i,j,val];
% 
%             for iB = 1:fem_specifics.range_rb_functions
%                 
%                 disp( strcat('Considering RB function - ', num2str(iB) ) );
% 
%                 pyorb_u = [fem_specifics.("rb_func_" + num2str(iB-1))'; zeros(n_nodes_p, 1)];
%                 C_1 = sparse(2*n_nodes_u+n_nodes_p, 2*n_nodes_u+n_nodes_p);
%                 C_1 = add_convective_term( C_1, pyorb_u, fespace_u );
% 
%                 [i,j,val] = find( C_1(1:2*n_nodes_u, 1:2*n_nodes_u) );
%                 array.("A" + num2str(iB+1) ) = [i,j,val];
%             end
%             
%             % needed when building the Jacobian
%             for iB = 1:fem_specifics.range_rb_functions
%                 
%                 disp( strcat('Considering RB function with flipped term - ', num2str(iB) ) );
% 
%                 pyorb_u = [fem_specifics.("rb_func_" + num2str(iB-1))'; zeros(n_nodes_p, 1)];
%                 C_1 = sparse(2*n_nodes_u+n_nodes_p, 2*n_nodes_u+n_nodes_p);
%                 C_1 = add_flipped_convective_term( C_1, pyorb_u, fespace_u );
% 
%                 [i,j,val] = find( C_1(1:2*n_nodes_u, 1:2*n_nodes_u) );
%                 array.("A" + num2str(fem_specifics.range_rb_functions+iB+1) ) = [i,j,val];
%             end
% 
%             disp( 'Finished to build fom NS affine components ' );
% 
%             return
% 
%         end
%         
%         if operator == 'f'
%             u_lifting_stokes = build_stokes_lifting( fespace_u, fespace_p, fem_specifics );
%             u_lifting_stokes(1+2*n_nodes_u:end) = 0;
%             dirichlet_functions = @(x) [0 0;0 0;0 0;0 0;0 0;0 0]';
% 
% 
%             [A_stokes, b_stokes] = assembler_steady_stokes( fespace_u, fespace_p, f, 1.0, dirichlet_functions, ...
%                                                             neumann_functions );
% 
%             b = b_stokes - A_stokes * u_lifting_stokes;
% 
%             b(1:n_nodes_u) = apply_dirichlet_bc_rhs( b(1:n_nodes_u), fespace_u, @(x) zeros(1, 2)*dirichlet_functions(x) );
%             b(n_nodes_u+1:2*n_nodes_u) = apply_dirichlet_bc_rhs( b(n_nodes_u+1:2*n_nodes_u), fespace_u, ... 
%                                                                  @(x) zeros(1, 2)*dirichlet_functions(x) );
% 
%             array.f0 = b(1:2*n_nodes_u);
% 
%             [A_ns, b_ns] = assembler_steady_navier_stokes( fespace_u, fespace_p, f, 0.0, dirichlet_functions, ...
%                                                            neumann_functions );
% 
%             b = b_ns - A_ns( u_lifting_stokes ) * u_lifting_stokes;
%             b(1:n_nodes_u) = apply_dirichlet_bc_rhs( b(1:n_nodes_u), fespace_u, @(x) zeros(1, 2)*dirichlet_functions(x) );
%             b(n_nodes_u+1:2*n_nodes_u) = apply_dirichlet_bc_rhs( b(n_nodes_u+1:2*n_nodes_u), fespace_u, ... 
%                                                                  @(x) zeros(1, 2)*dirichlet_functions(x) );
%             array.f1 = b(1:2*n_nodes_u);
% 
%             return
% 
%         end
% 
%     end
end


