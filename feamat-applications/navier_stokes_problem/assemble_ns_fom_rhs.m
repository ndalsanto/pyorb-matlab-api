function [array] = assemble_ns_fom_rhs( navier_stokes_fom_problem, param, varargin )

lifting_name = strcat( navier_stokes_fom_problem.fem_specifics.full_path, navier_stokes_fom_problem.fem_specifics.simulation_name );

u_lifting_stokes = build_stokes_lifting( navier_stokes_fom_problem.fespace_u, navier_stokes_fom_problem.fespace_p, lifting_name );
dirichlet_functions = @(x) [0 0; 0 0; 0 0; 0 0; 0 0; 0 0]';
n_nodes_u = size(navier_stokes_fom_problem.fespace_u.nodes,1);

if nargin == 2
    
    M_velocities = assemble_ns_non_affine_matrix( navier_stokes_fom_problem, param );
    M = sparse(M_velocities.A(:, 1), M_velocities.A(:, 2), M_velocities.A(:, 3), 2*n_nodes_u, 2*n_nodes_u);

    b = - M * u_lifting_stokes(1:2*n_nodes_u);
    b(1:n_nodes_u) = apply_dirichlet_bc_rhs( b(1:n_nodes_u), navier_stokes_fom_problem.fespace_u, @(x) zeros(1, 2)*dirichlet_functions(x) );
    b(n_nodes_u+1:2*n_nodes_u) = apply_dirichlet_bc_rhs( b(n_nodes_u+1:2*n_nodes_u), navier_stokes_fom_problem.fespace_u, @(x) zeros(1, 2)*dirichlet_functions(x) );
    array.f = b(1:2*n_nodes_u);
end

if nargin == 4
%     array.f = assemble_ns_fom_rhs_element_list( param, navier_stokes_fom_problem.fespace_u, varargin{1}, u_lifting_stokes(1:2*n_nodes_u) );
%     array.f(1:n_nodes_u) = apply_dirichlet_bc_rhs( array.f(1:n_nodes_u), navier_stokes_fom_problem.fespace_u, @(x) zeros(1, 2)*dirichlet_functions(x) );
%     array.f(n_nodes_u+1:2*n_nodes_u) = apply_dirichlet_bc_rhs( array.f(n_nodes_u+1:2*n_nodes_u), navier_stokes_fom_problem.fespace_u, @(x) zeros(1, 2)*dirichlet_functions(x) );
%     array.f = array.f(varargin{2});

    M_velocities = assemble_ns_non_affine_matrix( navier_stokes_fom_problem, param, varargin{1} );
    M = sparse(M_velocities.A(:, 1), M_velocities.A(:, 2), M_velocities.A(:, 3), 2*n_nodes_u, 2*n_nodes_u);
    array.f = - M(varargin{2}, :) * u_lifting_stokes(1:2*n_nodes_u);
end

end