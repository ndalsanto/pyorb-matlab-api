classdef navier_stokes_fom_problem < matlab_fom_problem
    properties
        fespace_u
        fespace_p
        mesh
        fem_specifics
    end
    
    methods
      function obj = set_matlab_fem_simulation( obj, fem_specifics )
        mesh_file = fem_specifics.mesh_name;
        obj.mesh = read_mesh( mesh_file );
        bc_flags    = [1 0 1 0 1 1];
        obj.fespace_u = create_fespace( obj.mesh, 'P2', bc_flags );
        obj.fespace_p = create_fespace( obj.mesh, 'P1', bc_flags );
        obj.fem_specifics = fem_specifics;
      end

      function [sol] = solve_parameter( obj, param )

        disp('NS solver with parameter')
        disp(param)

        f = [0; 0];
        mu = 0.01 * param(1);

        lifting_name = strcat( obj.fem_specifics.full_path, obj.fem_specifics.simulation_name );

        u_lifting_stokes = build_stokes_lifting( obj.fespace_u, obj.fespace_p, lifting_name );

        dirichlet_functions = @(x) [0 0; 0 0; 0 0; 0 0; 0 0; 0 0]';
        neumann_functions   = @(x) [0 0;0 0;0 0;0 0;0 0;0 0]';

        n_nodes_u = size(obj.fespace_u.nodes,1);
        n_nodes_p = size(obj.fespace_p.nodes,1);

        C_1 = sparse( 2*n_nodes_u+n_nodes_p, 2*n_nodes_u+n_nodes_p );
        C_1 = add_convective_term( C_1, u_lifting_stokes, obj.fespace_u );
        C_1 = add_flipped_convective_term( C_1, u_lifting_stokes, obj.fespace_u );

        B_zeros = sparse( n_nodes_p, 2*n_nodes_u );
        M_velocities = assemble_ns_non_affine_matrix( obj, param );
        M_velocities = M_velocities.A;
        M = [sparse(M_velocities(:, 1), M_velocities(:, 2), M_velocities(:, 3), 2*n_nodes_u, 2*n_nodes_u), (B_zeros'); ...
                B_zeros, sparse( n_nodes_p, n_nodes_p ) ];

        [A_no_lifting, b_no_lifting] = assembler_steady_navier_stokes( obj.fespace_u, obj.fespace_p, ...
                                                                       f, mu, dirichlet_functions, ...
                                                                       neumann_functions );

        % building rhs of the equation coming from lifting
        u_lifting_stokes(1+2*n_nodes_u:end) = 0;
        b = b_no_lifting - ( A_no_lifting(u_lifting_stokes)) * u_lifting_stokes;

        b = b - M * u_lifting_stokes;
        b(1:n_nodes_u) = apply_dirichlet_bc_rhs( b(1:n_nodes_u), obj.fespace_u, @(x) zeros(1, 2)*dirichlet_functions(x) );
        b(n_nodes_u+1:2*n_nodes_u) = apply_dirichlet_bc_rhs( b(n_nodes_u+1:2*n_nodes_u), obj.fespace_u, @(x) zeros(1, 2)*dirichlet_functions(x) );

        A = @(u) A_no_lifting(u) + C_1 + M;

        % solve stokes problem with u = 0 to get initial guess for newton's method
        u0 = zeros(size( obj.fespace_u.nodes, 1) * 2,1 );
        x0 = A(u0) \ b;

        % solve system with newton's method
        method.name = 'newton';
        method.f = @(u) A(u)*u-b;
        method.x0 = x0;
        method.jac = @(u) build_jac_navier_stokes( A, u, obj.fespace_u );
        method.tol = 1e-8;
        method.maxit = 100;

        [ns_sol,~,~] = solve_fluid_system( A, b, obj.fespace_u, obj.fespace_p, method) ;

        sol.u = [ns_sol.u1; ns_sol.u2];
        
      end
      
      function [array] = build_fom_affine_components( obj, operator, fem_specifics )

        n_nodes_u = size( obj.fespace_u.nodes, 1 );
        n_nodes_p = size( obj.fespace_p.nodes, 1 );
        dirichlet_stokes_functions    = @(x) [0 0;0 0;0 0;0 0;0 0;0 0]';
        neumann_functions             = @(x) [0 0;0 0;0 0;0 0;0 0;0 0]';
        f = [0; 0];

        if operator == 'A'

            disp( strcat( 'Retrieving operator A affine components for Navier Stokes - ', num2str(n_nodes_u) ) );
            
            lifting_name = strcat( obj.fem_specifics.full_path, obj.fem_specifics.simulation_name );
            u_lifting_stokes = build_stokes_lifting( obj.fespace_u, obj.fespace_p, lifting_name );
            u_lifting_stokes(1+2*n_nodes_u:end) = 0;

            C_1 = sparse( 2*n_nodes_u+n_nodes_p, 2*n_nodes_u+n_nodes_p );
            C_1 = add_convective_term( C_1, u_lifting_stokes, obj.fespace_u );
            C_1 = add_flipped_convective_term( C_1, u_lifting_stokes, obj.fespace_u );

            [i,j,val] = find( C_1(1:2*n_nodes_u, 1:2*n_nodes_u) );
            array.A0 = [i,j,val];

            [A_stokes, ~] = assembler_steady_stokes( obj.fespace_u, obj.fespace_p, f, 1., dirichlet_stokes_functions, ...
                                                     neumann_functions);

            [i,j,val] = find( A_stokes(1:2*n_nodes_u, 1:2*n_nodes_u) );
            array.A1 = [i,j,val];

            for iB = 1:obj.fem_specifics.range_rb_functions
                
                disp( strcat('Considering RB function - ', num2str(iB) ) );

                pyorb_u = [fem_specifics.("rb_func_" + num2str(iB-1))'; zeros(n_nodes_p, 1)];
                C_1 = sparse(2*n_nodes_u+n_nodes_p, 2*n_nodes_u+n_nodes_p);
                C_1 = add_convective_term( C_1, pyorb_u, obj.fespace_u );

                [i,j,val] = find( C_1(1:2*n_nodes_u, 1:2*n_nodes_u) );
                array.("A" + num2str(iB+1) ) = [i,j,val];
            end
            
            % needed when building the Jacobian
            for iB = 1:obj.fem_specifics.range_rb_functions
                
                disp( strcat('Considering RB function with flipped term - ', num2str(iB) ) );

                pyorb_u = [fem_specifics.("rb_func_" + num2str(iB-1))'; zeros(n_nodes_p, 1)];
                C_1 = sparse(2*n_nodes_u+n_nodes_p, 2*n_nodes_u+n_nodes_p);
                C_1 = add_flipped_convective_term( C_1, pyorb_u, obj.fespace_u );

                [i,j,val] = find( C_1(1:2*n_nodes_u, 1:2*n_nodes_u) );
                array.("A" + num2str( obj.fem_specifics.range_rb_functions+iB+1) ) = [i,j,val];
            end

            disp( 'Finished to build fom NS affine components ' );

            return

        end

        if operator == 'f'
            lifting_name = strcat( obj.fem_specifics.full_path, obj.fem_specifics.simulation_name );
            u_lifting_stokes = build_stokes_lifting( obj.fespace_u, obj.fespace_p, lifting_name );
            u_lifting_stokes(1+2*n_nodes_u:end) = 0;
            dirichlet_functions = @(x) [0 0;0 0;0 0;0 0;0 0;0 0]';


            [A_stokes, b_stokes] = assembler_steady_stokes( obj.fespace_u, obj.fespace_p, f, 1.0, dirichlet_functions, ...
                                                            neumann_functions );

            b = b_stokes - A_stokes * u_lifting_stokes;

            b(1:n_nodes_u) = apply_dirichlet_bc_rhs( b(1:n_nodes_u), obj.fespace_u, @(x) zeros(1, 2)*dirichlet_functions(x) );
            b(n_nodes_u+1:2*n_nodes_u) = apply_dirichlet_bc_rhs( b(n_nodes_u+1:2*n_nodes_u), obj.fespace_u, ... 
                                                                 @(x) zeros(1, 2)*dirichlet_functions(x) );

            array.f0 = b(1:2*n_nodes_u);

            [A_ns, b_ns] = assembler_steady_navier_stokes( obj.fespace_u, obj.fespace_p, f, 0.0, dirichlet_functions, ...
                                                           neumann_functions );

            b = b_ns - A_ns( u_lifting_stokes ) * u_lifting_stokes;
            b(1:n_nodes_u) = apply_dirichlet_bc_rhs( b(1:n_nodes_u), obj.fespace_u, @(x) zeros(1, 2)*dirichlet_functions(x) );
            b(n_nodes_u+1:2*n_nodes_u) = apply_dirichlet_bc_rhs( b(n_nodes_u+1:2*n_nodes_u), obj.fespace_u, ... 
                                                                 @(x) zeros(1, 2)*dirichlet_functions(x) );
            array.f1 = b(1:2*n_nodes_u);

            return

        end

      end
   
      function [array] = assemble_fom_matrix( obj, param, varargin )
         
        if nargin < 4
          array = assemble_ns_non_affine_matrix( obj, param );
        end
        if nargin >= 4
          array = assemble_ns_non_affine_matrix( obj, param, varargin{1}, varargin{2} );
        end
        
      end
      
      function [array] = assemble_fom_rhs( obj, param, varargin )
                  
        if nargin < 4
          array = assemble_ns_fom_rhs( obj, param );
        end
        if nargin >= 4
          array = assemble_ns_fom_rhs( obj, param, varargin{1}, varargin{2} );
        end 
          
      end
      
      function [elements] = find_mdeim_elements_fom_specifics( obj, indices_mat )
        n_nodes_u = size(obj.fespace_u.nodes,1);

        indices_mat_for_scalar_case = indices_mat;
        indices_mat_for_scalar_case(:, 1:2) = mod( indices_mat_for_scalar_case(:, 1:2), n_nodes_u );
        
        elements = find_elements_given_matrix_indices( obj.fespace_u, obj.fespace_u, indices_mat_for_scalar_case );
      end
      
      function elements = find_elements_for_deim_fom_specifics( obj, indices )
        n_nodes_u = size( obj.fespace_u.nodes, 1 );
        indices_for_scalar_case = indices;
        indices_for_scalar_case(:, 1) = mod( indices_for_scalar_case(:, 1), n_nodes_u );

        elements = find_elements_given_vector_indices( obj.fespace_u, indices_for_scalar_case );
      end
    end
end