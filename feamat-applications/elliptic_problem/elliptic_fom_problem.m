classdef elliptic_fom_problem < matlab_fom_problem

    properties
        fespace
        mesh
        fem_specifics
    end
    
    methods
      function [obj] = set_matlab_fem_simulation( obj, fem_specifics )

        disp( 'Initializing FEM elliptic problem' )
        obj.fem_specifics = fem_specifics;
        
        n_elements_x = fem_specifics.number_of_elements;
        poly_degree  = fem_specifics.polynomial_degree;
        n_elements_y = n_elements_x;

        bottom_left_corner_x = 0;
        bottom_left_corner_y = 0;

        L = 1.0;
        H = 1.0;

        mesh_name = fem_specifics.mesh_name;

        obj.mesh = create_mesh( bottom_left_corner_x, bottom_left_corner_y, ...
                                L, H, n_elements_x, n_elements_y, mesh_name );

        current_model = fem_specifics.model;
        current_dirichlet = fem_specifics.use_nonhomogeneous_dirichlet;

        bc_flags = define_bc_flags( current_model, current_dirichlet );

        obj.fespace = create_fespace( obj.mesh, poly_degree, bc_flags );

      end
      
      function [sol] = solve_parameter( obj, param )
        
        current_model = obj.fem_specifics.model;

        mu = obj.build_diffusion( param, current_model );
        [f, dirichlet_functions, neumann_functions] = obj.build_source_and_bc( param );

        if strcmp(current_model, 'nonaffine') && strcmp( obj.fem_specifics.use_nonhomogeneous_dirichlet, 'Y' )
            non_hom_dirichlet_functions = @(x) [1;0;0;0];
            [ A, b ] = assembler_poisson( obj.fespace, f, mu, non_hom_dirichlet_functions, neumann_functions );
            uL = b * 0.0;
            uL = apply_dirichlet_bc_rhs( uL, obj.fespace, non_hom_dirichlet_functions );
            b = b - A * uL;
            b = apply_dirichlet_bc_rhs( b, obj.fespace, dirichlet_functions );
            [ A, ~ ] = assembler_poisson( obj.fespace, f, mu, dirichlet_functions, neumann_functions );
        else
            [ A, b ] = assembler_poisson( obj.fespace, f, mu, dirichlet_functions, neumann_functions );
        end

        sol.u  = A \ b;

      end


      function [array] = assemble_fom_natural_norm_matrix( obj )

        f = @(x) 0.*x(1,:);
        dirichlet_functions = @(x) [0;0;0;0];
        neumann_functions   = @(x) [0;0;0;0];

        [ A, b ] = assembler_poisson( obj.fespace, f, 1, dirichlet_functions, neumann_functions );
        [i,j,val] = find( A );
        array.A = [i,j,val];

      end

      function [f, dirichlet_functions, neumann_functions] = build_source_and_bc( obj, param )

        current_model = obj.fem_specifics.model;
          
        if strcmp( current_model, 'thermal_block' )        
            dirichlet_functions = @(x) [0;0;0;0];
            neumann_functions = @(x) [1;0;0;0];
            f = @(x) 0*x(1,:);
        end

        if strcmp( current_model, 'nonaffine' )        
            f = @(x) 0*x(1,:) + 1;
            dirichlet_functions = @(x) [0;0;0;0];
            neumann_functions   = @(x) [0;0;0;0];
        end

        if strcmp( current_model, 'elliptic_example' )
            f = @(x) 0.*x(1,:) + param(2) * ( (x(1,:)-0.5).^2 + (x(2,:)-0.5).^2 < 0.01 );
            dirichlet_functions = @(x) [0;0;0;0];
            neumann_functions   = @(x) [0;0;param(3);0];
        end


    end

      
      function mu = build_diffusion( obj, param, current_model )
        if strcmp( current_model, 'elliptic_example' )
            mu = @(x) ( (x(1,:)-0.5).^2 + (x(2,:)-0.5).^2 < 0.01 ) ...
                      + param(1) * ( (x(1,:)-0.5).^2 + (x(2,:)-0.5).^2 >= 0.01 );
            return 
        end

        if strcmp( current_model, 'thermal_block' ) || strcmp( current_model, 'nonaffine_thermal_block')
            mu = @(x) param(1)*(x(1,:)<0.5).*(x(2,:)<0.5) ...
            + param(2)*(x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.0).*(x(2,:)<0.5) ...
            + param(3)*(x(1,:)>=0.0).*(x(1,:)<0.5).*(x(2,:)>=0.5).*(x(2,:)<1.0) ...
            + 1.0 * (x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.5).*(x(2,:)<1.0) ;
            return 
        end

        if strcmp( current_model, 'nonaffine' )
            mu = @(x) 1. * ( param(3) + ( 1. / param(3) ) ... 
               * exp( - ( ( x(1,:)-param(1) ) .* ( x(1,:)-param(1) ) ...
                        + ( x(2,:)-param(2) ) .* ( x(2,:)-param(2) ) ) / param(3) ) );

            return 
        end
      end
      
      function [array] = build_fom_affine_components( obj, operator, ~ )

        considered_model = obj.fem_specifics.model;

        if ( strcmp( operator, 'f' ) == 0 ) && ( strcmp( considered_model, 'thermal_block' ) == 0 && strcmp( considered_model, 'elliptic_example' ) == 0 )
            error('This operator for the chosen model is not supported');
        end

        trivial_parameter = zeros( 10 );
        [f, dirichlet_functions, neumann_functions] = obj.build_source_and_bc( trivial_parameter, considered_model );
        
        if strcmp( 'thermal_block', considered_model )

            if operator == 'A'

                mu = @(x) (x(1,:)<0.5).*(x(2,:)<0.5);
                [ A, ~ ] = assembler_poisson( obj.fespace,f,mu,dirichlet_functions,neumann_functions );

                [i,j,val] = find( A );
                array.A0 = [i,j,val];

                mu = @(x) (x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.0).*(x(2,:)<0.5);
                [ A, ~ ] = assembler_poisson( obj.fespace,f,mu,dirichlet_functions,neumann_functions );

                [i,j,val] = find( A );
                array.A1 = [i,j,val];

                mu = @(x) (x(1,:)>=0.0).*(x(1,:)<0.5).*(x(2,:)>=0.5).*(x(2,:)<1.0);
                [ A, ~ ] = assembler_poisson( obj.fespace,f,mu,dirichlet_functions,neumann_functions );

                [i,j,val] = find( A );
                array.A2 = [i,j,val];

                mu = @(x) 1.0 * (x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.5).*(x(2,:)<1.0) ;
                [ A, ~ ] = assembler_poisson( obj.fespace,f,mu,dirichlet_functions,neumann_functions );

                [i,j,val] = find( A );
                array.A3 = [i,j,val];

            end

            if operator == 'f'
                mu = @(x) 1.0 * (x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.5).*(x(2,:)<1.0) ;
                [ ~, b ] = assembler_poisson( obj.fespace,f,mu,dirichlet_functions,neumann_functions );
                array.f0 = b;

            end
        end
        
        if strcmp( 'elliptic_example', considered_model )
            
            if operator == 'A'
                mu = @(x) ( (x(1,:)-0.5).^2 + (x(2,:)-0.5).^2 < 0.01 );
                [ A, ~ ] = assembler_poisson( obj.fespace,f,mu,dirichlet_functions,neumann_functions );
                [i,j,val] = find( A );
                array.A0 = [i,j,val];

                mu = @(x) ( (x(1,:)-0.5).^2 + (x(2,:)-0.5).^2 >= 0.01 );
                [ A, ~ ] = assembler_poisson( obj.fespace,f,mu,dirichlet_functions,neumann_functions );
                [i,j,val] = find( A );
                array.A1 = [i,j,val];
            end
            
            if operator == 'f'
                mu = 0.0;
                f = @(x) 0.*x(1,:) + ( (x(1,:)-0.5).^2 + (x(2,:)-0.5).^2 < 0.01 );
                homogeneous_neumann_functions = @(x) [0;0;0;0];
                [ ~, b ] = assembler_poisson( obj.fespace, f, mu, dirichlet_functions, homogeneous_neumann_functions );
                array.f0 = b;
                
                ones_neumann_functions   = @(x) [0;0;1;0];
                f = @(x) 0.*x(1,:);
                [ ~, b ] = assembler_poisson( obj.fespace, f, mu, dirichlet_functions, ones_neumann_functions );
                array.f1 = b;
                
            end
        end
            
      end
      
      function [array] = assemble_fom_matrix( obj, param, varargin )

        dirichlet_functions = @(x) [0;0;0;0];
        neumann_functions = @(x) [1;0;0;0];

        % forcing term (not employed)
        f = @(x) 0*x(1,:);

        current_model = obj.fem_specifics.model;

        mu = obj.build_diffusion( param, current_model );

        if nargin < 4
            [ A, ~ ] = assembler_poisson( obj.fespace, f, mu,dirichlet_functions, neumann_functions );
        else
            element_list = varargin{1};
            [ A, ~ ] = assembler_poisson( obj.fespace, f, mu,dirichlet_functions, neumann_functions, element_list );
        end

        if nargin < 4
           [ i, j, val ] = find( A );
           array.A = [ i, j, val ];
        end
        if nargin >= 4
            indeces_list = varargin{2}; % supposed to be a matrix of size nb_of_indices x 2 
            elements_A = zeros( size( indeces_list, 1 ), 1 );
            for ii = 1:size( indeces_list , 1 )
                elements_A(ii) = A(indeces_list(ii, 1), indeces_list(ii, 2));
            end

            array.A = [indeces_list, elements_A];
        end
      end
      
    function [array] = assemble_fom_rhs( obj, param, varargin )
        if nargin < 4
            array = assemble_elliptic_fom_rhs( obj, param );
        end
        if nargin >= 4
            array = assemble_elliptic_fom_rhs( obj, param, varargin{1}, varargin{2} );
        end
    end
    
    function [elements] = find_mdeim_elements_fom_specifics( obj, indices_mat )
        elements= find_elements_given_matrix_indices( obj.fespace, obj.fespace, indices_mat );
    end
    
    function elements = find_elements_for_deim_fom_specifics( obj, indices )
       elements = find_elements_given_vector_indices( obj.fespace, indices ); 
    end
        
  end
end


