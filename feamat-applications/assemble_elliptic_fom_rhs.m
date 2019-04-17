function [array] = assemble_elliptic_fom_rhs( elliptic_problem, param, varargin )

    % forcing term
    f = @(x) 0*x(1,:);

    current_model = elliptic_problem.fem_specifics.model;
    current_dirichlet = elliptic_problem.fem_specifics.use_nonhomogeneous_dirichlet;

    mu = elliptic_problem.build_diffusion( param, current_model );
    
    if strcmp( current_model, 'nonaffine' )

        f = @(x) 0*x(1,:) + 1;
        
        dirichlet_functions = @(x) [0;0;0;0];
        neumann_functions   = @(x) [0;0;0;0];
        
    end

    if strcmp( current_dirichlet, 'Y' )
        non_hom_dirichlet_functions = @(x) [1;0;0;0];
        if nargin == 2
            [ A, b ] = assembler_poisson( elliptic_problem.fespace, f, mu, dirichlet_functions, neumann_functions );
            uL = b * 0.0;
            uL = apply_dirichlet_bc_rhs( uL, elliptic_problem.fespace, non_hom_dirichlet_functions );
            b = b - A * uL;
            b = apply_dirichlet_bc_rhs( b, elliptic_problem.fespace, dirichlet_functions );
        else
            uL = zeros( size(elliptic_problem.fespace.nodes, 1), 1 );
            uL = apply_dirichlet_bc_rhs( uL, elliptic_problem.fespace, non_hom_dirichlet_functions );
            b = assemble_rhs_elementlist( elliptic_problem.fespace, f, mu, varargin{1}, uL );
            b = apply_dirichlet_bc_rhs( b, elliptic_problem.fespace, dirichlet_functions );
        end
    else
        [ ~, b ] = assembler_poisson( elliptic_problem.fespace, f, mu, dirichlet_functions, neumann_functions );
    end

    if nargin > 2
        b = b(varargin{2});
    end

    array.f = b;

end