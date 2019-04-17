function [u_lifting_stokes] = build_stokes_lifting( fespace_u, fespace_p, lifting_name )

lifting_name = strcat( lifting_name, '_stokes_lifting.mat' );
file_already_exists = exist( lifting_name, 'file');

if file_already_exists == 2
    load( lifting_name, '-mat', 'u_lifting_stokes' );
else
    f = [0; 0];
    r = 0.4;
    r1 = -r;
    r2 =  r;
    U = 10;
    v_in_scalar = @(x)  - U * 4 * (r1 - x) .* (r2 - x) .* (r1 - x < 0) .* (r2 - x > 0);
    v_in = @(x) v_in_scalar(x(2));
    dirichlet_stokes_functions    = @(x) [0 0; 0 0; 0 0; 0 0; 0 0; v_in(x) 0]';
    neumann_functions             = @(x) [0 0;0 0;0 0;0 0;0 0;0 0]';
    [A_stokes, b_stokes] = assembler_steady_stokes( fespace_u, fespace_p, f, 1.0, dirichlet_stokes_functions, ...
                                                    neumann_functions);
    u_lifting_stokes = A_stokes \ b_stokes;

    save( lifting_name, 'u_lifting_stokes' );
end

end
