function bc_flags = define_bc_flags( current_model, current_dirichlet )

if strcmp( current_model, 'thermal_block' )
    bc_flags = [0 0 1 0];
end

if strcmp( current_model, 'nonaffine' )
    bc_flags = [1 1 1 1];
end

if strcmp( current_model, 'nonaffine' ) && strcmp( current_dirichlet, 'Y' )
    bc_flags = [1 1 0 1];
end

if strcmp( current_model, 'elliptic_example' )
    bc_flags = [1 1 0 1];
end

end
