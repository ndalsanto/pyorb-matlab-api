function array = assemble_fom_natural_norm_matrix( fem_specifics )

    fom_problem = load_fom_problem( fem_specifics );

    array = fom_problem.assemble_fom_natural_norm_matrix( );

end