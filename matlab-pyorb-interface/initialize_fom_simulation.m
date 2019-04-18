function initialize_fom_simulation( fom_specifics )

    disp('Initializing simulation ')

    fom_problem_name = strcat( fom_specifics.full_path, fom_specifics.simulation_name, '_matlab_fom' );
    file_name_fom_problem = strcat( fom_problem_name, '.mat');

    if strcmp( fom_specifics.model, 'nonaffine' ) || strcmp( fom_specifics.model, 'thermal_block' )
        disp( 'I am building the MATLAB fom problem since it does not exist from file ' );
        matlab_fom_problem_instance = elliptic_fom_problem;
    end
    
    if strcmp( fom_specifics.model, 'navier_stokes' )
        disp( 'I am building the MATLAB fom ns problem since it does not exist from file ' );
        matlab_fom_problem_instance = navier_stokes_fom_problem;
    end

    matlab_fom_problem_instance = set_matlab_fem_simulation( matlab_fom_problem_instance, fom_specifics );
    save( file_name_fom_problem, 'matlab_fom_problem_instance' );
end
