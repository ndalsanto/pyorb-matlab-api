function [matlab_fom_problem_instance] = initialize_fom_simulation( fom_specifics )
%SET_FEM_SIMULATION Summary of this function goes here
%   Detailed explanation goes here

    fom_problem_name = strcat( fom_specifics.full_path, fom_specifics.simulation_name, '_matlab_fom' );

    file_name_fom_problem = strcat( fom_problem_name, '.mat');
    file_already_exists = exist( file_name_fom_problem, 'file');

    if file_already_exists == 2
        disp( 'I am loading the MATLAB fom problem since it exists from file ' );
        disp(file_name_fom_problem)
        load( file_name_fom_problem, '-mat', 'matlab_fom_problem_instance' );
        return
    end

    if strcmp( fom_specifics.model, 'nonaffine' ) || strcmp( fom_specifics.model, 'thermal_block' )
        disp( 'I am building the MATLAB fom problem since it does not exist from file ' );
        matlab_fom_problem_instance = elliptic_fom_problem;
    end
    
    if strcmp( fom_specifics.model, 'navier_stokes' )
        disp( 'I am building the MATLAB fom problem since it does not exist from file ' );
        matlab_fom_problem_instance = navier_stokes_fom_problem;
    end

    matlab_fom_problem_instance = set_matlab_fem_simulation( matlab_fom_problem_instance, fom_specifics );
    save( file_name_fom_problem, 'matlab_fom_problem_instance' );
end
