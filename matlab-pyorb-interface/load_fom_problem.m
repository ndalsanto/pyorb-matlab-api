function [matlab_fom_problem_instance] = load_fom_problem( fom_specifics )
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
    
%     matlab_fom_problem_instance = initialize_fom_simulation( fom_specifics );
    
end