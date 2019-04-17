function [elements] = find_mdeim_elements_fom_specifics( fem_specifics, indices_mat )
% Wrapper for extracting the list of mesh elements contributing to
% populating the given indices of the stiffness matrix, used in MDEIM
% framework
% input=
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the chosen model
% output=
%           elements: array of elements

    disp('Finding MDEIM elements from MATLAB pyorb interface')

    fom_problem = initialize_fom_simulation( fem_specifics );
    
    elements = fom_problem.find_mdeim_elements_fom_specifics( indices_mat );
    
%     if strcmp( 'navier_stokes', considered_model )
%         [fespace_u, ~] = set_ns_simulation( fem_specifics );
% 
%         n_nodes_u = size(fespace_u.nodes,1);
% 
%         indices_mat_for_scalar_case = indices_mat;
%         indices_mat_for_scalar_case(:, 1:2) = mod( indices_mat_for_scalar_case(:, 1:2), n_nodes_u );
%         
%         elements = find_elements_given_matrix_indices( fespace_u, fespace_u, indices_mat_for_scalar_case );
%     end

end
