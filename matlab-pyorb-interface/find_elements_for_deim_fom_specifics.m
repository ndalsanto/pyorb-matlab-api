function elements = find_elements_for_deim_fom_specifics( fem_specifics, indices )
% Wrapper for extracting the list of mesh elements contributing to
% populating the given indices of the RHS, used in DEIM
% framework
% input=
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the chosen model
% output=
%           elements: array of elements

    disp('Finding DEIM elements from MATLAB pyorb interface')
    fom_problem = initialize_fom_simulation( fem_specifics );
    elements = fom_problem.find_elements_for_deim_fom_specifics( indices );

%     elements = find_elements_given_vector_indices( fespace, indices );

%     if strcmp( 'navier_stokes', considered_model )
% 
%         [fespace_u, ~] = set_ns_simulation( fem_specifics );
%         n_nodes_u = size(fespace_u.nodes,1);
%         indices_for_scalar_case = indices;
%         indices_for_scalar_case(:, 1) = mod( indices_for_scalar_case(:, 1), n_nodes_u );
% 
%         elements = find_elements_given_vector_indices( fespace_u, indices_for_scalar_case );
% 
%     end
end


