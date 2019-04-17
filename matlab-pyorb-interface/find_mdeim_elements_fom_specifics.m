function [elements] = find_mdeim_elements_fom_specifics( fem_specifics, indices_mat )
% Wrapper for extracting the list of mesh elements contributing to
% populating the given indices of the stiffness matrix, used in MDEIM
% framework
% input=
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the chosen model
% output=
%           elements: array of elements

    fom_problem = initialize_fom_simulation( fem_specifics );    
    elements = fom_problem.find_mdeim_elements_fom_specifics( indices_mat );

end
