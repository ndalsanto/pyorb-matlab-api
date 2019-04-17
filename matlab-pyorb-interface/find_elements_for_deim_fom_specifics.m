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

end


