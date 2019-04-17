function [array] = build_fom_affine_components( operator, fem_specifics )
% Assemble fom affine matrix for elliptic scalar problems
% input=
%           operator: operator corresponding to stiffness matrix or RHS 
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the chosen model
% output=
%           array: struct containing the affine stiffness matrices in COO format
    

    disp('Building affine FOM components from MATLAB pyorb interface')

    fom_problem = initialize_fom_simulation( fem_specifics );
    
    array = fom_problem.build_fom_affine_components( operator, fem_specifics );
end


