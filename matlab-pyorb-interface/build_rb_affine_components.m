function [array] = build_rb_affine_components( operator, fem_specifics )
% Assemble fom affine matrix for elliptic scalar problems
% input=
%           operator: operator corresponding to stiffness matrix or RHS 
%           fem_specifics: struct containing the information to build the
%           mesh, the fespace and the chosen model
% output=
%           array: struct containing the affine stiffness matrices in COO format

    fom_problem = load_fom_problem( fem_specifics );
    
    array = fom_problem.build_rb_affine_components( operator, fem_specifics );
end


