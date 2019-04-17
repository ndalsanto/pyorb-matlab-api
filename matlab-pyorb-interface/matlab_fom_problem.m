classdef matlab_fom_problem
   properties

   end
   methods
      function obj = set_matlab_fem_simulation( obj, ~ )
         disp('Please call set_matlab_fem_simulation in a child of this class');
      end
      
      function [] = solve_parameter( ~, ~ )
         disp('Please call solve_parameter in a child of this class');
      end
      
      function [] = build_fom_affine_components( ~, ~ )
         disp('Please call build_fom_affine_components in a child of this class');
      end
   
      function [] = assemble_fom_matrix( ~, ~ )
         disp('Please call assemble_fom_matrix in a child of this class');
      end
      
      function [] = assemble_fom_rhs( ~, ~ )
         disp('Please call assemble_fom_rhs in a child of this class');
      end
      
      function [] = find_mdeim_elements_fom_specifics( ~, ~ )
         disp('Please call find_mdeim_elements_fom_specifics in a child of this class');
      end
      
      function [] = find_elements_for_deim_fom_specifics( ~, ~ )
         disp('Please call find_elements_for_deim_fom_specifics in a child of this class');
      end
end
end