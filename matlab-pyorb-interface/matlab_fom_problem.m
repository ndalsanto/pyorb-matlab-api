classdef matlab_fom_problem
   properties

   end
   methods
      function obj = set_matlab_fem_simulation( obj, ~ )
         disp('Please call a child of this class');
      end
      
      function [] = solve_parameter( ~ )
         disp('Please call a child of this class');
      end
end
end