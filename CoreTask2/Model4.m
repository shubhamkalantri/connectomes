classdef Model4 < FunctionalConnectomeModel
    methods
        function obj = Model4()
            % `f_ij = a_ij + b_ij * t_ij + g_ij * t_ij^2`
            obj.Name = "quadratic in indirect";
            obj.NumEquationParams = 3;
        end

        function [X, y, missing] = BuildXy(~, functional, ~, indirect)
            N = size(functional, 1);
            rest = size(functional, 2:3);
            X = cat(2, ...
                ones(                [N  1  rest]), ...
                reshape(indirect,    [N  1  rest]), ...
                reshape(indirect.^2, [N  1  rest]));
            X = reshape(X,            N, 3, []);   % NumModelParams == 3
            y = reshape(functional,   N, 1, []);
            missing = zeros(rest);
        end
    end
end
