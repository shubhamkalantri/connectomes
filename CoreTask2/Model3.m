classdef Model3 < FunctionalConnectomeModel
    methods
        function obj = Model3()
            % `f_ij = a_ij + b_ij * t_ij`
            obj.Name = "linear in indirect";
            obj.NumEquationParams = 2;
        end

        function [X, y, missing] = BuildXy(~, functional, ~, indirect)
            % Build design matrix and y.
            N = size(functional, 1);
            rest = size(functional, 2:3);
            X = cat(2, ...
                ones(               [N 1 rest]), ...
                reshape(indirect,   [N 1 rest]));
            X = reshape(X,           N, 2, []); % NumModelParams == 2
            y = reshape(functional,  N, 1, []);
            missing = zeros(rest);
        end
    end
end
