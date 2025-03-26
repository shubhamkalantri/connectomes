classdef SigmoidStructuralIndirect < FunctionalConnectomeModel
    properties
        Alpha (1,1) double
        Theta (1,1) double
    end

    methods
        function obj = SigmoidStructuralIndirect(alpha, theta)
            % `f_ij = a_ij
            %       + b_ij * [1 + exp(-α * (s_ij - θ))]^(-1)
            %       + g_ij * t_ij`
            obj.Name = sprintf(...
                "sigmoid structural (α=%.2f, θ=%.2f) + indirect", alpha, theta);
            obj.NumEquationParams = 3;
            obj.Alpha = alpha;
            obj.Theta = theta;
        end

        function [X, y, missing] = BuildXy(obj, functional, structural, indirect)
            % Build design matrix X, and y.
            N = size(functional, 1);
            rest = size(functional, 2:3);
            X = cat(2, ...
                ones(                                  [N 1 rest]), ...
                reshape( ...
                1 ./ (exp(-obj.Alpha .* (structural - obj.Theta))), ...
                [N 1 rest]), ...
                reshape(indirect, [N 1 rest]));
            y = reshape(functional,                    [N 1 rest]);
            % "for each edge where there is data"
            missing = ~any(structural, 1);
            missing = squeeze(missing);
        end
    end
end
