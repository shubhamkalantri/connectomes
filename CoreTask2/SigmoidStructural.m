classdef SigmoidStructural < FunctionalConnectomeModel
    properties
        Alpha (1,1) double
        Theta (1,1) double
    end

    methods
        function obj = SigmoidStructural(alpha, theta)
            % `f_ij = a_ij
            %       + b_ij * [1 + exp(-α * (s_ij - θ))]^(-1)`
            obj.Name = sprintf("sigmoid structural (α=%f, θ=%f)", alpha, theta);
            obj.NumEquationParams = 2;
            obj.Alpha = alpha;
            obj.Theta = theta;
        end

        function [X, y, missing] = BuildXy(obj, functional, structural, ~)
            % Build design matrix X, and y.
            N = size(functional, 1);
            rest = size(functional, 2:3);
            X = cat(2, ...
                ones(                                  [N 1 rest]), ...
                reshape( ...
                1 ./ (exp(-obj.Alpha .* (structural - obj.Theta))), ...
                [N 1 rest]));
            y = reshape(functional,                    [N 1 rest]);
            % "for each edge where there is data"
            missing = ~any(structural, 1);
            missing = squeeze(missing);
        end
    end
end
