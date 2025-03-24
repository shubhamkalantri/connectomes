classdef Model2 < FunctionalConnectomeModel
    methods
        function obj = Model2()
            % `f_ij = a_ij + b_ij * s_ij + g_ij * s_ij^2`
            obj.Name = "quadratic in structural";
            obj.NumEquationParams = 3;
        end

        function [X, y, missing] = BuildXy(~, functional, structural, ~)
            % Build design matrix and y.
            N = size(functional, 1);
            rest = size(functional, 2:3);
            X = cat(2, ...
                ones(                  [N 1 rest]), ...
                reshape(structural,    [N 1 rest]), ...
                reshape(structural.^2, [N 1 rest]));
            y = reshape(functional,    [N 1 rest]);
            % "for each edge for which there is data"
            missing = ~any(structural, 1);
            missing = squeeze(missing);
        end
    end
end
