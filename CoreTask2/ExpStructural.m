classdef ExpStructural < FunctionalConnectomeModel
    properties
        Alpha (1,1) double
    end

    methods
        function obj = ExpStructural(alpha)
            % `f_ij = a_ij
            %       + b_ij * exp(-α * s_ij)`
            obj.Name = sprintf("exp(-%.2f * structural)", alpha);
            obj.NumEquationParams = 2;
            obj.Alpha = alpha;
        end

        function [X, y, missing] = BuildXy(obj, functional, structural, ~)
            % Build design matrix X, and y.
            N = size(functional, 1);
            rest = size(functional, 2:3);
            X = cat(2, ...
                ones(                                  [N 1 rest]), ...
                reshape(exp(-obj.Alpha * structural),  [N 1 rest]));
            y = reshape(functional,                    [N 1 rest]);
            % "for each edge where there is data"
            missing = ~any(structural, 1);
            missing = squeeze(missing);
        end
    end
end
