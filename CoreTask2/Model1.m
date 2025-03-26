classdef Model1 < FunctionalConnectomeModel
    methods
        function obj = Model1()
            % `f_ij = a_ij + b_ij * s_ij`
            obj.Name = "linear in structural";
            obj.NumEquationParams = 2;
        end

        function [X, y, missing] = BuildXy(~, functional, structural, ~)
            % Build design matrix X, and y.
            N = size(functional, 1);
            rest = size(functional, 2:3);
            X = cat(2, ...
                ones(               [N 1 rest]), ...
                reshape(structural, [N 1 rest]));
            y = reshape(functional, [N 1 rest]);
            % "for each edge where there is data"
            missing = ~any(structural, 1);
            missing = squeeze(missing);
        end
    end
end
