classdef LogStructural < FunctionalConnectomeModel
    methods
        function obj = LogStructural()
            % `f_ij = a_ij
            %       + b_ij * log(1 + s_ij)`
            obj.Name = "log linear in structural";
            obj.NumEquationParams = 2;
        end

        function [X, y, missing] = BuildXy(~, functional, structural, ~)
            % Build design matrix X, and y.
            N = size(functional, 1);
            rest = size(functional, 2:3);
            X = cat(2, ...
                ones(                         [N 1 rest]), ...
                reshape(log(1 + structural),  [N 1 rest]));
            y = reshape(functional,           [N 1 rest]);
            % "for each edge where there is data"
            missing = ~any(structural, 1);
            missing = squeeze(missing);
        end
    end
end
