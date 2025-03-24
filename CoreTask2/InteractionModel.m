classdef InteractionModel < FunctionalConnectomeModel
    methods
        function obj = InteractionModel()
            % `f_ij = a_ij
            %       + b_ij * s_ij
            %       + g_ij * t_ij
            %       + d_ij * (s_ij * t_ij)`
            obj.Name = "interaction model";
            obj.NumEquationParams = 4;
        end

        function [X, y, missing] = BuildXy(~, functional, structural, indirect)
            % Build design matrix X, and y.
            N = size(functional, 1);
            rest = size(functional, 2:3);
            interaction = structural .* indirect;
            X = cat(2, ...
                ones(                [N 1 rest]), ...
                reshape(structural,  [N 1 rest]), ...
                reshape(indirect,    [N 1 rest]), ...
                reshape(interaction, [N 1 rest]));
            y = reshape(functional,  [N 1 rest]);
            % "for each edge where there is data"
            missing = ~any(structural, 1);
            missing = squeeze(missing);
        end
    end
end
