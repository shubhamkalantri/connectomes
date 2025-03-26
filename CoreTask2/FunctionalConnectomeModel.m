classdef (Abstract) FunctionalConnectomeModel < handle
    properties (SetAccess = protected)
        % NAME  Model's friendly name, purely for debugging purposes
        Name (1,1) string
        % NUMEQUATIONPARAMS    Number of parameters in the model's equation
        NumEquationParams (1,1) {mustBeNonnegative, mustBeInteger}
        % COEFFICIENTESTIMATES  Model's coefficients, populated by concrete
        %                       implementations
        CoefficientEstimates (:, 68, 68) double
    end

    methods (Abstract)
        % BUILDXY   Computes the design matrix X and target variable y for
        %           this model, given the functional and structural
        %           connectomes and the indirect connectivity matrix.
        %           For inputs/outputs with shape (N, A, B), we have:
        %               X's shape: (N, NumEquationParams, A, B)
        %               y's shape: (N, 1,                 A, B)
        % 
        %           The values of A and B depend on whether the model is
        %           fitting each individual edge or some aggregation.
        %           When fitting each edge, (A, B) == (68, 68).
        [X, y, missing] = BuildXy(obj, functional, structural, indirect)
    end

    methods (Static)
        function res = Residuals(y_true, y_pred)
            % RESIDUALS Computes the residuals between the true and
            %           predicted value.
            res = y_true - y_pred;
        end

        function sse = SSE(y_true, y_pred)
            % SSE   Computes the sum of squared errors between the true and
            %       predicted value, omitting NaNs from the sum.
            %       NaNs are omitted to make AIC and BIC values useful.
            residuals = FunctionalConnectomeModel.Residuals(y_true, y_pred);
            sse = sum(residuals(:).^2, "omitnan");
        end
    end
    
    methods
        function ResetCoefficientEstimates(obj)
            % RESETCOEFFICIENTESTIMATES Resets coefficient estimates to NaN
            %                           matrix of appropriate dimensions.
            obj.CoefficientEstimates = nan(obj.NumEquationParams, 68, 68);
        end

        function k = NumParams(obj, fitMethod)
            % NUMPARAMS     Returns the number of model parameters for a
            %               given fitting method. Used when calculating AIC
            %               and BIC.
            arguments
                obj FunctionalConnectomeModel
                fitMethod char {mustBeMember(fitMethod, ...
                    {'Fit', 'FitPooled'})}
            end
            if strcmp(fitMethod, 'Fit')
                % 'Fit' fits each edge, so the number of parameters is the
                % number of non-NaN coefficient estimates.
                % NOTE: We're counting the number of non-NaN coefficients
                %       in the lower (wlog) triangular portion of the
                %       coefficient matrix because it is symmetric.
                k = 0;
                for p = 1:obj.NumEquationParams
                    C = squeeze(~isnan(obj.CoefficientEstimates(p,:,:)));
                    k = k + sum(tril(C), "all");
                end
            elseif strcmp(fitMethod, 'FitPooled')
                % 'FitPooled' fits the model to the aggregration of the
                % edge data, treating all edge data as a single observation
                % so the number of parameters is the number of coefficients
                % in the model's equation.
                k = obj.NumEquationParams;
            else
                error("unreachable");
            end
        end

        function aic = AIC(obj, y_true, y_pred, fitMethod)
            % AIC       Computes the Akaike Information Criterion.
            n = numel(y_true);
            sse = obj.SSE(y_true, y_pred);
            aic = n * log(sse / n) + 2 * obj.NumParams(fitMethod);
        end
        
        function bic = BIC(obj, y_true, y_pred, fitMethod)
            % BIC       Computes the Bayesian Information Criterion.
            n = numel(y_true);
            sse = obj.SSE(y_true, y_pred);
            bic = n * log(sse / n) + obj.NumParams(fitMethod) * log(n);
        end

        function Fit(obj, functional, structural, indirect)
            % FIT   Fits the model's linear model relating the functional
            %       connectome and the structural and indirect structural
            %       connectivity connectomes for each individual edge.
            obj.ResetCoefficientEstimates();
            [X, y, missing] = obj.BuildXy(functional, structural, indirect);
            X = X(:, :, ~missing);
            y = y(:, :, ~missing);
            coeffs = pagemldivide(X, y);
            coeffs = squeeze(coeffs);
            obj.CoefficientEstimates(:, ~missing) = coeffs;
        end

        function FitPooled(obj, functional, structural, indirect)
            % FITPOOLED Combines each subject's edge data into a single
            %           observation and fits the model's linear model to
            %           that pooled vector. This results in a single set of
            %           coefficients.
            %           NOTE: since we're treating all available edge data
            %           as a single observation, I don't think it makes
            %           sense to then still make the distinction of which
            %           edges are missing data when updating the model's
            %           coefficients.
            obj.ResetCoefficientEstimates();
            [X, y, missing] = obj.BuildXy(functional, structural, indirect);
            X = X(:, :, ~missing);
            y = y(:, :, ~missing);
            X_pooled = reshape(X, [], obj.NumEquationParams);
            y_pooled = reshape(y, [], 1);
            coeffs = pagemldivide(X_pooled, y_pooled);
            % X_pooled has shape: (NumSubjects * 68 * 68 - NaNs, k)
            % y_pooled has shape: (NumSubjects * 68 * 68 - NaNs, 1)
            % coeffs has shape:   (k, 1)
            for p = 1:obj.NumEquationParams
                obj.CoefficientEstimates(p, :) = coeffs(p);
            end
        end

        function f_pred = Predict(obj, structural, indirect)
            % PREDICT   Outputs the model's prediction given a structural
            %           and indirect structural connectivity connectomes.
            %           NOTE: coefficients that weren't estimated are NaN.
            
            % Obtain and reshape design matrix.
            N = size(structural, 1);
            dummy = zeros(N, 68, 68);
            X = obj.BuildXy(dummy, structural, indirect);
            k = obj.NumEquationParams;
            X = reshape(X, [N k 68 68]);

            % Multiply design matrix by coefficients.
            C = reshape(obj.CoefficientEstimates, [k 1 68 68]);
            f_pred = pagemtimes(X, C);
            f_pred = reshape(f_pred, [N 68 68]);
        end
    end
end
