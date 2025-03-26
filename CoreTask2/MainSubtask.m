function MainSubtask(models, functional, structural, indirect, opts)
%MAINSUBTASK    Answer to questions 2/3 and 4 of core task 2. Consists of:
%               - k-fold cross-validation on AIC, BIC, SSE,
%               - generating plots of variations of parameter estimates and
%                 fit quality across different edges,
%               - plots the relationship described by the N best fitting
%                 models for top– and bottom-k edges
    arguments
        models cell
        functional (:, 68, 68) double
        structural (:, 68, 68) double
        indirect   (:, 68, 68) double
        opts.bestModelsN {mustBeInteger, mustBePositive}
        opts.sortModelsByMeanOf char ...
            {mustBeMember(opts.sortModelsByMeanOf,{'aic','bic','sse'})}
        opts.numEdges {mustBeInteger, mustBePositive}
        opts.cvMethod char ...
            {mustBeMember(opts.cvMethod, {'Kfold', 'LeaveOneOut'})}
        opts.numFolds {mustBeInteger, mustBePositive}
        opts.fitMethod char ...
            {mustBeMember(opts.fitMethod, {'Fit', 'FitPooled'})}
        opts.randomSeed = 0 % default random seed
    end

    opts.bestModelsN = min(length(models), opts.bestModelsN);
    if strcmp(opts.cvMethod, 'LeaveOneOut')
        opts.numFolds = size(functional, 1);
    elseif ~isfield(opts, 'numFolds')
        error("specify numFolds for k-fold cross-validation");
    end

    % Q2/3 (a) Cross-validation on AIC, BIC, SSE
    AICs = zeros(length(models), opts.numFolds);
    BICs = zeros(size(AICs));
    SSEs = zeros(size(AICs));
    RESs = cell(length(models), 1);
    CoefficientVariations = cell(size(RESs));
    TestFolds = cell(size(RESs));
    for m = 1:numel(models)
        model = models{m};
        [aic, bic, sse, res, coeffs, testFold] = CrossValidation( ...
            model, functional, structural, indirect, ...
            "cvMethod", opts.cvMethod, ...
            "numFolds", opts.numFolds, ...
            "fitMethod", opts.fitMethod, ...
            "randomSeed", opts.randomSeed);
        fprintf("-- %s\n", model.Name);
        fprintf("   AIC: %9.2f ± %9.2f | BIC: %9.2f ± %9.2f | SSE: %9.2f ± %9.2f\n", ...
            mean(aic), std(aic), mean(bic), std(bic), mean(sse), std(sse));
        AICs(m, :) = aic;
        BICs(m, :) = bic;
        SSEs(m, :) = sse;
        RESs{m} = squeeze(mean(cell2mat(res), "includenan")); % RES test folds average
        CoefficientVariations{m} = squeeze(std(coeffs, 0, 1));
        TestFolds{m} = testFold;
    end

    % Q2/3 (b)  Consider the variation in parameter estimates and fit
    %           quality across different edges.
    %           NOTE: plots are only useful when fitting each edge
    %           individually.
    if strcmp(opts.fitMethod, 'Fit')
        % Variation in parameter estimates.
        ParameterNames = {"α", "β", "γ", "δ", "ε", "ζ", "η", "θ", "ι", ...
            "κ", "λ", "μ", "ν", "ξ", "ο", "π", "ρ", "σ", "τ", "υ", "φ", ...
            "χ", "ψ", "ω"};
        for m = 1:numel(models)
            coeffVar = CoefficientVariations{m};
            NumCoefficients = size(models{m}.CoefficientEstimates, 1);
            figure;
            tiledlayout(1, NumCoefficients);
            for c = 1:NumCoefficients
                nexttile
                cStd = squeeze(coeffVar(c, :, :));
                imagesc(cStd, "AlphaData", ~isnan(cStd));
                colorbar;
                title(sprintf("std(%s)", ParameterNames{c}));
                axis equal square
                xlim([1 68])
                ylim([1 68])
            end
            sgtitle(models{m}.Name);
            set(gcf, 'Name', "Standard Deviations of Coefficients");
            set(gcf, 'Position', [100, 100, NumCoefficients * 300 + 100, 400]);
        end
        % Fit quality across different edges.
        nrows = 4;
        figure;
        tiledlayout(nrows, ceil(numel(models) / nrows));
        for m = 1:numel(models)
            nexttile
            coeffVar = RESs{m} .^ 2;
            imagesc(coeffVar, "AlphaData", ~isnan(coeffVar));
            colorbar;
            title(models{m}.Name);
            axis equal square
            xlim([1 68])
            ylim([1 68])
        end
        set(gcf, 'Name', "SSE Test Fold Average");
    end

    % Q2/3 (c) Plot the relationship described by the best-fitting models for
    %          a small number of individual connections where it holds
    %          particularly strongly or weakly.
    % Sort models by mean criterion over test folds, look at best N.
    switch opts.sortModelsByMeanOf
        case "aic"
            criterion = AICs;
        case "bic"
            criterion = BICs;
        case "sse"
            criterion = SSEs;
        otherwise
            error("b r u h");
    end
    [~, argsort] = sort(mean(criterion, 2));
    bestModels = models(argsort);
    bestModels = bestModels(1:opts.bestModelsN);
    bestTestFolds = TestFolds(argsort);
    bestTestFolds = bestTestFolds(1:opts.bestModelsN);
    % Plot relationships for best and worst edges.
    for rank = 1:opts.bestModelsN
        figure;
        tiledlayout(1, 2);
        % Average test fold predictions and true values.
        allTestFolds = cat(2, bestTestFolds{rank}{:});
        averageTestFolds = squeeze(mean(allTestFolds, 2));    
        f_true_avg = squeeze(averageTestFolds(1,:,:));
        f_pred_avg = squeeze(averageTestFolds(2,:,:));
        % Sort edges by their average squared error, discarding indices of
        % missing edges.
        squared_error_avg = (f_true_avg - f_pred_avg) .^ 2;
        [~, indices] = sort( ...
            squared_error_avg(:), "ascend", "MissingPlacement", "last");
        missingEdges = sum(isnan(squared_error_avg), "all");
        indices = indices(1:end-missingEdges);
        % Obtain best and worst `numEdges` edges' indices.
        bestEdgesIndices = indices(1:opts.numEdges);
        worstEdgesIndices = indices(end-opts.numEdges+1:end);

        % Plot predictions for best edges.
        f_true = squeeze(f_true_avg(bestEdgesIndices));
        f_pred = squeeze(f_pred_avg(bestEdgesIndices));
        mse = mean((f_true - f_pred) .^ 2);
        nexttile
        plot(f_true, f_pred, 'go', "linewidth", 2);
        hold on
        yEqualsX = refline(1, 0);
        yEqualsX.Color = 'blue';
        yEqualsX.LineStyle = '--';
        xlabel('True Value');
        ylabel('Predicted Value');
        title(sprintf("%d Best-Fitting Edges (MSE = %3.2e)", ...
            opts.numEdges, mse));
        axis equal square

        % Plot predictions for worst edges.
        f_true = squeeze(f_true_avg(worstEdgesIndices));
        f_pred = squeeze(f_pred_avg(worstEdgesIndices));
        mse = mean((f_true - f_pred) .^ 2);
        nexttile
        plot(f_true, f_pred, 'rx', "linewidth", 2);
        hold on
        yEqualsX = refline(1, 0);
        yEqualsX.Color = 'blue';
        yEqualsX.LineStyle = '--';
        xlabel('True Value');
        ylabel('Predicted Value');
        title(sprintf("%d Worst-Fitting Edges (MSE = %3.2e)", ...
            opts.numEdges, mse));
        axis equal square

        sgtitle(opts.fitMethod);
        set(gcf, 'Position', [100 + 300 * rank, 800 - 300 * rank, 800, 400]);
        set(gcf, 'Name', ...
            sprintf("[%s] P%d) %s", ...
            opts.sortModelsByMeanOf, rank, bestModels{rank}.Name));
    end
end