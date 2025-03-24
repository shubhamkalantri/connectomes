function [AIC, BIC, SSE, RES, CoefficientEstimates, TestFolds] = ...
    KFoldCrossValidation(model, functional, structural, indirect, opts)
%KFOLDCROSSVALIDATION   Runs k-fold cross-validation on a
%                       FunctionalConnectomeModel, returning the model's
%                       AIC, BIC, SSE, and RES(iduals) over each test fold,
%                       as well as the test folds (f_true, f_pred).
    arguments
        model FunctionalConnectomeModel
        functional  (:, 68, 68) double
        structural  (:, 68, 68) double
        indirect    (:, 68, 68) double
        opts.cvFolds {mustBeInteger, mustBePositive}
        opts.fitMethod char ...
            {mustBeMember(opts.fitMethod, ...
            {'Fit', 'FitPooled'})}
        opts.randomSeed = 0 % default random seed
    end
    rng(opts.randomSeed);
    k = opts.cvFolds;
    indices = crossvalind("Kfold", size(functional, 1), k);
    AIC = zeros(k, 1);
    BIC = zeros(k, 1);
    SSE = zeros(k, 1);
    RES = cell(k, 1);
    CoefficientEstimates = zeros([k model.NumEquationParams 68 68]);
    TestFolds = cell(k, 1);
    for i = 1:k
        test_idx = (indices == i);
        train_idx = ~test_idx;
        % Train.
        f_train = functional(train_idx, :, :);
        s_train = structural(train_idx, :, :);
        t_train = indirect(train_idx, :, :);
        model.(opts.fitMethod)(f_train, s_train, t_train);
        CoefficientEstimates(k, :, :, :) = model.CoefficientEstimates;
        % Evaluate.
        f_true = functional(test_idx, :, :);
        s_test = structural(test_idx, :, :);
        t_test = indirect(test_idx, :, :);
        f_pred = model.Predict(s_test, t_test);
        TestFolds{i} = cat(1, ...
            reshape(f_true, [1 size(f_true)]), ...
            reshape(f_pred, [1 size(f_pred)]));
        RES{i} = model.Residuals(f_true, f_pred);
        AIC(i) = model.AIC(f_true, f_pred, opts.fitMethod);
        BIC(i) = model.BIC(f_true, f_pred, opts.fitMethod);
        SSE(i) = model.SSE(f_true, f_pred);
    end
end