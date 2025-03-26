function [AIC, BIC, SSE, RES, CoefficientEstimates, TestFolds] = ...
    CrossValidation(model, functional, structural, indirect, opts)
%CROSSVALIDATION    Runs k-fold/leave-one-out cross-validation on a
%                   FunctionalConnectomeModel, returning the model's AIC,
%                   BIC, SSE, and RES(iduals) over each test fold, as well
%                   as the test folds (f_true, f_pred).
    arguments
        model FunctionalConnectomeModel
        functional  (:, 68, 68) double
        structural  (:, 68, 68) double
        indirect    (:, 68, 68) double
        opts.cvMethod char ...
            {mustBeMember(opts.cvMethod, {'Kfold', 'LeaveOneOut'})}
        opts.numFolds {mustBeInteger, mustBePositive}
        opts.fitMethod char ...
            {mustBeMember(opts.fitMethod, {'Fit', 'FitPooled'})}
        opts.randomSeed = 0 % default random seed
    end
    rng(opts.randomSeed);
    N = size(functional, 1);
    switch opts.cvMethod
        case 'Kfold'
            numFolds = opts.numFolds;
            indices = crossvalind("Kfold", N, numFolds);
        case 'LeaveOneOut'
            numFolds = N;
            indices = 1:numFolds;
        otherwise
            error("br uh");
    end
    AIC = zeros(numFolds, 1);
    BIC = zeros(numFolds, 1);
    SSE = zeros(numFolds, 1);
    RES = cell(numFolds, 1);
    CoefficientEstimates = zeros([numFolds model.NumEquationParams 68 68]);
    TestFolds = cell(numFolds, 1);
    for i = 1:numFolds
        test_idx = (indices == i);
        train_idx = ~test_idx;
        % Train.
        f_train = functional(train_idx, :, :);
        s_train = structural(train_idx, :, :);
        t_train = indirect(train_idx, :, :);
        model.(opts.fitMethod)(f_train, s_train, t_train);
        CoefficientEstimates(i, :, :, :) = model.CoefficientEstimates;
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