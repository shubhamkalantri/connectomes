% x is the data, y is the target and k is the number of parameters
function [CVerr, AIC_all, BIC_all]= leave1OutCV(x,y,k)
CVerr = NaN(68, 68);        % Store CV MSE per edge
AIC_all = NaN(68, 68);      % Store average AIC per edge
BIC_all = NaN(68, 68); 

for i = 1:68
    for j = 1:68
        x_edge = squeeze(x(:,:, i, j));
        y_edge = squeeze(y(:, i, j));
        if x_edge == 0
            continue;
        end
        aic_total =0;
        bic_total =0;
        sse = 0;
        for fold = 1:19
            test_idx = fold;
            train_idx = setdiff(1:19, test_idx);
            x_train = x_edge(train_idx,:);
            y_train = y_edge(train_idx);
            x_test = x_edge(test_idx,:);
            y_test = y_edge(test_idx);

            coeffs = x_train \ y_train;
            yhat_test =x_test * coeffs;
            % disp(size(yhat_test))
            % Update test error (squared error)
            sse = sse + (yhat_test - y_test)^2;
            % Compute AIC & BIC on training data
            yhat_train = x_train * coeffs;
            residuals = y_train - yhat_train;
            sigma2 = mean(residuals.^2);
            L = - (length(train_idx)/2) * log(2 * pi * sigma2) - ...
                (1 / (2 * sigma2)) * sum(residuals.^2);
            aic_total = aic_total + (2 * k - 2 * L);
            bic_total = bic_total + (k * log(length(train_idx)) - 2 * L);
        end
        CVerr(i, j) = sse/ 19;
        AIC_all(i, j) = aic_total / 19;
        BIC_all(i, j) = bic_total / 19;
    end
   
end