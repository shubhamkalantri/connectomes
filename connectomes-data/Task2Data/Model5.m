% f = a + b s +y t^2
function [alpha, beta, gamma] = Model5(rsf_list, wfa_list)
alpha = zeros(68, 68);
beta = zeros(68, 68);
gamma = zeros(68,68);

for i = 1:68
    for j = 1:68
        S = squeeze(wfa_list(:,i,j));
        T = squeeze(rsf_list(:,i,j));
        % if s_ij is 0
        if wfa_list(:,i,j) ==0
            continue 
        end
        S = (S - mean(S)) / std(S);
        T = (T - mean(T)) / std(T);

       
        % Extract x (T) and y (F) values for this pair
        x = [ones(length(S)) S T];
        y = squeeze(rsf_list(:,i, j));
        coeffs = x \ y;

        gamma(i,j) = coeffs(3);
        beta(i, j) = coeffs(2);   % Slope (β)
        alpha( i, j) = coeffs(1);  % Intercept (α)
    end
end