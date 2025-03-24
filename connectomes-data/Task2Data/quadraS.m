% f = a + b s +y s^2
function [alpha, beta, gamma] = quadraS(rsf_list, wfa_list)
alpha = zeros(68, 68);
beta = zeros(68, 68);
gamma = zeros(68,68);

for i = 1:68
    for j = 1:68
        % if s_ij is 0
        if wfa_list(:,i,j) ==0
            continue 
        end
       
        % Extract x (T) and y (F) values for this pair
        x = squeeze(wfa_list(:, i, j));
        y = squeeze(rsf_list(:,i, j));
        % ✅ Fix: Ensure at least two data points
        coeffs = polyfit(x, y, 2); 

        gamma(i,j) = coeffs(1);
        beta(i, j) = coeffs(2);   % Slope (β)
        alpha( i, j) = coeffs(3);  % Intercept (α)
    end
end