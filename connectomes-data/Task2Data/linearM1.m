function [alpha, beta] =  linearM1(rsf_list,wfa_list)
alpha = zeros(68, 68);
beta = zeros(68, 68);
for i = 1:68
    for j = 1:68
        if wfa_list(:,i,j) ==0
            continue 
        end
        x = squeeze(wfa_list(:,i,j));
        y = squeeze(rsf_list(:,i,j));
        % if S ==0
        %     continue
        % end

        
        % Fit a linear model: y = alpha + beta * x
        coeffs = polyfit(x,y,1);
        beta( i, j) = coeffs(1);   % Slope (β)
        alpha( i, j) = coeffs(2);  % Intercept (α)
    end
end

end