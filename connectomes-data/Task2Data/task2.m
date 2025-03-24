% task 2 
%  wfa names are structual, rsf are function 
wfa_list = zeros(19,68,68);
rsf_list = zeros(19,68,68);
for i= 1:19
    name = i+31;
    filename_wfa = sprintf('%d_WFA_68.csv', name);
    filename_rsf = sprintf('%d_rsfMRI_68.csv', name);
    wfa_list(i,:,:) = csvread(filename_wfa,1,0);
    rsf_list(i,:,:) = csvread(filename_rsf,1,0);
end

%%
% t_ij = max_k { min { s_ik, s_kj } } s.t. sik, skj ≠ 0.
T = zeros(19, 68, 68); % Preallocate matrix for t_ij

for wfa_index=1:19
    S= squeeze(wfa_list(wfa_index,:,:));
    for i = 1:68
        for j = 1:68
            min_values = zeros(1, 68); % Store min(s_ik, s_kj) for all k    
            for k = 1:68
                if S(i,k)==0 ||S(k, j)==0
                    continue
                end
                min_values(k) = min(S(i, k), S(k, j)); % Compute min { s_ik, s_kj }
            end
            T(wfa_index,i, j) = max(min_values); % Take max over k
        end
    end
end

%% fit the lienar models to the data 
%  f_ij = alpha + beta_ij t_ij
% Initialize matrices to store α and β coefficients
alpha = zeros(68, 68);
beta = zeros(68, 68);


for i = 1:68
    for j = 1:68
        % Extract x (T) and y (F) values for this pair
        x = squeeze(T(:, i, j));
        y = squeeze(rsf_list(:, i, j));
       
        % Fit a linear model: y = alpha + beta * x
        coeffs = polyfit(x,y ,1);
        beta(i, j) = coeffs(1);   % Slope (β)
        alpha(i, j) = coeffs(2);  % Intercept (α)
    end
end


% Display results
disp('Alpha coefficients (Intercepts):');
disp(alpha);
disp('Beta coefficients (Slopes):');
disp(beta);

%% 4 f_ij = aij + beta_ij t_ij + gamma_ij t_ij^2
alpha = zeros(68, 68);
beta = zeros(68, 68);
gamma = zeros(68,68);

for i = 1:68
    for j = 1:68
        % Extract x (T) and y (F) values for this pair
        t = squeeze(T(:, i, j));
        x = [ones(length(t)) t t.^2 ];
        y = squeeze(rsf_list(:,i, j));
        % [coeffs,~,~] = polyfit(x, y, 2); 
        coeffs = x \y;
        gamma(i,j) = coeffs(3);
        beta(i, j) = coeffs(2);   % Slope (β)
        alpha( i, j) = coeffs(1);  % Intercept (α)
    end
end 

%% t2 linear model 1
[alpha, beta]=linearM1(rsf_list,wfa_list);

%% t2 model 2 
[alpha, beta, gamma] = quadraS(rsf_list, wfa_list);
 
%% t2 model 5 
[alpha, beta, gamma] = Model5(rsf_list, wfa_list);

%% task 3 cross validation for 5 models d
% model 1 :x =[1, S_ij]
% modeling the relationship between structural and functional connectivity for each individual edge across subjects.
% Leave-one-out is great when:
% You have very limited data (like 19 subjects).
% You want to evaluate generalization on each individual.
% You care about stability across folds (since every subject gets tested).
x = zeros(19,2,68,68);
y= rsf_list;
for i=1:68
    for j =1:68
        x(:,1,i,j)=1;
        x(:,2,i,j)=wfa_list(:,i,j);
    end
end
% for model 1: fij = αij + βij sij
[mse, aic,bic]=leave1OutCV(x,y,2);

%% model 2 fij = αij + βij sij + γij sij2
x= zeros(19,3,68,68);
y= rsf_list;
for i=1:68
    for j =1:68
        x(:,1,i,j)=1;
        x(:,2,i,j)=wfa_list(:,i,j);
        x(:,3,i,j) = wfa_list(:,i,j).^2;
    end
end
% for model 1: fij = αij + βij sij
[mse, aic,bic]=leave1OutCV(x,y,3);

%% model 3 fij = αij + βij tij







