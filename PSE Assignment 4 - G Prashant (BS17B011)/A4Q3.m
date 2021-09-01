load stack_loss.mat

%% Part a - Determing Correlations between regressors and dependent variable
corr_psi1 = corr(psi1vec,yvec)
corr_psi2 = corr(psi2vec,yvec)
corr_psi3 = corr(psi3vec,yvec)

%% Part b and c - Fitting Linear Model

phi_mat = [psi1vec,psi2vec,psi3vec];
lin_reg = fitlm(phi_mat,yvec)6+

%% Part d - Confidence Interval
N = length(yvec);
alpha = 0.05;
t_alpha = tinv(alpha/2,N-4);
psi1 = 80;
psi2 = 25;
psi3 = 90;
psi_vec = [1 psi1 psi2 psi3];
phi_mat = [ones(N,1) psi1vec psi2vec psi3vec];
y_hat = psi_vec*lin_reg.Coefficients.Estimate;
se = sqrt(sum((yvec-lin_reg.predict).^2)/(N-4));  %Standard Error

lb = y_hat + t_alpha*se*(sqrt(psi_vec*inv(phi_mat'*phi_mat)*(psi_vec')))
ub = y_hat - t_alpha*se*(sqrt(psi_vec*inv(phi_mat'*phi_mat)*(psi_vec')))


%% Part e - Prediction interval

lb = y_hat + t_alpha*se*(sqrt(psi_vec*inv(phi_mat'*phi_mat)*(psi_vec')+1))
ub = y_hat - t_alpha*se*(sqrt(psi_vec*inv(phi_mat'*phi_mat)*(psi_vec')+1))