%% Question 4 - Determining optimal MAE predictor of chi square distribution

% xf(x)
func1 = @(x) x.*chi2pdf(x,10);
% E(|X-Xhat|)
func2 = @(x_hat) 2*x_hat.*chi2cdf(x_hat,10)-x_hat - integral(func1,0,x_hat)+integral(func1,x_hat,inf);
% Numerical Optimization to find MAE with starting value 10
[x_star,fval] = fminsearch(func2,10);

fprintf("Optimal MAE predictor value X* = %f\n",x_star)
fprintf("The avaerage absolute error at the optimal value = %f\n",fval)