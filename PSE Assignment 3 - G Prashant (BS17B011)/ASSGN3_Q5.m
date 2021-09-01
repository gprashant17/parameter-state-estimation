true_theta = 5;
rng(42)
N = 1000;
R = 5000;

%Random samples from exponential distribution after shifting mean by theta
y_mat = exprnd(1,R,N)+true_theta;

%MVUE Estimator
Tn = min(y_mat')-1/N;
var_Tn = var(Tn);  %%Variance Lower Bound

%Theta hat 1 = sample mean -1
theta_hat1 = mean(y_mat,2)-1;
var_theta_hat_1 = var(theta_hat1);
eff_theta_hat_1 = var_Tn/var_theta_hat_1   %Efficiency of theta hat 1

%Theta hat 2 = sample median -ln(2)
theta_hat2 = median(y_mat,2)-log(2);
var_theta_hat_2 = var(theta_hat2);
eff_theta_hat_2 = var_Tn/var_theta_hat_2   %Efficiency of theta hat 2

if eff_theta_hat_1>eff_theta_hat_2
    fprintf("\nTheta Hat 1 = Sample Mean - 1 is more efficient\n");
else
    fprintf("\nTheta Hat 2 = Sample Median - ln(2) is more efficient\n")
end
