%% Sampling 100 points and estimating sample covariance
rng(20)
mu_x = 1;  %Mean = 1
std_x = 2;   %Standard deviation = 2
N = 100;
X_sample = mu_x + randn(N,1).*std_x;
Y_sample = 3*X_sample.^2 + 5*X_sample;

%% Estimating Sample Covariance Matrix with code
sigma_hat_code = sample_cov_matrix(X_sample,Y_sample);
fprintf("The estimated Covariance of Y and X with 100 samples is \n")
sigma_hat_code

%% Estimating Sample Covariance Matrix with cov function in MATLAB
sigma_hat_cov = cov(X_sample,Y_sample);
fprintf("The estimated Covariance of Y and X with 100 samples, calculated using cov function is \n")
sigma_hat_cov

%% Simulating with increasing sample sizes
rng(21)
sample_sizes = 100:100:1e6;
cov_vals = zeros(length(sample_sizes),1);

i = 1;
for N = sample_sizes
    X_sample = mu_x + randn(N,1).*std_x;
    Y_sample = 3*X_sample.^2 + 5*X_sample;
    sigma_y_x_hat = sum((X_sample-mean(X_sample)).*(Y_sample-mean(Y_sample)))/N;
    cov_vals(i) = sigma_y_x_hat;
    i = i+1;
end
%Theoretical value of Covariance of X and Y = 44
theoretical_cov = 44;
%% Plotting
figure;
plot(sample_sizes,cov_vals)
line([min(sample_sizes),max(sample_sizes)],[theoretical_cov,theoretical_cov],'Color','red')
xlabel("Sample Size");
legend(["Estimated Covariance","Theoretical Covariance"])
ylabel("Estimate of Covariance of X and Y");
title("Estimated Sample Covariance of X and Y with increasing Sample Size")
%% Function for determining Covariance Matrix of two Samples
function cov_mat = sample_cov_matrix(X,Y)
    N = length(X);
    cov_mat = zeros(2,2);
    %Variance of X
    cov_mat(1,1) = sum((X-mean(X)).*(X-mean(X)))/N;
    %Covariance of X and Y
    cov_mat(1,2) = sum((X-mean(X)).*(Y-mean(Y)))/N;
    cov_mat(2,1) = sum((Y-mean(Y)).*(X-mean(X)))/N;
    %Variance of Y
    cov_mat(2,2) = sum((Y-mean(Y)).*(Y-mean(Y)))/N;
end
