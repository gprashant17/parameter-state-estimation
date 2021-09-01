%% Finding credible interval of theta using Bayesian estimate
rng(42)
N = 10000; %Number of observations

alpha = 0.05; %for 95% probability

theta_0 = 10; %true value of theta - mean of Poisson distbn
yk = poissrnd(theta_0,1,N);  %Sampling N data points from Poisson
theta_b = mean(yk)+1/(2*N)  %Bayesian Estimate (MMSE)of theta

V = 2*N*mean(yk)+1;  %dof for chi^2 distbn
a = chi2inv(alpha/2,V)/(2*N)  %lower bound
b = chi2inv(1-alpha/2,V)/(2*N)  %upper bound 