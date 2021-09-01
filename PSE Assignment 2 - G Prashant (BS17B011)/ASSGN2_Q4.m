%% MLE and Linear Regression - Question 4
rng(1)
%Number of observations
N = 100;
%Randomly generating x vector
xk = 3 + 0.3*randn(N,1);
%True values of parameters a and b (a=2,b=3)
a_true = 2;
b_true = 3;
sig_e = 1;
%Random GWN signal
epsilon = randn(N,sig_e^2);
%Determining true values of y
yk = a_true*xk+b_true+epsilon;

%Defining Log Likelihood function
loglik = @(a,b) sum((yk-a.*xk-b).^2)/(-2*sig_e^2);

%Finding minima point of negative likelihood using fminsearch with 0,0
%initialization
params = fminsearch(@(p) -loglik(p(1),p(2)),[0;0]);
a_mle_emp = params(1);
b_mle_emp = params(2);

%MLE estimates of a and b obtained through analytical approach
a_mle_analytical = cov(xk,yk)/var(xk);
a_mle_analytical = a_mle_analytical(1,2);
b_mle_analytical = mean(yk) - a_mle_analytical*mean(xk);

fprintf("True Values of parameters\n")
fprintf("a = %f\n",a_true)
fprintf("b = %f\n\n",b_true)

fprintf("MLE estimates of parameters obtained through MATLAB simulation\n")
fprintf("a = %f\n",a_mle_emp)
fprintf("b = %f\n\n",b_mle_emp)

fprintf("MLE estimates of parameters obtained through analytical approach\n")
fprintf("a = %f\n",a_mle_analytical)
fprintf("b = %f\n",b_mle_analytical)

%% Plotting 
a_vals = linspace(-10,10,100);
b_vals = linspace(-10,10,100);
[a_vals,b_vals] = meshgrid(a_vals,b_vals);
likelihoods = zeros(size(a_vals));
s = size(likelihoods);
for i = 1:s(1)
    for j=1:s(2)
        likelihoods(i,j) = loglik(a_vals(i,j),b_vals(i,j));
    end
end
figure;
surf(a_vals,b_vals,likelihoods)
xlabel("a")
ylabel("b")
zlabel("log likelihood")
title("Likelihood as a function of parameters a and b")
hold on

%Locating MLE point on the surface with a red dot
h = plot3(a_mle_emp,b_mle_emp,loglik(a_mle_emp,b_mle_emp),"ro");
legend(["Likelihood Surface","MLE Estimate"])
set(h, 'MarkerFaceColor', get(h,'Color'));

%MLE Line fitting
figure;
plot(xk,yk,"o")
xlabel("xk")
ylabel("yk")
x = linspace(min(xk),max(xk),5);
y = a_mle_analytical*x+b_mle_analytical;
hold on
plot(x,y)
legend(["Data","MLE Linear Regression Fit"])
title("Linear Fit using MLE")