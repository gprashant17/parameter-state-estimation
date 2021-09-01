%% Simulation
theta_0 = 3;  %True value -> 3
med_0 = theta_0/sqrt(2);
cdf_ = @(y) y.^2/theta_0.^2;   %CDF of the Density function

%Generating samples from the distribution
full = linspace(0,theta_0,1e6);
full_prob = cdf_(full);
sum(full_prob)
N_vals = 10:10:1000;
med_est = zeros(length(N_vals),1);
z=1;
for N = N_vals
    r = rand(N,1);
    yk = zeros(N,1);
    j=1;
    for i=1:N
        ind = find(r(i)<full_prob,1);
        yk(j) = full(ind(1));
        j = j+1;
    end
    med_est(z) = (1/sqrt(2))*prod(yk)^(1/N)*((2*N+1)/(2*N))^N; %Estimating Median
    z = z+1;
end

%% Plotting estimates with increasing N
figure;
plot(N_vals,med_est,"-")
hold on
plot(N_vals,med_0*ones(1,length(N_vals)),'-r')
xlabel("N")
ylabel("Estimate of Median")
title("Consistency of ML estimate of Median")
legend(["ML Median Estimate","True Value"])

%% Testing asymptotic Gaussian of median estimate
R = 100; %Realizations
N = 10000;  %Large number of Observations
med_r = zeros(1,R);
for z = 1:R
    r = rand(N,1);
    yk = zeros(N,1);
    parfor (i=1:N,4)
        ind = find(r(i)<full_prob,1);
        yk(i) = full(ind(1));
    end
    med_r(z) = (1/sqrt(2))*prod(yk.^(1/N))*((1+(1/(2*N)))^N);
end
%% Histogram of estimates for N = 10000
hist(med_r,20)
xlabel("Median Estimate")
ylabel("Frequency")
title("Histogram of Bayesian Estimate of Meidan - N=10000")
