%% Generate R realizations of a Random Walk Process process with GWN, simulated for N timepoints

% Sample sizes and realizations
N = 500; R = 300;
%Time points
time_points = 0:(N-1);
%Realization X Observation matrix
data_mat = zeros(R,N);
%Storing variances across all realizations for different time points
var_vector = zeros(1,N);

for i = 2:N
    %Simulating Random Walk Process with GWN 
    data_mat(:,i) = data_mat(:,i-1)+ randn(R,1);
    %Calculating Variance
    var_vector(i) = var(data_mat(:,i));
end

%Plotting Variance
plot(time_points,var_vector,"-")
xlabel("Time Instance k")
ylabel("Variance Across Realizations")
title("Variance across 300 Realizations at different time points of a Random Walk Process")