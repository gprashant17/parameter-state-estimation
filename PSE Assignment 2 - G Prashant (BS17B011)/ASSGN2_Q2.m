%% Generate R realizations of a Random Walk Process process with GWN, simulated for N timepoints

% Sample sizes and realizations
N = 500; R = 500;
%Time points
time_points = 0:(N-1);
%Realization X Observation matrix - ystar and y
y_star_mat = zeros(R,N);
%y matrix
y_mat = zeros(R,N);

%Variance od u
var_u = 2;

%Parameter constants b and f
b = 2;
f = 0.5;

%Variance of e such that signal to noise ratio is 10
var_e = 16/15;

u_mat = zeros(R,1);

for i = 1:(N-1)
    %GWN u
    u_mat(:,i) = sqrt(var_u)*randn(R,1);
    %Calculating ystar(k+1)
    y_star_mat(:,(i+1)) = u_mat(:,i)*b - f*y_star_mat(:,i);
    %Calculating y(k+1)
    y_mat(:,(i+1)) = y_star_mat(:,(i+1)) + sqrt(var_e)*randn(R,1);
end

theor_var_y = b^2 * var_u/(1-f^2) + var_e^2;
fprintf("Theoretical Variance of y = %f\n",theor_var_y)
%Calculating variance of y[300]
empirical_var_y_300 = var(y_mat(:,300));
fprintf("Empirical Variance of y for k=300 = %f\n",empirical_var_y_300)

theor_acvf_y_1 = -f*b^2 * var_u/(1-f^2);
fprintf("Theoretical ACVF of y with lag 1 = %f\n",theor_acvf_y_1)
empirical_acvf_y_1 = cov(y_mat(:,300),y_mat(:,299));
empirical_acvf_y_1 = empirical_acvf_y_1(1,2);
fprintf("Empirical ACVF of y with lag 1 = %f\n",empirical_acvf_y_1)

theor_ccvf_yu_1 = 0;
fprintf("Theoretical CCVF of y and u with lag 1 = %f\n",theor_ccvf_yu_1)

% as index 300th column of u matrix is equivalent to k=299 as per the
% generating fungtion
empirical_ccvf_yu_1 = cov(y_mat(:,300),u_mat(:,300));
empirical_ccvf_yu_1 = empirical_ccvf_yu_1(1,2);

fprintf("Empirical CCVF of y and u with lag 1 = %f\n",empirical_ccvf_yu_1)

theor_ccvf_yu_2 = b*var_u;
fprintf("Theoretical CCVF of y and u with lag 2 = %f\n",theor_ccvf_yu_2)

% as index 299th column of u matrix is equivalent to k=298 as per the
% generating fungtion
empirical_ccvf_yu_2 = cov(y_mat(:,300),u_mat(:,299));
empirical_ccvf_yu_2 = empirical_ccvf_yu_2(1,2);
fprintf("Empirical CCVF of y and u with lag 2 = %f\n",empirical_ccvf_yu_2)

fprintf("\nIt can be observed that the theoretical values and the empirical values were closely equal but did not match exactly due to less number of realizations\n")
