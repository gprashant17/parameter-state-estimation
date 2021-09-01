%% Load Data
load NMRlogWell.mat

%% Initialize Hyperparameters and Call BOCD Function
%Initializing hyperparameters of prior distribution
mu0 = 1.15; k = 0.01; alph = 20; bet = 2; lam_CP = 250;
% Call BOCD Function
[R,max_R,pmean] = bocd(mu0,k,alph,bet,lam_CP,y);
%Rotate R matrix and remove the top portion for better vi
R = rot90(R);
R = R(700:end,:);

%% Plotting
% Observations, Predicted means and Changepoints
T = length(y);
figure;
x1 = subplot(2,1,1);
plot(1:T,y);
xlabel('Time instant k')
ylabel('y[k]')
title('Plot of Data y[k]')
xlim([1 1100])
hold on
plot(1:T,pmean,'k-','LineWidth',1.8)
ylim([10 15])
hold on
for xl = [85,220, 452, 812, 873, 941,999]
xline(xl,'g--','LineWidth',1)
end
legend(['Observation'],['Predicted Mean'],['Detected Changepoint'],'Location','north')

% Maximum run length posterior probability plot
x2 = subplot(2,1,2);
plot(1:T,max_R(2:end),'r-')
xlabel('Time instant k')
ylabel('run length')
title('Maximum probability run length')
xlim([1 1100])

% Intensity plot of run length posterior probability distribution
figure;
imshow(imcomplement(flipud(R)*255));
axis('on', 'image')
set(gca,'YDir','normal')
xlabel('Time instant k')
ylabel('run length')
title('Posterior predictive probability map of run length')
xlim([1 1100])

%% Function - BOCD Algorithm
function [R,max_R,pmean] = bocd(mu0,k,alph,bet,lam_CP,y)
%{
Function that implements Bayesian Online Changepoint Detection Algorithm -
Adam and MacKay

Arguments:
mu0 - prior mean hyperparamter
k - Scaling factor prior
alph - initial prior alpha of gamma distribution
bet - initial prior bet of gamma distribution

Outputs:
R - Run Length Posterior Probability matrix
max_R - Run Length Maximum Posterior Probability vector
pmean - predicted mean vector 
%}

% Initializing 
T = length(y);
R = zeros(T+1,T+1);
R(1,1) = 1;
max_R    = zeros(1,T+1);
max_R(1) = 0;
mu_vals = [mu0];
k_vals = [k];
alph_vals = [alph];
bet_vals = [bet];
H = 1/lam_CP;
pmean = zeros(1,T);

% For each incoming data point
for t = 1:T
    x = y(t);
    pi_r = [];
    log_R = log(R);
    % Calculating predictive mean
    pmean(t) = sum(exp(log_R(t, 1:t)) .* mu_vals(1:t));
    % Calculating UPM Probability
    for i = 1:t
        LAM = alph_vals(i)*k_vals(i)/(bet_vals(i)*(k_vals(i)+1));
        stat = (x-mu_vals(i))*sqrt(LAM);
        pi_r = [pi_r sqrt(LAM)*tpdf(stat,2*alph_vals(i))];
    end
    
    % Growth probability
    R(t+1, 2:t+1) = R(t, 1:t) .* pi_r * (1-H);
    % Changepoint Probability
    R(t+1, 1) = sum(R(t, 1:t) .* pi_r * H);
    % Normalization
    R(t+1,:) = R(t+1,:)/sum(R(t+1,:));
    
    % Storing Run length corresponding to maximum probability
    max_R(t+1) = find(R(t+1,:) == max(R(t+1,:)))-1;
    
    % Updating Hyperparameters
    bet_vals = [bet  bet_vals + (k_vals.*(x-mu_vals).^2)./(2*(k_vals+1))];
    mu_vals = [mu0 (mu_vals.*k_vals+x)./(k_vals+1)];
    k_vals = [k k_vals+1];
    alph_vals = [alph alph_vals+1/2];
    
end
log_R = log(R);
pmean(T+1) = sum(exp(log_R(T+1, 1:T+1)) .* mu_vals(1:T+1));
R = R(2:T+1,:);
pmean = pmean(2:end);
end