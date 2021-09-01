%% Load Historic and New Data
load historic.mat
yh = y';
load new.mat
yn = y';
clear y;

%% ACF/PACF plots
% Code adopted from MATLAB Live Script
acf_yh = autocorr(yh,'NumLags',20);
pacf_yh = parcorr(yh,'NumLags',20);
% Plot sample ACF and sample PACF
figure
subplot(211)
autocorr(yh,'NumLags',20);
ylabel('Sample ACF');
xlabel('')
set(gca,'fontsize',12,'fontweight','bold');
hca = gca;
set(hca.Children(4),'LineWidth',2,'Color','m')
hca.YLim(2) = 1.1;
box off
subplot(212)
parcorr(yh,'NumLags',20);
ylabel('Sample PACF');
set(gca,'fontsize',12,'fontweight','bold');
hca = gca;
set(hca.Children(4),'LineWidth',2)
box off

%% Fitting an AR Model of order 1
% Code adopted from MATLAB Live Script
mod_ar1 = arima(1,0,0);
mod_ar1.Constant = 0;
mod_ar1est = estimate(mod_ar1,yh);
[res_ar1,~,logL] = infer(mod_ar1est,yh);
% ACF of residuals
figure;
autocorr(res_ar1,'NumLags',20)
ylabel('Sample ACF');
set(gca,'fontsize',12,'fontweight','bold');
title('ACF of residuals from AR(1) model')
hca = gca;
set(hca.Children(4),'LineWidth',2)
box off
% Whiteness test of residuals
[ht_resar3,pval] = lbqtest(res_ar1)

%% Using Normal Equation to determine AR(1) coefficient (Optional)
%Uncomment Below Lines
% reg = [0 yh(1:end-1)];
% out = yh;
% theta0 = reg\out;
%% RLS + BOCD
N = length(yn);
% Initial AR(1) Coefficient
theta0 = mod_ar1est.AR{1};
% Initial inverse of Covariance
P0 = 1/(yh'*yh);
% RLS Object
rls_obj = recursiveLS (1,'InitialParameters',theta0 ,'InitialParameterCovariance',P0);

T = N - 200;

%Initializing hyperparameters of prior distribution
mu0 = 0; k = 1; alph = 10; bet = 1; lam_CP = 50; H = 1/lam_CP;
% Initializing 
R = zeros(T+1,T+1);
R(1,1) = 1;
max_R    = zeros(1,T+1);
max_R(1) = 1;
mu_vals = [mu0];
k_vals = [k];
alph_vals = [alph];
bet_vals = [bet];
flag = 0;
pmean = zeros(1,T);
theta_vals = zeros(1,N);
ypred = [];
innovat = [];
for i = 1:N
    if i == 1
        x = yh(end);
    else
        x = yn(i-1);
    end
    %RLS Coef and output estimate
    [theta,yhat] = rls_obj(yn(i),x);
    theta_vals(i) = theta;
    ypred(i) = yhat;
    %Calculating Innovation
    e = yn(i) - yhat;
    innovat(i) = e;
    if i > 200 && flag == 0
        xt = e;
        pi_r = [];
        t = i-200;
        log_R = log(R);
        % Calculating predictive mean
        pmean(t) = sum(exp(log_R(t, 1:t)) .* mu_vals(1:t)); 
        % Calculating UPM Probability
        for j = 1:t
            LAM = alph_vals(j)*k_vals(j)/(bet_vals(j)*(k_vals(j)+1));
            stat = (x-mu_vals(j))*sqrt(LAM);
            pi_r = [pi_r sqrt(LAM)*tpdf(stat,2*alph_vals(j))];
        end
        % Growth probability
        R(t+1, 2:t+1) = R(t, 1:t) .* pi_r * (1-H);
        % Changepoint Probability
        R(t+1, 1) = sum(R(t, 1:t) .* pi_r * H);
        % Normalization
        R(t+1,:) = R(t+1,:)/sum(R(t+1,:));
        
        % Storing Run length corresponding to maximum probability
        max_R(t+1) = find(R(t+1,:) == max(R(t+1,:)))-1;
        % Terminate BOCD if CP is detected
        if max_R(t+1) - max_R(t) < max_R(t+1) - 5 && max_R(t+1) < 5
            CP_Time = t+200;
            fprintf('Changepoint Found at %d\n',CP_Time)
            flag = 1;
        end
        % Updating Hyperparameters
        bet_vals = [bet  bet_vals + (k_vals.*(x-mu_vals).^2)./(2*(k_vals+1))];
        mu_vals = [mu0 (mu_vals.*k_vals+x)./(k_vals+1)];
        k_vals = [k k_vals+1];
        alph_vals = [alph alph_vals+1/2];
    end
    % Store Coeff for 1st part
    if flag == 1
        model1_theta = rls_obj.Parameters
        model1_cov = rls_obj.ParameterCovariance;
        % Initiate another RLS
        rls_obj = recursiveLS (1 ,'InitialParameters' , rls_obj.Parameters ,'InitialParameterCovariance',rls_obj.ParameterCovariance);
        flag = 2;
    end
end

log_R = log(R);
pmean(t+1) = sum(exp(log_R(t+1, 1:t+1)) .* mu_vals(1:t+1));
R = R(2:T+1,:);
pmean = pmean(2:end);

% Store Coeff for 2nd part
model2_theta = rls_obj.Parameters
model2_cov = rls_obj.ParameterCovariance;

%% Plot
%Streaming Data
figure;
subplot(2,1,1)
plot(1:N,yn)
hold on
% plot(1:N,ypred,'k-')
xline(CP_Time,'r--','LineWidth',1.5)
xlabel('Time instant k')
ylabel('y_n[k]')
title('Plot of streaming data y_n[k]')
legend(['Incoming data'],['Detected Changepoint'],'Location','northwest')

%Innovation
subplot(2,1,2)
plot(1:N,innovat)
xline(CP_Time,'r--','LineWidth',1.5)
xlabel('Time instant k')
ylabel('e[k]')
title('Plot of Innovations e[k] = y_n[k] - yhat[k]')
legend(['Innovation (residuals)'],['Detected Changepoint'],'Location','southwest')

%Estimated Theta, with time
figure;
plot(1:N,theta_vals)
hold on
xline(CP_Time,'r--','LineWidth',1.5)
xlabel('Time instant k')
ylabel('theta[k]')
title('Plot of AR(1) coefficient against time')
legend(['Parameter estimate'],['Detected Changepoint'],'Location','northwest')
