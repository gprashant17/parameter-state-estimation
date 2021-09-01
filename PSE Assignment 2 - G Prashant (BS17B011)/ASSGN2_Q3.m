%% Loading data
load a2_q3.mat

%% Exploratory Analysis of the data
N = length(vk);
figure();
%Plotting vk
plot(vk,"-")
xlabel("k")
ylabel("v[k]")
title("Realization of Non Stationary Random Process v")

%Calculating ACF and PACF - commands adopted from MATLAB live script
acf_vk = autocorr(vk,'NumLags',20);
pacf_vk = parcorr(vk,'NumLags',20);
% Plot sample ACF and sample PACF
figure
subplot(211)
autocorr(vk,'NumLags',20);
ylabel('Sample ACF');
xlabel('')
set(gca,'fontsize',12,'fontweight','bold');
hca = gca;
set(hca.Children(4),'LineWidth',2,'Color','m')
hca.YLim(2) = 1.1;
box off

subplot(212)
parcorr(vk,'NumLags',20);
ylabel('Sample PACF');
set(gca,'fontsize',12,'fontweight','bold');
hca = gca;
set(hca.Children(4),'LineWidth',2)
box off

%% Systematic determination of ARIMA model parameters

%%Performing ADF test to confirm the presence of integrating effect
[ht_adf,pval] = adftest(vk)
if pval>0.05
    fprintf("Presence of Random Walk effect\n")
end

%% Difference series analysis
vdk = diff(vk);
% Plot the differenced series
figure
plot((2:N),vdk,'linewidth',2)
xlabel('Time Instance k')
ylabel('v[k]-v[k-1]')
set(gca,'fontsize',12,'fontweight','bold')
title('Differenced series','fontsize',13)
box off

%% SAMPLE ACF AND PACF FOR DIFFERENCE SERIES vdk
figure
subplot(211)
autocorr(vdk,'NumLags',20);
ylabel('Sample ACF');
xlabel('')
title('Sample ACF and PACF of differenced series')
set(gca,'fontsize',12,'fontweight','bold');
hca = gca;
set(hca.Children(4),'LineWidth',2,'Color','m')
hca.YLim(2) = 1.1;
box off
subplot(212)
parcorr(vdk,'NumLags',20);
ylabel('Sample PACF');
set(gca,'fontsize',12,'fontweight','bold');
hca = gca;
set(hca.Children(4),'LineWidth',2)
box off

%% Performing ADF test to confirm the presence of integrating effect
[ht_adf,pval] = adftest(vdk);
pval
if pval>0.05
    fprintf("Presence of Random Walk/Integrating effect in the difference series\n")
else
    fprintf("No presence of Random Walk/Integrating effect in the difference series\n")
end

%% Determining P and M of ARMA model with the difference series

%%Build an AR(2,1,1) model
mod_arima211 = arima(2,1,1);
mod_arima211.Constant = 0;
mod_arima211est = estimate(mod_arima211,vk);
result = summarize(mod_arima211est);

%AIC for ARIMA(2,1,1)
result.AIC

res_arima211 = infer(mod_arima211est,vk);

% ACF of residuals plot
figure
autocorr(res_arima211,'NumLags',20)
ylabel('Sample ACF');
set(gca,'fontsize',12,'fontweight','bold');
title('ACF of residuals from AIRMA(2,1,1) model')
hca = gca;
set(hca.Children(4),'LineWidth',2)
box off

%Ljung Box Test
[ht_resarma11,pval] = lbqtest(res_arima211)

%%The codes for performing the analyses were adopted frorm the live script