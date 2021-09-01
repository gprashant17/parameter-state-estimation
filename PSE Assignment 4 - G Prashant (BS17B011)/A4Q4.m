%% Loading data
load engine_thrust.mat

%% Fitting Linear Regression Model on all regressors including intercept
lin_model = fitlm(Phi,y)

%% Plotting Residual vs order
figure;
plot(1:40',lin_model.Residuals.Raw,'-o')
xlabel("k")
ylabel("Residual e[k]")
yline(0);
title("Plot of residues for all observations")
%% Histogram of residuals
histogram(lin_model.Residuals.Raw,'Normalization','pdf','NumBins',10)
xlabel("Residuals");
ylabel("Normalized Frequency")
title("Histogram of Residuals")

%% 
Phi1 = Phi(:,1);
Phi2 = Phi(:,2);
Phi3 = Phi(:,3);
Phi4 = Phi(:,4);
Phi5 = Phi(:,5);
Phi6 = Phi(:,6);

% Plotting residual plot - change for different regressor
figure;
plot(Phi2,lin_model.Residuals.Raw,'o')
corr1 = corr(lin_model.Residuals.Raw,Phi2);
title(["Residual vs Regressor psi2"]+["; Corr(e,psi2) = "+string(corr1)])
xlabel("psi2[k]")
ylabel("Residual[k]")
hold on 
yline(0);
legend(["Residual","x=0"])
%% Part (b) 
% We notice that insignificant coefs from previous part are c0,c2,c3,c4
lin_model_2 = fitlm(Phi(:,[1,5,6]),y,'Intercept',false)

%% Part (c) - Stepwise Linear Regression
step_lin = stepwiselm(Phi,y,'PEnter',0.1,'PRemove',0.15)

%% Part (d)
% [sorted_res, idx] = sort(lin_model_2.Residuals.Raw)
Phi1 = Phi(:,1);
Phi5 = Phi(:,5);
Phi6 = Phi(:,6);

[sorted_reg1, idx] = sort(Phi1);
[sorted_reg5, idx] = sort(Phi5);
[sorted_reg6, idx] = sort(Phi6);


% sorted_regr = Phi1(idx);
figure;
plot(sorted_reg1,lin_model_2.Residuals.Raw(idx),'o')
corr1 = corr(lin_model_2.Residuals.Raw,Phi1);
title(["Residual vs Regressor psi1"]+["; Corr(e,psi1) = "+string(corr1)])
xlabel("psi1[k]")
ylabel("Residual[k]")
hold on 
yline(0);
legend(["Residual","x=0"])

figure;
plot(sorted_reg5,lin_model_2.Residuals.Raw(idx),'o')
corr5 = corr(lin_model_2.Residuals.Raw,Phi5);

title(["Residual vs Regressor psi5"]+["; Corr(e,psi5) = "+string(corr5)])
xlabel("psi5[k]")
ylabel("Residual[k]")
hold on 
yline(0);

figure;
plot(sorted_reg6,lin_model_2.Residuals.Raw(idx),'o')
corr6 = corr(lin_model_2.Residuals.Raw,Phi6);
title(["Residual vs Regressor psi6"]+["; Corr(e,psi6) = "+string(corr6)])
xlabel("psi6[k]")
ylabel("Residual[k]")
hold on 
yline(0);