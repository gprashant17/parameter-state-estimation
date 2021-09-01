%% Load Data
load qdist_data.mat

%% Histogram
N = length(Ydata);
figure;
histfit(Ydata,10,"kernel");
xlabel('y')
ylabel('Freq')
title('Histogram fit of given data - 10 bins')

%% Cullen-Frey Graph
figure;
cfgraph(Ydata)
%% AD Test for multiple distributions
[hnorm,pnorm] = adtest(Ydata,'Distribution','normal');
[hexp,pexp] = adtest(Ydata,'Distribution','exp');
[hln,plm] = adtest(Ydata,'Distribution','lognormal');
[hwei,pwei] = adtest(Ydata,'Distribution','weibull');

%% Fit Distribution GUI
%writetable(Ydata,"Q1.txt")
Main_FitDistribution_GUI