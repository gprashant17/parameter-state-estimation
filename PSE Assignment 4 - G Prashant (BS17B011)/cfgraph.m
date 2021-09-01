function stats_sk = cfgraph(data,sksq_lim,kurt_lim)
% FUNCTION FOR GENERATING CULLEN-FREY GRAPH 
% Works for Univariate data only
%
% Requires the Statistics and ML Toolbox
%
% For CH5115: Parameter and State Estimation
%
% Inputs:
%           data: Observations (univariate)
% Optional:
%           sksq_lim: [Lower Upper] limits of Sq. Skewness for plotting
%           kurt_lim: [Lower Upper] limits of Sq. Skewness for plotting
% Outputs:
%           stats_sk: Observed and bootstrapped skewness and kurtosis
%                         (structure object)
%
% Arun K. Tangirala
% November 25, 2020

if (nargin == 1)
    setxy_lim = 0; 
elseif (nargin == 2)
    error('Supply limits for both squared skewness and kurtosis');
elseif (nargin == 3)
    setxy_lim = 1;
end
    
% Sample size
data = data(:);
N = length(data);

% Compute skewness and kurtosis
skew_stat = skewness(data);
kurt_stat = kurtosis(data,0); % Unbiased estimate

% Bootstrapped statistics
B = min(N,50);
skew_bsp = bootstrp(B,@skewness,data);
kurt_bsp = bootstrp(B,@kurtosis,data);

% Maximum skewness for plotting
if setxy_lim == 0
    sksq_max = max(2*skew_stat^2,12);
    sksq_lim = [-0.25 sksq_max];
    kurt_min = 0; kurt_max = max(2*kurt_stat,19);
    kurt_lim = [kurt_min kurt_max];
end

% Determine points corresponding to standard distributions
sk_std = [];
sk_std(1,:) = [0 3];    % Gaussian
sk_std(2,:) = [0 1.8]; % Uniform
sk_std(3,:) = [0 4.2]; % Logistic
sk_std(4,:) = [4 9];    % Exponential

% Squared skewness and Kurtosis for Gamma distribution
sksq_max = sksq_lim(2);
sksqgamma_vec = (0:0.1:sksq_max)';
kugamma_vec = 1.5*sksqgamma_vec + 3;

% Squared skewness and Kurtosis for lognormal distribution
sigma = (0.01:0.01:1)';
sksq_logn = (exp(sigma.^2) - 1).*(exp(sigma.^2) + 2).^2;
ind_vec = find(sksq_logn <= sksq_max);
kuvec_logn = exp(4*sigma.^2) + 2*exp(3*sigma.^2) + 3*exp(2*sigma.^2) - 3;

% Determine skewness and kurtosis bounds for Beta distribution
kurtvec1 = sksqgamma_vec +1;
kurtvec2 = 1.5*sksqgamma_vec + 3;

% Draw closed polygon
sksq_min = min(sksqgamma_vec); 
kurt_min1 = kurtvec1(1); kurt_min2 = kurtvec2(1);
kurt_max1 = kurtvec1(end); kurt_max2 = kurtvec2(end);

% Set up the figure for drawing
figure
set(gcf,'Color',[1 1 1]);


% Draw closed polygon for beta distribution
X = [sksq_min sksq_min sksq_max sksq_max sksq_min];
Y = [kurt_min1 kurt_min2 kurt_max2 kurt_max1 kurt_min1];
h_fill = fill(X,Y,[0.8 0.8 0.8],'FaceAlpha',0.5);
h_fill.EdgeColor = 'none';
set(gca,'fontsize',12,'fontweight','bold');
hold on
grid on
box off

% Figure adjustments
set(gca,'Xlim',sksq_lim,'Ylim',kurt_lim);
set(gca,'YDir','reverse');
xlabel('Squared Skewness');
ylabel('Kurtosis');

plot(sk_std(1,1),sk_std(1,2),'*','MarkerSize',10,'linewidth',2);
plot(sk_std(2,1),sk_std(2,2),'d','Color',[0 0.4470 0.7410],'MarkerSize',10,'linewidth',2);
plot(sk_std(3,1),sk_std(3,2),'x','MarkerSize',10,'linewidth',2);
plot(sk_std(4,1),sk_std(4,2),'+','MarkerSize',10,'linewidth',2);

plot(sksqgamma_vec,kugamma_vec,'--','linewidth',2);
plot(sksq_logn(ind_vec),kuvec_logn(ind_vec),':','linewidth',2);

% Plot the statistic and bootstrapped statistics
plot(skew_stat^2,kurt_stat,'ko','MarkerSize',10,'linewidth',2,'MarkerFaceColor','black');
s = scatter(skew_bsp.^2,kurt_bsp,10,[0.9290 0.6940 0.1250],'filled');
s.MarkerFaceAlpha = 0.5;

% Legend
h_leg = legend({'Beta'; 'Gaussian'; 'Uniform'; 'Logistic'; 'Exponential'; 'Gamma'; 'Lognormal';'Observed' ; 'Bootstraps';});
legpos = h_leg.Position;
h_text = text();
h_text.Units = 'normalized';
h_text.Position = [0.95*legpos(1) 0.95*legpos(2) 0];
h_text.String = 'Weibull is close to Gamma and Lognormal';
h_text.FontSize = 9;

% Title
title('Cullen-Frey Graph','fontsize',13,'fontweight','bold');

% Return the output
stats_sk = struct('Skewness',skew_stat,'Kurtosis',kurt_stat,'bootstat',[skew_bsp kurt_bsp]);
