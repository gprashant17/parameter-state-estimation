function y = liqlmeas_fun(x)

% Measurement function for liquid level system
% To be used with EKF

% For CH5115: Parameter and State Estimation
% Arun K. Tangirala
% January 08, 2021


% Get parameters from Model Workspace
% coder.e   xtrinsic('get_param');
% coder.extrinsic('get');
% 
% hws = get_param('extkalmfilt_simdemo','modelworkspace')

%Cv = hws.getVariable('Cv');
Cv = 0.8;
if x(2)<0
    x(2) = 0;
end
y = Cv*sqrt(x(2));
