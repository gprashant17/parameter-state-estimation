function x = liqlstatetrans_fun(x,u)

% Discrete-time approximation of the state transition function 
% for liquid level system
% Example state transition function for discrete-time nonlinear state
% estimators.
%
% Inputs:
%    x - Present state x[k]
%    u - Present input
%
% Outputs:
%   Propagated states x[k+1]
%
% See also extendedKalmanFilter, unscentedKalmanFilter

% Get parameters from Model Workspace
% coder.extrinsic('get_param')
% coder.extrinsic('get');
% 
% hws = get_param('extkalmfilt_simdemo','modelworkspace');

% Qtrue = hws.getVariable('Qtrue');
%Ts = hws.getVariable('Ts');
%Cv = hws.getVariable('Cv');
%Ac = hws.getVariable('Ac');

Ac = 1.2; Cv = 0.8; Ts = 0.1;

% Euler integration of continuous-time dynamics x'=f(x) with sample time Ts
x(1) = x(1) + LiqLevelContinuous1(x(1),u,Cv,Ac)*Ts;
x(2) = x(2) + LiqLevelContinuous2(x,Cv,Ac)*Ts;

end
function dx1dt = LiqLevelContinuous1(x,u,Cv,Ac)
% ODE for liquid level dynamics
if x < 0
    x = 0;
end
dx1dt = (u/Ac) - (Cv/Ac)*sqrt(x);
end

function dx2dt = LiqLevelContinuous2(x,Cv,Ac)

% ODE for liquid level dynamics
if x(2) < 0
    x(2) = 0;
end

if x(1) < 0
    x(1) = 0;
end

dx2dt = (Cv/Ac)*sqrt(x(1)) - (Cv/Ac)*sqrt(x(2));
end