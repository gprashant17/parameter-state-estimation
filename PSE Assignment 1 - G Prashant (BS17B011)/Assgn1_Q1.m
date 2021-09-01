%% Integrating e^(-y).e^(-1/y) from 0.2 to 0.4
format long
func = @(y) exp(-y).*exp(-1./y);
integral(func,0.2,0.4)