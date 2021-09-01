%% State space system
A = [0.9 0 0 ; 1 1.2 -0.5916 ; 0 0.5916 0];
B = [1;0;0];
C = [2 0.8 -0.6761];
D = 0;

%% Observability Matrix
O = [C;C*A;C*A^2];
% rank(O)

%% Controllability Matrix
Con = [B, A*B,(A^2)*B];
% rank(Con)

%% Minimal Realization
[ABAR,BBAR,CBAR,T,K] = obsvf(A,B,C);
A0 = ABAR(2:end,2:end);
B0 = BBAR(2:end);
C0 = CBAR(2:end);
D0 = 0;

%% Part (b) - Kalman Gain
p = [0.2 0.1];
K = place(A0',C0',p)';
%% Part (c) - Designing and evaluating Observer without noise
%Choosing the eigenvalues as 0.2 and 0.1
p = [0.2 0.1];
K = place(A0',C0',p)';
plant_nf = ss(A0,B0,C0,D0,1,'InputName','u');
est_plant = estim(plant_nf,K,1,1);

%True Data Generation
N = 2046;
rng(1000);
ukvec = idinput(N,'prbs',[0 0.4],[-1 1]);
orig_plant = ss(A,B,C,D,1,'InputName','u');
[ykvec,~,~] = lsim(orig_plant,[ukvec]);

%Estimation using minimal realization
[ykhat,~,xk] = lsim(est_plant,[ukvec ykvec]);
plot(1:length(ykvec),ykvec);
hold on
plot(1:length(ykvec),ykhat(:,1),'-');
xlabel('Time instant k');
ylabel('y[k]');
title('True and Predicted y[k]');
legend(['True'],['Predicted']);
xlim([0,2050])
var_ = var(ykvec - ykhat(:,1))
%% Part (d) - Designing and evaluating Observer with noise
%Choosing the eigenvalues as 0.2 and 0.1
 
plant_nf = ss(A0,B0,C0,D0,1,'InputName','u');
est_plant = estim(plant_nf,K,1,1);

%True Data Generation
N = 2046;
Q = 0.1 ; R = 1;
rng(2131)
wkvec = randn(N,1)*sqrt(Q);
rng(21)
vkvec = randn(N,1)*sqrt(R);
rng(1000)
ukvec = idinput(N,'prbs',[0 0.4],[-1 1]);
orig_plant = ss(A,[B ones(3,1) zeros(3,1)],C,[0 0 1],1,'InputName',{'u','w','v'},'OutputName',{'ym'});
[ykvec,~,~] = lsim(orig_plant,[ukvec,wkvec,vkvec]);

%Estimation using minimal realization
[ykhat,~,xkhat] = lsim(est_plant,[ukvec ykvec]);
plot(1:length(ykvec),ykvec);
hold on
plot(1:length(ykvec),ykhat(:,1));
xlabel('Time instant k');
ylabel('y[k]');
title('True and Predicted y[k]');
legend(['True'],['Predicted']);
xlim([0,2050])
var_state = sum(var(xk - xkhat))
var_meas = var(ykvec - ykhat(:,1))


subplot(2,1,1)
plot(1:length(ykvec),xk(:,1));
hold on
plot(1:length(ykvec),xkhat(:,1));
xlim([0,2050])
title('State x1')
xlabel('Time instant k');
ylabel('x1[k]');
legend(['True'],['Predicted']);

subplot(2,1,2)
plot(1:length(ykvec),xk(:,2));
hold on
plot(1:length(ykvec),xkhat(:,2));
xlim([0,2050])
title('State x2')
xlabel('Time instant k');
ylabel('x2[k]');
legend(['True'],['Predicted']);

%% Part (d) - Eigenvalue analysis
%Choosing the eigenvalues - change values to check variance 
p = [0.9 0.1];
%Placing A-KC
K = place(A0',C0',p)';
Q = 0.1 ; R = 1;
rng(2131)
wkvec = randn(N,1)*sqrt(Q);
rng(21)
vkvec = randn(N,1)*sqrt(R);
plant_wn = ss(A0,[B0 ones(2,1) zeros(2,1)],C0,[0 0 1],1,'InputName',{'u','w','v'},'OutputName',{'ym'});
[yk,~,xk] = lsim(plant_wn,[ukvec wkvec vkvec]);

est_plant_wn = estim(plant_wn,K,1,1);
[ykhat,~,xkhat] = lsim(est_plant_wn,[ukvec ykvec]);

subplot(2,1,1)
plot(1:length(ykvec),xk(:,1));
hold on
plot(1:length(ykvec),xkhat(:,1));
xlim([0,2050])
title('State x1')
xlabel('Time instant k');
ylabel('x1[k]');
var_x1 = var((xk(:,1) - xkhat(:,1)))
subplot(2,1,2)
plot(1:length(ykvec),xk(:,2));
hold on
plot(1:length(ykvec),xkhat(:,2));
xlim([0,2050])
var_x2 = var((xk(:,2) - xkhat(:,2)))
title('State x2')
xlabel('Time instant k');
ylabel('x2[k]');

var_x1