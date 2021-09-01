%% Assigning values to parameters
a1 = 1.2; a2 = 1.2;
cv1 = 0.8;cv2 = 0.8;
Fis = 2;
Ts = 0.1;
%% Part (b) - Discretize model

% Continuous model
A = [-cv1/a1 0 ; cv1/a1 -cv2/a2];
B = [1/a1 ; 0];
C = [0 cv2];
D = 0;
sys_con = ss(A,B,C,D,'InputName',{'u'});

%Discretizing model with sampling interval of Ts = 0.1
sys_dis = c2d(sys_con,Ts)

%% Part (d) - Simulating data with noise
%Generating input output data
N = 1275;
rng(2121)
ukvec = idinput(N,'prbs',[0 0.2],[-1 1]); %input

Q = diag([0.2 0.1]); %Process noise variance
R = 0.1; %Measurement Noise Variance
%Process Noise Generation
rng(712)
wkvec1 = sqrt(0.2)*randn(N,1);  
rng(362)
wkvec2 = sqrt(0.1)*randn(N,1);
%Measurement Noise Generation
rng(213)
vkvec = sqrt(R)*randn(N,1);

sys_noise = ss(sys_dis.A,[sys_dis.B , [1 ; 0] ,[0 ; 1] zeros(2,1)],sys_dis.C,[0 0 0 1],Ts,'inputname',{'u' 'w1','w2','v'},'outputname',{'ym'});
[ykvec,kvec,xkmat] = lsim(sys_noise,[ukvec wkvec1 wkvec2 vkvec]);

%Plotting
figure;
plot(kvec,ykvec);
xlabel('Time Instant k')
ylabel('y[k]');
title('y[k]')
figure;
subplot(2,1,1)
plot(kvec,xkmat(:,1));
title('State h1')
xlabel('Time instant k');
ylabel('x1[k]');
subplot(2,1,2)
plot(kvec,xkmat(:,2));
title('State h2')
xlabel('Time instant k');
ylabel('x2[k]');

%% Part (e) - Determining optimal Q
rng(2020);
Sys = ss(sys_dis.A,[sys_dis.B , [1 ; 0] ,[0 ; 1]],sys_dis.C,[0 0 0],Ts,'inputname',{'u' 'w1','w2'},'outputname',{'ym'});
Q_vals1 = linspace(0,0.5,51);
Q_vals2 = linspace(0,0.5,51);

y_var = [];
x_var = [];
Rf = 1;
for Qf1 = Q_vals1
    for Qf2 = Q_vals2
        Qf = [Qf1 0 ; 0 Qf2];
        % Parallel Kalman Filter
        [kalmf,Kp,P,Kf,Pf,Ky] = kalman(Sys,Qf,Rf); % Kp and Kf, the predictor and filter Kalman gain
        Process1 = ss(sys_dis.A,[sys_dis.B , [1 ; 0] ,[0 ; 1] zeros(2,1)],sys_dis.C,[0 0 0 1],0.1,'inputname',{'u' 'w1','w2','v'},'outputname',{'ym'},'statename',{'x1','x2'});
        Process1Kf = parallel(Process1,kalmf,1,1,[],[]);
        SimTotal1 = feedback(Process1Kf,1,5,1,1);
        SimTotal1 = SimTotal1((1:4),[1 2 3 4]);
        [outmat1,kvec] = lsim(SimTotal1,[wkvec1 wkvec2 vkvec ukvec]);

        var_xkhat = sum(var(outmat1(:,3:4) - xkmat(:,1:2)));
        var_ykhat = var(outmat1(:,1) - outmat1(:,2));
        x_var = [x_var var_xkhat];
        y_var = [y_var var_ykhat];
    end
end

%Plotting
x_var = reshape(x_var,51,51)';
y_var = reshape(y_var,51,51)';
figure;
surf(Q_vals1,Q_vals2,x_var);
xlabel('Q1');
ylabel('Q2');
zlabel('Q3');
title('Variance in State estimates - Tuning Q')
figure;
surf(Q_vals1,Q_vals2,y_var);
xlabel('Q1');
ylabel('Q2');
zlabel('Q3');
title('Variance in measurement estimates - Tuning Q')

[Q1i,Q2j] = find(x_var==min(min(x_var)));
Q1best = Q_vals1(Q1i);
Q2best = Q_vals2(Q2j);