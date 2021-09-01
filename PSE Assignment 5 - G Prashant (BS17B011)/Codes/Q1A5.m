    load ltvdata1.mat

%% Part (b) - Optimizing for forgetting factor
%Forgetting factors
lambda_vals = linspace(0.95,1,10);
var_vals = [];
M = inf;

%0th instant input
u0 = mean(uk);

for lambda = lambda_vals 
    % RLS object for Model 1
    obj_p1 = recursiveLS (3,'InitialParameters' ,0, 'InitialParameterCovariance' ,0.1,'ForgettingFactor',lambda);
    % RLS object for Model 2
    obj_p2 = recursiveLS (3,'InitialParameters' ,0, 'InitialParameterCovariance' ,0.1, 'ForgettingFactor',lambda);

    thetaest_vec1 = [];
    switch_track = [];  %Tracking switch
    y_est = [];   %Measured 

    for k = 3:length(uk)
        % Regressor vector
        if k == 3
            X1 = [yk(k-1), uk(k-1), uk(k-2)];
            X2 = [yk(k-1), yk(k-2), u0];
        else
            X1 = [yk(k-1), uk(k-1), uk(k-2)];
            X2 = [yk(k-1), yk(k-2), uk(k-3)];
        end
        
        %Making a copy of RLS objects
        if isempty(obj_p1.ParameterCovariance)
            obj_p1_b4 =  recursiveLS (3,'InitialParameters' ,obj_p1.Parameters, 'InitialParameterCovariance' ,0.1,'ForgettingFactor',lambda);
        else
            obj_p1_b4 =  recursiveLS (3,'InitialParameters' ,obj_p1.Parameters, 'InitialParameterCovariance' ,obj_p1.ParameterCovariance,'ForgettingFactor',lambda);
        end 
       if isempty(obj_p2.ParameterCovariance)
            obj_p2_b4 = recursiveLS (3,'InitialParameters' ,obj_p2.Parameters, 'InitialParameterCovariance' ,0.1,'ForgettingFactor',lambda);
       else
            obj_p2_b4 = recursiveLS (3,'InitialParameters' ,obj_p2.Parameters, 'InitialParameterCovariance' ,obj_p2.ParameterCovariance,'ForgettingFactor',lambda);
       end
       
       %Storing current value of b2
       if isempty(obj_p1.Parameters)
        b2 = obj_p1.InitialParameters(3);
       else
        b2 = obj_p1.Parameters(3);
       end
       
        %Output variables for both models
        yk1 = yk(k);
        yk2 = yk(k) - b2*uk(k-2);
        %Updating parameter for both models
        [theta1 ,ykhat1] = obj_p1(yk1,X1);
        [theta2 ,ykhat2] = obj_p2(yk2,X2);
        %Calculating residuals for both models
        e1 = (ykhat1-yk1)^2;
        e2 = (ykhat2-yk2)^2;
        
        %Noting model that gave minimum residual
        if e1<e2
            switch_track(k-2) = 1;
            y_est(k-2) = ykhat1;
            obj_p2 = obj_p2_b4;
        else
            switch_track(k-2) = 2;
            y_est(k-2) = ykhat2+b2*uk(k-2);
            obj_p1 = obj_p1_b4;
        end
    end
    
    % Estimating variance in parameter estimates
    var_ = var((y_est'-yk(3:end)));
    var_vals = [var_vals var_];
    
    %Best Lambda
    if var_<M
        M = var_;
        lambda_best = lambda;
        theta1_best = theta1;
        theta2_best = theta2;
        switch_best = switch_track;
        y_est_best = y_est;
    end
end

%Plotting
figure;
plot(lambda_vals,var_vals);
xlabel('lambda')
ylabel('Variance')
title('Variance in parameter estimates vs Forgetting Factor')
figure;
plot(3:length(uk),yk(3:end))
hold on
plot(3:length(uk),y_est_best,'-')
xlabel('Time instant k')
ylabel('y[k]')
title('y[k] - True and estimated')
legend(["True",['Estimate']])
figure;
plot(3:length(uk),switch_best);
xlabel('Time instant k')
ylabel('Model State')
title('Switching Time Identification')
%% Part (c) - Optimizing for forgetting factor
% clear;
%Window lengths
win_lengths = [10, 50,100,150,200,250,300,400,500,1000];
var_vals = [];
sw_var_vals = [];
M = inf;

%0th instant input
u0 = mean(uk);

%Initial Regressors for both models
X_m1 = zeros(length(uk)-2,3);
X_m1(:,1) = yk(2:(end-1));
X_m1(:,2) = uk(2:(end-1));
X_m1(:,3) = uk(1:(end-2));
X_m2 = zeros(length(uk)-2,3);
X_m2(:,1) = yk(2:(end-1));
X_m2(:,2) = yk(1:(end-2));
X_m2(2:end,3) = uk(1:(end-3));
X_m2(1,3) = u0;

for win_size = win_lengths 
    % RLS object for Model 1
    obj_p1 = recursiveLS (3,'History','Finite','WindowLength',win_size,'InitialOutputs',yk(3:(win_size+2)),'InitialRegressors',X_m1(1:win_size,:));
    % RLS object for Model 2
    obj_p2 = recursiveLS (3,'History','Finite','WindowLength',win_size,'InitialOutputs',yk(3:(win_size+2)),'InitialRegressors',X_m2(1:win_size,:));

    thetaest_vec1 = [];
    switch_track = [];
    y_est = [];
    
    % Regressor vector
    for k = 3:length(uk)
        if k == 3
            X1 = [yk(k-1), uk(k-1), uk(k-2)];
            X2 = [yk(k-1), yk(k-2), u0];
        else
            X1 = [yk(k-1), uk(k-1), uk(k-2)];
            X2 = [yk(k-1), yk(k-2), uk(k-3)];
        end
        
        %Making a copy of RLS objects
        obj_p1_b4 =  recursiveLS (3,'InitialParameters' ,obj_p1.Parameters,'History','Finite','WindowLength',win_size,'InitialOutputs',obj_p1.InitialOutputs,'InitialRegressors',obj_p1.InitialRegressors);
        obj_p2_b4 = recursiveLS (3,'InitialParameters' ,obj_p2.Parameters,'History','Finite','WindowLength',win_size,'InitialOutputs',obj_p2.InitialOutputs,'InitialRegressors',obj_p2.InitialRegressors);
       
        %Storing current value of b2
       if isempty(obj_p1.Parameters)
        b2 = obj_p1.InitialParameters(3);
       else
        b2 = obj_p1.Parameters(3);
       end
      
        yk1 = yk(k);
        yk2 = yk(k) - b2*uk(k-2);
        [theta1 ,ykhat1] = obj_p1(yk1,X1);
        [theta2 ,ykhat2] = obj_p2(yk2,X2);

        e1 = (ykhat1-yk1)^2;
        e2 = (ykhat2-yk2)^2;

        if e1<e2
            switch_track(k-2) = 1;
            y_est(k-2) = ykhat1;
            obj_p2 = obj_p2_b4;
        else
            switch_track(k-2) = 2;
            y_est(k-2) = ykhat2+b2*uk(k-2);
            obj_p1 = obj_p1_b4;
        end
    end
    var_ = var((y_est'-yk(3:end)));
    var_vals = [var_vals var_];
    var_sw = var(switch_track);
    sw_var_vals = [sw_var_vals var_sw];
    metric = var_sw*var_/(var_sw+var_);
    if metric<M
        M = metric;
        var_best = var_;
        win_best = win_size;
        theta1_best = theta1;
        theta2_best = theta2;
        switch_best = switch_track;
        y_est_best = y_est;
    end
end 

figure;
plot(win_lengths,var_vals)
hold on
plot(win_lengths,sw_var_vals)
xlabel('Window Size')
ylabel('Variance')
title('Optimizing window length')
legend('Parameter variance','switch variance')
figure;
plot(3:length(uk),yk(3:end))
hold on
plot(3:length(uk),y_est_best,'-')
xlabel('Time instant k')
ylabel('y[k]')
title('y[k] - True and estimated')
legend(["True",['Estimate']])
figure;
plot(3:length(uk),switch_best);
xlabel('Time instant k')
ylabel('Model State')
title('Switching Time Identification')
