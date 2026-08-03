function [avg_fixed] = TVAR_fixedirf(draw,priors,data)

N = priors.N;
p = priors.Lag;
H = priors.horizon;
B = 200;
%yt = data.ylevel(1,N-1); % GDP level at t
%gt = data.ylevel(1,1); % Govt Spending comp level at t
ylevel = data.ylevel; % Real VAR variables in level at t (1 x N)

phi = draw.phi; % N x (1+Np)

R = chol(draw.Q); % N x N

base = zeros(H,N,B);
shock = zeros(H,N,B);
girf = zeros(H,N,B);
cumul_girf = zeros(H,N,B);
level_girf = zeros(H,N,B);

%u1_shock = (1*yt/gt)/R(1,1);
%shock_dollars = (1*yt/gt) * gt;

u1_shock = 1/R(1,1);

% Baseline simulation

for b = 1:B
    
    % fixed structural shock path
    u_draw = randn(H,N);
    u_shock = u_draw;
    u_shock(1,1) = u1_shock; % replace first structural gov shock
    
    % separate histories for baseline and shocked paths
    Xb = data.X;
    Xs = data.X;
    
    for h = 1:H

        yhat_b = Xb*phi'; % 1 x N Conditional Mean
        eps_b = u_draw(h,:)*R; 
        y_b = yhat_b + eps_b;
    
        base(h,:,b) = y_b;
    
        update_lags = Xb(:,2:1+N*(p-1));
        Xb = [1 y_b update_lags]; % Keep updating X as we move along H
    
% Shock simulation
        
        yhat_s = Xs*phi'; % 1 x N Conditional Mean
        eps_s = u_shock(h,:)*R; % Same draw of H x N reduced-form shocks
        
        y_s = yhat_s + eps_s;
    
        shock(h,:,b) = y_s;
    
        update_lags = Xs(:,2:1+N*(p-1));
        Xs = [1 y_s update_lags]; % Keep updating X as we move along H
        
    end
    
    girf(:,:,b) = shock(:,:,b)-base(:,:,b);
    cumul_girf(:,:,b) = cumsum(girf(:,:,b));
    
    for n = 1:N-1
        level_girf(:,n,b) = (cumul_girf(:,n,b)*ylevel(1,n))/ylevel(1,1);
    end
    %level_girf(:,N-2,b) = girf(:,N-2,b); % deficit/GDP - level, no conversion
    level_girf(:,N,b) = girf(:,N,b); % threshold variable - report as is
end

% GIRF = shocked path - baseline path
avg_fixed = mean(level_girf, 3);

end


    