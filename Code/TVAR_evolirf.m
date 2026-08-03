function [avg_evol] = TVAR_evolirf(draw,priors,data)

N = priors.N;
d = priors.d;
p = priors.Lag;
H = priors.horizon;
B = 200;
yt = data.ylevel(1,N-1); % GDP level at t
gt = data.ylevel(1,1); % Govt Spending comp level at t
ylevel = data.ylevel; % VAR variables in level at t (1 x N)

c = data.deutil_avg; % threshold value

phi1 = draw.phi1; % Low regime N x (1+Np)
phi2 = draw.phi2; % High regime N x (1+Np)

R1 = chol(draw.Q1); % Low regime N x N
R2 = chol(draw.Q2); % High regime N x N

base = zeros(H,N,B);
shock = zeros(H,N,B);
girf = zeros(H,N,B);
cumul_girf = zeros(H,N,B);
level_girf = zeros(H,N,B);

q_idx = N; % threshold variable is last variable
q_col = 1 + (d - 1)*N + q_idx; % column in X for q_{t-d}

if data.X(1,q_col) <= c
    R = R1;
else
    R = R2;
end

u1_shock = log(1 + 0.01*yt/gt)/R(1,1);
shock_dollars = log(1 + 0.01*yt/gt) * gt;

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
        
        qb = Xb(1,q_col);
        
        if qb <= c
            
            yhat_b = Xb*phi1'; % 1 x N Conditional Mean
            eps_b = u_draw(h,:)*R1; 
            
        else
            
            yhat_b = Xb*phi2'; % 1 x N Conditional Mean
            eps_b = u_draw(h,:)*R2; 
            
        end
        
        y_b = yhat_b + eps_b;
        
        base(h,:,b) = y_b;
        
        update_lags = Xb(:,2:1+N*(p-1));
        Xb = [1 y_b update_lags]; % Keep updating X as we move along H
    
% Shock simulation
    
        qs = Xs(1,q_col);
    
        if qs <= c
            
            yhat_s = Xs*phi1'; % 1 x N Conditional Mean
            eps_s = u_shock(h,:)*R1;
        
        else
            
            yhat_s = Xs*phi2'; % 1 x N Conditional Mean
            eps_s = u_shock(h,:)*R2;
            
        end
            
        y_s = yhat_s + eps_s;
        
        shock(h,:,b) = y_s;
        
        update_lags = Xs(:,2:1+N*(p-1));
        Xs = [1 y_s update_lags]; % Keep updating X as we move along H    
        
    end
    
    girf(:,:,b) = shock(:,:,b)-base(:,:,b);
    cumul_girf(:,:,b) = cumsum(girf(:,:,b));
    
    for n = [1:N-3, N-1]
        level_girf(:,n,b) = (cumul_girf(:,n,b)*ylevel(1,n))/shock_dollars;
    end
    level_girf(:,N-2,b) = girf(:,N-2,b); % deficit/GDP - level, no conversion
    level_girf(:,N,b) = girf(:,N,b); % threshold variable - report as is
end

% GIRF = shocked path - baseline path
avg_evol = mean(level_girf, 3);

end

    