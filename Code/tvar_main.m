%% Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bayesian TVAR Example (Fiscal Policy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;

%% Bayesian TVAR

% Setup model parameters

priors.Lag = 4; % number of lags for FAVAR
priors.horizon = 21; % horizon for IRF
priors.keep = 15000; % Iterations to be considered
priors.burn = 15000; % Iterations to be discarded
priors.display = 1; % display iterations in command window
priors.K = 1; % No of govt spend components
priors.M = 3; % No of other macroeconomic variables in VAR (T,Y)
priors.N = priors.K + priors.M;
priors.constant = 1; % zero if no constant in the VAR
priors.FPS = 1; % Fiscal Policy Shock to the variable in VAR
priors.trend = 0; % Trend component
priors.d = 1; % Lag of capacity utilization

% Import Data

data.raw1 = readtable('Data_components_July.xlsx','Sheet',2);
data.raw = table2array(data.raw1(81:264,:)); 

% Transforming data

data.GC = data.raw(:,4         ); % Gov Cons
data.GI = data.raw(:,5); % Gov Inv
data.tax = data.raw(:,32); % tax rate
data.GDP = data.raw(:,12); % GDP
data.PGDP = data.raw(:,13); % GDP Def
data.POP = data.raw(:,14); % Population
data.GREV = data.raw(:,15); % Tax Revenue
data.FISREC = data.raw(:,39); % Federal Govt Receipts
data.FISEXP = data.raw(:,40); % Federal Govt Exp
data.FISDEF = (data.FISEXP - data.FISREC)./data.GDP; % Fiscal Def as % of GDP
data.UTIL = data.raw(:,60); % Capacity utilization
data.UNEMP = data.raw(:,61); % Unemployment rate
data.DEUTIL_RAW = data.raw(:,62); % Mean-adj Capacity utilization

% Defining growth rates

data.taxgrowth = diff(log(data.GREV./(data.PGDP))); % Growth rate real taxes
data.ggrowth = diff(log((data.GI+data.GC)./(data.PGDP))); % Growth rate real total gov spend
data.gdpgrowth = diff(log(data.GDP./(data.PGDP))); % Growth rate real GDP

data.rtax = data.GREV*100./data.PGDP;
data.rg = (data.GI+data.GC)*100./data.PGDP;
data.rgdp = data.GDP*100./data.PGDP;

%% Unit root diagnostics for multiple series (ADF + KPSS)

unitroot.alpha = 0.05;
unitroot.seriesList = {'ggrowth','taxgrowth','gdpgrowth','DEUTIL_RAW'};
unitroot.n = numel(unitroot.seriesList);

% Preallocate
unitroot.ADF_p  = NaN(unitroot.n,1);
unitroot.ADF_h  = false(unitroot.n,1);
unitroot.KPSS_p = NaN(unitroot.n,1);
unitroot.KPSS_h = false(unitroot.n,1);

for i = 1:unitroot.n
    unitroot.nm = unitroot.seriesList{i};

    % Clean numeric column vector (drop NaN/Inf)
    unitroot.y = data.(unitroot.nm);
    unitroot.y = unitroot.y(:);
    unitroot.y = unitroot.y(isfinite(unitroot.y));

    % ADF: H0 = unit root; h=1 => reject => stationary
    [unitroot.h_adf, unitroot.p_adf] = adftest(unitroot.y, 'Model','ARD', 'Alpha',unitroot.alpha);

    % KPSS: H0 = stationarity; h=0 => fail to reject => stationary
    [unitroot.h_kpss, unitroot.p_kpss] = kpsstest(unitroot.y, 'Trend',false, 'Alpha',unitroot.alpha);

    unitroot.ADF_h(i)  = unitroot.h_adf;
    unitroot.ADF_p(i)  = unitroot.p_adf;
    unitroot.KPSS_h(i) = unitroot.h_kpss;
    unitroot.KPSS_p(i) = unitroot.p_kpss;
end

unitroot.ADF_h  = unitroot.ADF_h(:);
unitroot.KPSS_h = unitroot.KPSS_h(:);

unitroot.stationary_adf  = repmat("No",  numel(unitroot.ADF_h),  1);
unitroot.stationary_kpss = repmat("No",  numel(unitroot.KPSS_h), 1);

unitroot.stationary_adf(unitroot.ADF_h == 1)     = "Yes";   % ADF: h=1 => stationary
unitroot.stationary_kpss(unitroot.KPSS_h == 0)   = "Yes";   % KPSS: h=0 => stationary

unitroot.unitroot_results = table( ...
    string(unitroot.seriesList(:)), unitroot.ADF_p(:), unitroot.KPSS_p(:), ...
    unitroot.stationary_adf, unitroot.stationary_kpss, ...
    'VariableNames', {'series','p_value_adf','p_value_kpss','stationary_adf','stationary_kpss'} );

%disp(unitroot.unitroot_results);

%% Setting up the data

% Defining the two regimes

data.DEUTIL = data.DEUTIL_RAW(1+priors.Lag:end-priors.d);
data.deutil_avg = mean(data.DEUTIL);
data.pct_above = sum(data.DEUTIL > data.deutil_avg)/size(data.DEUTIL,1);

data.regime1_idx = (data.DEUTIL <= data.deutil_avg);
data.regime2_idx = (data.DEUTIL > data.deutil_avg);

% Given the data from 1968Q1 to 2022Q4, approximately 55% of the observations are
% above the mean of the mean-adjusted capacity

%data.Yall = [data.defcgrowth data.defigrowth data.ndefcgrowth data.ndefigrowth data.slgg data.FISDEF(2:end,:) data.gdpgrowth data.DEUTIL_RAW(2:end,:)]; 
data.Yall = [data.ggrowth*100 data.taxgrowth*100 data.gdpgrowth*100 data.DEUTIL_RAW(2:end,:)]; 
data.Tall = size(data.Yall,1);
data.Xvar = ones(data.Tall,1);

for p = 1:priors.Lag
    data.Xvar = [data.Xvar,add_lag(data.Yall,p)];
end

data.Xvar = data.Xvar(1+priors.Lag:end,:); 
data.Yvar = data.Yall(1+priors.Lag:end,:); 
data.T = size(data.Xvar,1);

data.Yreg1 = data.Yvar(data.regime1_idx,:);
data.Xreg1 = data.Xvar(data.regime1_idx,:);

data.Yreg2 = data.Yvar(data.regime2_idx,:);
data.Xreg2 = data.Xvar(data.regime2_idx,:);

% 99 observations in regime 1 and 120 observations in regime 2

% Data for computing dollar values of IRFs

data.Ylevel = [data.rg data.rtax data.rgdp data.DEUTIL_RAW];
data.Ylevel = data.Ylevel(2+priors.Lag:end,:);
data.Ylevel_reg1 = data.Ylevel(data.regime1_idx,:);
data.Ylevel_reg2 = data.Ylevel(data.regime2_idx,:);

%% Bayesian Estimation (Minnesota Prior)

% No of coefficients (constant + 4 lags of four variables in each equation)

priors.no_coef = priors.N * (1 + priors.N * priors.Lag);

% Setting up starting values and priors - see notes (tvar.pdf)

priors.phi0 = zeros(priors.N,priors.no_coef/priors.N); % Following Bernanke et al (2005)
priors.phi0 = priors.phi0(:); % vector form

priors.index = zeros(priors.N,priors.Lag);
for i=1:priors.N
    priors.index(i,:) = priors.constant+i:priors.N:(priors.no_coef/priors.N);
end
priors.sigma_reg1 = zeros(1*priors.N,1);
priors.sigma_reg2 = zeros(1*priors.N,1);

for q = 1:priors.N
    data.y1 = data.Yreg1(:,q); % dependent variable
    data.x1 = data.Xreg1(:,[1 priors.index(q,:)]); % constant + lags of dependent variable only
    priors.areg1 = data.x1 \ data.y1; % regression coefficient
    priors.sigma_reg1(q,1) = ((data.y1 - data.x1*priors.areg1)'*(data.y1 - data.x1*priors.areg1))/(size(data.x1,1)-size(data.x1,2)); % error variance
end

for q = 1:priors.N
    data.y2 = data.Yreg2(:,q); % dependent variable
    data.x2 = data.Xreg2(:,[1 priors.index(q,:)]); % constant + lags of dependent variable only
    priors.areg2 = data.x2 \ data.y2; % regression coefficient
    priors.sigma_reg2(q,1) = ((data.y2 - data.x2*priors.areg2)'*(data.y2 - data.x2*priors.areg2))/(size(data.x2,1)-size(data.x2,2)); % error variance
end

priors.Q0_reg1 = diag(priors.sigma_reg1); % VAR var cov matrix for regime 1
priors.Q0_reg2 = diag(priors.sigma_reg2); % VAR var cov matrix for regime 1

% Note: With univariate AR regressions, K < T and therefore we can use the
% above formula for variance computation (e'e/T-K)

priors.t0 = priors.N+2;

% Compute the prior hyperparameter H (var cov matrix of coefficients)

priors.lambda1 = 1;
priors.lambda2 = 1;
priors.lambda3 = 1;
priors.lambda4 = 1;

priors.H_i1 = zeros(priors.no_coef/priors.N,priors.N);
priors.H_i2 = zeros(priors.no_coef/priors.N,priors.N);

% Create an array of dimensions (N+N^2*lag) x N, which will contain the (N+N^2*lag) diagonal
% elements of the covariance matrix, in each of the N equations.

for i = 1:priors.N  % for each i-th equation
    for j = 1:(priors.no_coef/priors.N)   % for each j-th RHS variable in ith equation
        if ismember(j, 1) % j==1 if no trend
            priors.H_i1(j,i) = (priors.lambda4^2)*priors.sigma_reg1(i,1); % variance on constant 
        elseif find(j==priors.index(i,:))>0
            ll = find(j==priors.index(i,:),2);
            priors.H_i1(j,i) = (priors.lambda1/(ll^priors.lambda3))^2; % variance on own lags
        else
            for kj = 1:priors.N
                if find(j==priors.index(kj,:))>0
                    mm = kj;
                    nn = find(j==priors.index(kj,:),2);
                end
            end
            priors.H_i1(j,i) = ((priors.lambda1*priors.lambda2)^2*priors.sigma_reg1(i,1))/(nn^(2*priors.lambda3)*priors.sigma_reg1(mm,1)); % variance on lags of other variables
        end
    end
end

for i = 1:priors.N  % for each i-th equation
    for j = 1:(priors.no_coef/priors.N)   % for each j-th RHS variable in ith equation
        if ismember(j, 1) % j==1 if no trend
            priors.H_i2(j,i) = (priors.lambda4^2)*priors.sigma_reg2(i,1); % variance on constant 
        elseif find(j==priors.index(i,:))>0
            ll = find(j==priors.index(i,:),2);
            priors.H_i2(j,i) = (priors.lambda1/(ll^priors.lambda3))^2; % variance on own lags
        else
            for kj = 1:priors.N
                if find(j==priors.index(kj,:))>0
                    mm = kj;
                    nn = find(j==priors.index(kj,:),2);
                end
            end
            priors.H_i2(j,i) = ((priors.lambda1*priors.lambda2)^2*priors.sigma_reg2(i,1))/(nn^(2*priors.lambda3)*priors.sigma_reg2(mm,1)); % variance on lags of other variables
        end
    end
end

% We get the prior variance covariance matrix H
priors.H0_reg1 = diag(priors.H_i1(:));
priors.H0_reg2 = diag(priors.H_i2(:));

% Initial values

draw.phivar_reg1 = zeros(priors.N,1+priors.N*priors.Lag);
priors.bols1 = (data.Xreg1\data.Yreg1)';
priors.u_reg1 = data.Yreg1 - data.Xreg1*priors.bols1';

draw.phivar_reg1 = priors.bols1;
draw.Qvar_reg1 = (priors.u_reg1'*priors.u_reg1)/(size(data.Xreg1,1)-size(data.Xreg1,2));

draw.phivar_reg2 = zeros(priors.N,1+priors.N*priors.Lag);
priors.bols2 = (data.Xreg2\data.Yreg2)';
priors.u_reg2 = data.Yreg2 - data.Xreg2*priors.bols2';

draw.phivar_reg2 = priors.bols2;
draw.Qvar_reg2 = (priors.u_reg2'*priors.u_reg2)/(size(data.Xreg2,1)-size(data.Xreg2,2));

% Storage of results

results.phi_reg1 = zeros(priors.N,(priors.no_coef/priors.N),priors.keep); % each phi drawn will be N x (1+N*lag)
results.Q_reg1 = zeros(priors.N,priors.N,priors.keep); % each SIGMA drawn will be N x N

results.phi_reg2 = zeros(priors.N,(priors.no_coef/priors.N),priors.keep); 
results.Q_reg2 = zeros(priors.N,priors.N,priors.keep); 

% Fixed States IRF

results.fixed_IRF_reg1 = zeros(priors.horizon,priors.N,priors.keep); 
results.fixed_IRF_reg2 = zeros(priors.horizon,priors.N,priors.keep);

% Fixed States GDP Multiplier

results.fixed_gdpmult_reg1 = zeros(priors.horizon,1,priors.keep); 
results.fixed_gdpmult_reg2 = zeros(priors.horizon,1,priors.keep);

% Evolving States IRF

results.evol_IRF_reg1 = zeros(priors.horizon,priors.N,priors.keep); 
results.evol_IRF_reg2 = zeros(priors.horizon,priors.N,priors.keep);

% Evol States GDP Multiplier

results.evol_gdpmult_reg1 = zeros(priors.horizon,1,priors.keep); 
results.evol_gdpmult_reg2 = zeros(priors.horizon,1,priors.keep);

%% Gibbs Sampling

% REGIME 1 

draw.phivar = draw.phivar_reg1;
draw.Qvar   = draw.Qvar_reg1;

priors.H0 = priors.H0_reg1;
priors.Q0 = priors.Q0_reg1;

data.Yvar = data.Yreg1;
data.Xvar = data.Xreg1;

c = 0;
for i = 1:priors.burn+priors.keep
    
    [draw.phivar] = TVAR_phi(draw,priors,data);
    [draw.Qvar] = TVAR_Sigma(draw,priors,data);
    
    if i > priors.burn
        c = c+1;
        results.phi_reg1(:,:,c) = draw.phivar;
        results.Q_reg1(:,:,c) = draw.Qvar;
    end
    
    if (mod(i,100)==0) && (priors.display == 1)
        disp(sprintf(' Regime 1: Replication %s of %s' , num2str(i), num2str(priors.burn+priors.keep)));       
    end
end

%% REGIME 2

draw.phivar = draw.phivar_reg2;
draw.Qvar   = draw.Qvar_reg2;

priors.H0 = priors.H0_reg2;
priors.Q0 = priors.Q0_reg2;

data.Yvar = data.Yreg2;
data.Xvar = data.Xreg2;

c = 0;
for i = 1:priors.burn+priors.keep
    
    [draw.phivar] = TVAR_phi(draw,priors,data);
    [draw.Qvar] = TVAR_Sigma(draw,priors,data);
    
    if i > priors.burn
        c = c+1;
        results.phi_reg2(:,:,c) = draw.phivar;
        results.Q_reg2(:,:,c) = draw.Qvar;
    end
    
    if (mod(i,100)==0) && (priors.display == 1)
        disp(sprintf(' Regime 2: Replication %s of %s' , num2str(i), num2str(priors.burn+priors.keep)));       
    end
end

%% Fixed States GIRF Computation

for i = 1:priors.keep
    for r = 1:2
        if r==1
            draw.phi = results.phi_reg1(:,:,i);
            draw.Q = results.Q_reg1(:,:,i);
            
            draw.fixed_low_all = zeros(priors.horizon,priors.N,size(data.Xreg1,1));
    
            for j = 1:size(data.Xreg1,1)
                data.X = data.Xreg1(j,:);
                data.ylevel = data.Ylevel_reg1(j,:);
                
                [draw.fixed_low] = TVAR_fixedirf(draw,priors,data);
                draw.fixed_low_all(:,:,j) = draw.fixed_low;
                
            end
    
            results.fixed_IRF_reg1(:,:,i) = mean(draw.fixed_low_all,3); % Average over diff histories
            results.fixed_gdpmult_reg1(:,:,i) = results.fixed_IRF_reg1(:,priors.N-1,i)./results.fixed_IRF_reg1(:,1,i);
        
        else
            draw.phi = results.phi_reg2(:,:,i);
            draw.Q = results.Q_reg2(:,:,i);
            
            draw.fixed_high_all = zeros(priors.horizon,priors.N,size(data.Xreg2,1));
            
            for j = 1:size(data.Xreg2,1)
                data.X = data.Xreg2(j,:);
                data.ylevel = data.Ylevel_reg2(j,:);
                
                [draw.fixed_high] = TVAR_fixedirf(draw,priors,data);
                draw.fixed_high_all(:,:,j) = draw.fixed_high;
                
            end
    
            results.fixed_IRF_reg2(:,:,i) = mean(draw.fixed_high_all,3); 
            results.fixed_gdpmult_reg2(:,:,i) = results.fixed_IRF_reg2(:,priors.N-1,i)./results.fixed_IRF_reg2(:,1,i);
            
        end
    end
    
    if (mod(i,100)==0) && (priors.display == 1)
        disp(sprintf(' Fixed IRF: Replication %s of %s' , num2str(i), num2str(priors.keep)));       
    end
end

%% Evolving States GIRF Computation

for i = 1:priors.keep
    
    draw.phi1 = results.phi_reg1(:,:,i);
    draw.phi2 = results.phi_reg2(:,:,i);
    draw.Q1 = results.Q_reg1(:,:,i);
    draw.Q2 = results.Q_reg2(:,:,i);
    
    for r = 1:2
        if r==1
            
            draw.evol_low_all = zeros(priors.horizon,priors.N,size(data.Xreg1,1));
    
            for j = 1:size(data.Xreg1,1)
                data.X = data.Xreg1(j,:);
                data.Y = data.Yreg1(j,:);
                data.ylevel = data.Ylevel_reg1(j,:);
                
                [draw.evol_low] = TVAR_evolirf(draw,priors,data);
                draw.evol_low_all(:,:,j) = draw.evol_low;
                
            end
    
            results.evol_IRF_reg1(:,:,i) = mean(draw.evol_low_all,3);
            results.evol_gdpmult_reg1(:,:,i) = results.evol_IRF_reg1(:,priors.N-1,i)./results.evol_IRF_reg1(:,1,i);
        
        else
            
            draw.evol_high_all = zeros(priors.horizon,priors.N,size(data.Xreg2,1));
           
            for j = 1:size(data.Xreg2,1)
                data.X = data.Xreg2(j,:);
                data.Y = data.Yreg2(j,:);
                data.ylevel = data.Ylevel_reg2(j,:);
                
                [draw.evol_high] = TVAR_evolirf(draw,priors,data);
                draw.evol_high_all(:,:,j) = draw.evol_high;
                
            end
    
            results.evol_IRF_reg2(:,:,i) = mean(draw.evol_high_all,3);
            results.evol_gdpmult_reg2(:,:,i) = results.evol_IRF_reg2(:,priors.N-1,i)./results.evol_IRF_reg2(:,1,i);
            
        end
    end
    
    if (mod(i,100)==0) && (priors.display == 1)
        disp(sprintf(' Evolving IRF: Replication %s of %s' , num2str(i), num2str(priors.keep)));       
    end
end

%% Computing the posterior median and confidence intervals

posterior.Q1 = prctile(results.Q_reg1,[50 16 84],3);
posterior.Q2 = prctile(results.Q_reg2,[50 16 84],3);

posterior.phi1 = prctile(results.phi_reg1,[50 16 84],3);
posterior.phi2 = prctile(results.phi_reg2,[50 16 84],3);

posterior.IRF_Fixed1 = prctile(results.fixed_IRF_reg1,[50 16 84],3);
posterior.IRF_Fixed2 = prctile(results.fixed_IRF_reg2,[50 16 84],3);

posterior.IRF_Evol1 = prctile(results.evol_IRF_reg1,[50 16 84],3);
posterior.IRF_Evol2 = prctile(results.evol_IRF_reg2,[50 16 84],3);

posterior.GDPmult_Fixed1 = prctile(results.fixed_gdpmult_reg1,[50 16 84],3);
posterior.GDPmult_Fixed2 = prctile(results.fixed_gdpmult_reg2,[50 16 84],3);

posterior.GDPmult_Evol1 = prctile(results.evol_gdpmult_reg1,[50 16 84],3);
posterior.GDPmult_Evol2 = prctile(results.evol_gdpmult_reg2,[50 16 84],3);

%% Graphing IRFs

% Plot the GIRF for same variable together

data.names = {'Government Spending','Taxes', 'GDP', 'Mean-Adjusted Capacity Utilization'};

[posterior.minmax_var] = TVAR_minmax(priors,posterior);

priors.green_line = [0 0.5 0];
priors.green_band = [0.6 0.85 0.6];

priors.red_line   = [0.8 0 0];
priors.red_band   = [1 0.7 0.7];

for i = 1:priors.N
    
    figure(i);
    clf;
    set(gcf,'Color','white');
    
    %---------------- Top row ----------------%
    subplot(2,2,1)
    tvar_plotirf(squeeze(posterior.IRF_Fixed1(:,i,:)), posterior.minmax_var(i,:), priors.green_line, priors.green_band);
    title('LOW STATE','FontWeight','normal','FontSize',16);
    ylabel('FIXED STATE RESPONSES','FontSize',12);
    
    subplot(2,2,2)
    tvar_plotirf(squeeze(posterior.IRF_Fixed2(:,i,:)), posterior.minmax_var(i,:), priors.red_line, priors.red_band);
    title('HIGH STATE','FontWeight','normal','FontSize',16);
    
    %--------------- Bottom row --------------%
    %subplot(2,2,3)
    %tvar_plotirf(squeeze(posterior.IRF_Evol1(:,i,:)), posterior.minmax_var(i,:), priors.green_line, priors.green_band);
    %ylabel({'EVOLVING STATE RESPONSES'; 'AVERAGE OVER ALL HISTORIES'},'FontSize',12);
    
    %subplot(2,2,4)
    %tvar_plotirf(squeeze(posterior.IRF_Evol2(:,i,:)), posterior.minmax_var(i,:), priors.red_line, priors.red_band);
    
    % Main title
    if iscell(data.names)
        sgtitle(data.names{i},'FontSize',20);
    else
        sgtitle(char(data.names(i)),'FontSize',20);
    end
    
end
