function [SIGMA] = TVAR_Sigma(draw,priors,data)

Y = data.Yvar; % T x N
X = data.Xvar; % T x (No_coef/N)

e = Y - X * draw.phivar'; % draw.phi is N x (No_coef/N)
dof = size(Y,1) + priors.t0;
scale = e'*e + priors.Q0;

SIGMA = inv(wishrnd(inv(scale),dof));

end