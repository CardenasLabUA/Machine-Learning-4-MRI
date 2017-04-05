function [FittingPars,Signal_hat] =  fit_MSME(TE,Signal)
% Estimated T2 relaxation time using variable TR data
%
%% Syntaxis
%
%   [FittingPars,Ypred,Yci,vtrModel] =  fit_MSME(TE,Signal)
%
%% Inputs
%   TE= N x 1 vector of repetition time for VTR experiment. Units= seconds
%   Signal= N x 1 Matrix of signal for the VTR experiment. Units= AU
%
%% Ouputs
%   FittingPars= 5 X 1 numeric array with the following elements
%       FittingPars(1)= T2 time (estimated)
%       FittingPars(2)= standard eror of T2 time
%       FittingPars(3)= Lower 95% confience interval for the estimated T2
%       FittingPars(4)= Upper 95% confience interval for the estimated T2
%       FittingPars(5)= Adjusted Rsquared for the model 
%%  Author
% Julio Cárdenas-Rodríguez
% University of Arizona
% Tucson, AZ
% cardenaslab.org


% Adjust orientation
p=size(TE);
if p(1) < p(2)
TE=TE';
end

p=size(Signal);
if p(1) < p(2)
Signal=Signal';
end

Signal=Signal./max(Signal);

% Initial Guess
x0=[1.1,2.0];
T2vTR_func=@(pars,xdata) pars(1) * exp (-xdata./pars(2)) ;
lb=[0.5,0.010]';
ub=[1.5,1]';

% Fit variable TR model
options = optimoptions('lsqcurvefit');
options.Display='off';
[Coefficients,~,residual,~,~,~,jacobian]=lsqcurvefit(T2vTR_func,x0,TE,Signal,lb,ub,options);

conf95 = nlparci(Coefficients,residual,'jacobian',jacobian);


% Allocate Variables

T2pred=Coefficients(2);
T2ci=conf95(2,:);
Signal_hat=T2vTR_func(Coefficients,TE);
Rsquared=rsquare(Signal,Signal_hat);

FittingPars=[T2pred,T2ci,Rsquared]'; 