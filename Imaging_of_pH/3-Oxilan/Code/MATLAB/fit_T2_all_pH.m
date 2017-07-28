function [T2,T2_error,T2_rsq,Signal,TE]=fit_T2_all_pH( folder_path )
%% T2
% [Results,Signal,TE]=fit_T2( folder_path )
%
here=pwd;
cd(folder_path)
pH_folders=dir('pH*');
n_pH=length(pH_folders);

% Extract CEST data
last_n=0;
for k=1:n_pH
    cd([folder_path,pH_folders(k).name,'\T2']);
    TE=importdata('TE_ms.txt')/1000;
    %TE=TE(2:end);
    conc_files=dir('*mM*.txt');
    n_conc=length(conc_files);

     for q=1:n_conc
        s= importdata(conc_files(q).name);
        Signal(q+last_n,:)=s;
        Fittedpars=fit_MSME(TE,s./max(s))'; %#ok<*AGROW>
        T2(q+last_n,1)=Fittedpars(1);
        T2_error(q+last_n,:)=Fittedpars(2:3);
        T2_rsq(q+last_n,:)=Fittedpars(4);
     end
     last_n=last_n+n_conc;
 end
 

%Results=array2table([Concentration,Fittedpars],'VariableNames',{'Concentration','Estimated_T2','Lower_Bound','Upper_Bound','Rsquared'});

 cd(here)
% %% plot T2 Error bar
% figure()
% y=Fittedpars(:,1);
% ub=Fittedpars(:,3)-Fittedpars(:,1);
% lb=Fittedpars(:,1)-Fittedpars(:,2);
%  errorbar(Concentration,y,lb,ub,'-o','LineWidth',2,'MarkerSize',10);
%  xlabel('Concentration','FontSize',20');
%   ylabel('Estimated T2 (seconds)','FontSize',20');
  
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


function [r2, rmse] = rsquare(y,f,varargin)
% Compute coefficient of determination of data fit model and RMSE
%
% [r2 rmse] = rsquare(y,f)
% [r2 rmse] = rsquare(y,f,c)
%
% RSQUARE computes the coefficient of determination (R-square) value from
% actual data Y and model data F. The code uses a general version of 
% R-square, based on comparing the variability of the estimation errors 
% with the variability of the original values. RSQUARE also outputs the
% root mean squared error (RMSE) for the user's convenience.
%
% Note: RSQUARE ignores comparisons involving NaN values.
% 
% INPUTS
%   Y       : Actual data
%   F       : Model fit
%
% OPTION
%   C       : Constant term in model
%             R-square may be a questionable measure of fit when no
%             constant term is included in the model.
%   [DEFAULT] TRUE : Use traditional R-square computation
%            FALSE : Uses alternate R-square computation for model
%                    without constant term [R2 = 1 - NORM(Y-F)/NORM(Y)]
%
% OUTPUT 
%   R2      : Coefficient of determination
%   RMSE    : Root mean squared error
%
% EXAMPLE
%   x = 0:0.1:10;
%   y = 2.*x + 1 + randn(size(x));
%   p = polyfit(x,y,1);
%   f = polyval(p,x);
%   [r2 rmse] = rsquare(y,f);
%   figure; plot(x,y,'b-');
%   hold on; plot(x,f,'r-');
%   title(strcat(['R2 = ' num2str(r2) '; RMSE = ' num2str(rmse)]))
%   
% Jered R Wells
% 11/17/11
% jered [dot] wells [at] duke [dot] edu
%
% v1.2 (02/14/2012)
%
% Thanks to John D'Errico for useful comments and insight which has helped
% to improve this code. His code POLYFITN was consulted in the inclusion of
% the C-option (REF. File ID: #34765).

if isempty(varargin); c = true; 
elseif length(varargin)>1; error 'Too many input arguments';
elseif ~islogical(varargin{1}); error 'C must be logical (TRUE||FALSE)'
else c = varargin{1}; 
end

% Compare inputs
if ~all(size(y)==size(f)); error 'Y and F must be the same size'; end

% Check for NaN
tmp = ~or(isnan(y),isnan(f));
y = y(tmp);
f = f(tmp);

if c; r2 = max(0,1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
else r2 = 1 - sum((y(:)-f(:)).^2)/sum((y(:)).^2);
    if r2<0
    % http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
        warning('Consider adding a constant term to your model') %#ok<WNTAG>
        r2 = 0;
    end
end

rmse = sqrt(mean((y(:) - f(:)).^2));
end

end

  
  
end
 



 
 


