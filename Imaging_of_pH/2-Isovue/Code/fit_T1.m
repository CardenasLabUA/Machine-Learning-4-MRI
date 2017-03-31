function [Results,Signal,TR]=fit_T1( folder_path )
%% T1
here=pwd;
cd(folder_path)
file_names=dir('*mM*.txt');
number_of_samples=length(file_names);
TR=importdata('TR_ms.txt')./1000;

 for q=1:number_of_samples
    s= importdata(file_names(q).name);
    Signal(:,q)=s;
    Concentration(q,1)=str2double(file_names(q).name(1:5));
[FittingPars(:,q)] =  fit_variable_TR(TR,s./max(s)); %#ok<*AGROW>

 end

Estimated_T1=FittingPars(1,:)';
LowerBound=FittingPars(2,:)';
UpperBound=FittingPars(3,:)';
R_squared=FittingPars(4,:)';
Results=table(Concentration,Estimated_T1,LowerBound,UpperBound,R_squared);

cd(here)

%% plot T1 Error bar
figure()
y=FittingPars(1,:)';
ub=FittingPars(3,:)-FittingPars(1,:); ub=ub';
lb=FittingPars(1,:)-FittingPars(2,:); lb=lb';
 errorbar(Concentration,y,lb,ub,'-o','LineWidth',2,'MarkerSize',10);
 xlabel('Concentration','FontSize',20');
  ylabel('Estimated T1 (seconds)','FontSize',20');
  
function [FittingPars,Signal_hat] =  fit_variable_TR(RepetionTime,Signal)
% Estimated T1 relaxation time using variable TR data
%
%% Syntaxis
%
%   [FittingPars,Ypred,Yci,vtrModel] =  fit_T1_vtr(RepetionTime,Signal)
%
%% Inputs
%   RepetionTime= N x 1 vector of repetition time for VTR experiment. Units= seconds
%   Signal= N x 1 Matrix of signal for the VTR experiment. Units= AU
%
%% Ouputs
%   FittingPars= 5 X 1 numeric array with the following elements
%       FittingPars(1)= T1 time (estimated)
%       FittingPars(2)= standard eror of T1 time
%       FittingPars(3)= Lower 95% confience interval for the estimated T1
%       FittingPars(4)= Upper 95% confience interval for the estimated T1
%       FittingPars(5)= Adjusted Rsquared for the model 
%
%% Example
%
%   TR=0:.1:6;      
%   RepetionTime=TR';         
%   T1=2.7;
%   S=100 * ( 1 -exp (-TR./T1))';   Signal=awgn(S,25,'measured');
%   [FittingPars,Ypred,Yci,vtrModel] =  fitT1vtr(RepetionTime,Signal);
%  
%     plot(RepetionTime,Signal,'or',RepetionTime,Ypred,'-',RepetionTime,Yci,'--k');
%     xlabel('Repetition Time (sec)');    ylabel('Signal (au)');
%     title(['Estimated T1= ',num2str(FittingPars(1)), '   R-squared= ', num2str(FittingPars(end))])
%     legend({'Data','Predicted','Upper Bound','Lower Bound'},'Location','Best')
%
%%  Author
% Julio Cárdenas-Rodríguez
% University of Arizona
% Tucson, AZ
% cardenaslab.org


% Adjust orientation
p=size(RepetionTime);
if p(1) < p(2)
RepetionTime=RepetionTime';
end

p=size(Signal);
if p(1) < p(2)
Signal=Signal';
end

Signal=Signal./max(Signal);

% Initial Guess
x0=[1.1,2.0];
T1vTR_func=@(pars,xdata) pars(1) * ( 1 -exp (-xdata./pars(2)));
lb=[0.5,0.1]';
ub=[1.5,10]';

% Fit variable TR model
options = optimoptions('lsqcurvefit');
options.Display='off';
[Coefficients,~,residual,~,~,~,jacobian]=lsqcurvefit(T1vTR_func,x0,RepetionTime,Signal,lb,ub,options);

conf95 = nlparci(Coefficients,residual,'jacobian',jacobian);


% Allocate Variables

T1pred=Coefficients(2);
T1ci=conf95(2,:);
Signal_hat=T1vTR_func(Coefficients,RepetionTime);
Rsquared=rsquare(Signal,Signal_hat);

FittingPars=[T1pred,T1ci,Rsquared]'; 


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
 



 
 


