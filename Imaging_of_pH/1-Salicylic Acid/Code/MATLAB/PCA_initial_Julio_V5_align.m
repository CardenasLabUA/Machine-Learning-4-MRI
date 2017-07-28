clear all; 
% pH6.25		pH6.77		pH7.19		pH8.00
% pH6.02		pH6.60		pH7.03		pH7.40

pHL={'pH= 6.02','pH= 6.60','pH= 6.77','pH =7.03','pH= 7.19','pH= 7.40','pH= 8.00'};
 ppmexpected=linspace(-5,20,101);
 
[signal{1,1},ppm{1,1},B0{1,1}]=fit_CEST( '../data/pH6.02/CEST/' ); 
[signal{2,1},ppm{2,1},B0{2,1}]=fit_CEST( '../data/pH6.60/CEST/' ); 
[signal{3,1},ppm{3,1},B0{3,1}]=fit_CEST( '../data/pH6.77/CEST/' ); 
[signal{4,1},ppm{4,1},B0{4,1}]=fit_CEST( '../data/pH7.03/CEST/' ); 
[signal{5,1},ppm{5,1},B0{5,1}]=fit_CEST( '../data/pH7.19/CEST/' ); 
[signal{6,1},ppm{6,1},B0{6,1}]=fit_CEST( '../data/pH7.40/CEST/' ); 
[signal{7,1},ppm{7,1},B0{7,1}]=fit_CEST( '../data/pH8.00/CEST/' ); 
%% Prepare for PCA
X1=(cell2mat(signal));  X1=X1(:,5:end);
X2=(cell2mat(ppm));     X2=X2(:,5:end);
B0_inhom=cell2mat(B0);
% Concentration
Conc(1,1)=13.42;
for q=2:9
  Conc(q,1)  =Conc(q-1,1) .* (1/0.80);
end
C=Conc * ones(1,7);
C=C(:);
% pH
pH=[6.02,6.60,6.77,7.03,7.19,7.40,8.00]';
Y=pH * ones(1,9);
Y=Y(:);
%% Align data
[X1, intervals, indexes] = icoshift (1-X1(1,:), 1-X1,'whole');
% %% Normalize
%  Xnormalized=X1;
%   for q=1:size(X1,1)
%   Xnormalized(q,:)=X1(q,:)./norm(X1(q,:));
%   end
%%

xdata=ppm{1,1}(1,5:end)';

 Xnormalized=X1;
%   for q=1:size(X1,1)
%   Xnormalized(q,:)=X1(q,:)./norm(X1(q,:));
%   end

   for q=1:size(X1,1)
   Xnormalized(q,:)=X1(q,:)  ./  trapz(xdata,Xnormalized(q,:)');
   end

%% PLS Z spectra only
N=5;
 RESP=[Y,C./10];
 Increment=7;
for k=1:8
Ydata=     RESP(k:Increment:end,:);
    Xdata= Xnormalized(k:Increment:end,:);
    
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(Xdata,Ydata,N);
% Predict Reponse
yfit = [ones(size(Xdata,1),1) Xdata]*BETA;


subplot(2,4,k);
scatter(Ydata(:,2),yfit(:,2),100,Ydata(:,1),'filled') ; hold all;
colorbar; colormap('jet'); xlabel('Measured conc.'); ylabel('Predicted conc');
lsline
end
%% LASSO Concentration from PLS Conc
Cscaled=C./10;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(Xnormalized,RESP,10);
[~,FitInfo] = lasso(XS,Cscaled,'CV',10);
i=FitInfo.IndexMinMSE;
[B,FitInfo] = lasso(XS,Cscaled,'CV',10,'Lambda',FitInfo.Lambda(i));
Bhat=[B;FitInfo.Intercept];
Cscaledhat=[XS,ones(size(XS,1),1)] * Bhat;
figure() ; scatter(Cscaled,Cscaledhat,100,Y,'filled')

%% LASSO pH from PLS pH
[~,~,XS] = plsregress(Xnormalized,Y,10);
[~,FitInfo] = lasso(XS,Y,'CV',10);
i=FitInfo.IndexMinMSE;
[B,FitInfo] = lasso(XS,Y,'CV',10,'Lambda',FitInfo.Lambda(i));
Bhat=[B;FitInfo.Intercept];
Yhat=[XS,ones(size(XS,1),1)] * Bhat;
figure() ; scatter(Y,Yhat,100,C,'filled')
%%
%% LASSO pH from PLS Both
[~,~,XS] = plsregress(Xnormalized,RESP,10);
[~,FitInfo] = lasso(XS,Y,'CV',10);
i=FitInfo.IndexMinMSE;
[B,FitInfo] = lasso(XS,Y,'CV',10,'Lambda',FitInfo.Lambda(i));
Bhat=[B;FitInfo.Intercept];
Yhat=[XS,ones(size(XS,1),1)] * Bhat;
figure() ; scatter(Y,Yhat,100,C,'filled')
%%
%% LASSO pH from PCA
[~,PCAscores] = pca(Xnormalized,'NumComponents',10);
[~,FitInfo] = lasso(PCAscores,Y,'CV',10);
i=FitInfo.IndexMinMSE;
[B,FitInfo] = lasso(PCAscores,Y,'CV',10,'Lambda',FitInfo.Lambda(i));
Bhat=[B;FitInfo.Intercept];
Yhat=[PCAscores,ones(size(PCAscores,1),1)] * Bhat;
figure() ; scatter(Y,Yhat,100,C,'filled')

%%
myX=rand(65,101);
myB(1,randperm(63,10))=randn(1,10)';
mY=myX * myB';
[B,FitInfo] = lasso(myX,mY,'CV',10); %lassoPlot(B,FitInfo,'PlotType','CV');
i=FitInfo.IndexMinMSE;
[B,FitInfo] = lasso(myX,mY,'CV',10,'Lambda',FitInfo.Lambda(i));
myBhat=[B;FitInfo.Intercept];
mYhat=[myX,ones(size(myX,1),1)] * myBhat;
figure();
plot(mY); hold all; plot(mYhat,'s'); hold all; 
%%
rng default
Mdl = fitrsvm(XS,Y,'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus'))
