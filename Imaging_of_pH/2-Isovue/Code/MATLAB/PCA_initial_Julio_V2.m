clear all; close all;
% pH6.25		pH6.77		pH7.19		pH8.00
% pH6.02		pH6.60		pH7.03		pH7.40

pHL={'pH6.02','pH6.60','pH6.77','pH7.03','pH7.19','pH7.40','pH8.00'};
 ppmexpected=linspace(-5,20,101);
 
[signal{1,1},ppm{1,1},B0{1,1}]=fit_CEST( '../data/pH6.02/CEST/' ); 
[signal{2,1},ppm{2,1},B0{2,1}]=fit_CEST( '../data/pH6.60/CEST/' ); 
[signal{3,1},ppm{3,1},B0{3,1}]=fit_CEST( '../data/pH6.77/CEST/' ); 
[signal{4,1},ppm{4,1},B0{4,1}]=fit_CEST( '../data/pH7.03/CEST/' ); 
[signal{5,1},ppm{5,1},B0{5,1}]=fit_CEST( '../data/pH7.19/CEST/' ); 
[signal{6,1},ppm{6,1},B0{6,1}]=fit_CEST( '../data/pH7.40/CEST/' ); 
[signal{7,1},ppm{7,1},B0{7,1}]=fit_CEST( '../data/pH8.00/CEST/' ); 

%%
figure(1);
for q=1:7
subplot(3,3,q);
plot(ppm{q,1}(:,5:end)',100 .* signal{q,1}(:,5:end)','.-'); xlim([ -5 20]);
title(pHL{q});
xlabel('Saturation Offset (ppm)'); 
ylabel('CEST contrast (%)');
ylim([0 101]);
end
%% Preapre for PCA
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
[X1, intervals, indexes] = icoshift (1-X1(50,:), 1-X1,'whole');

%% Do PCA on signal only
N=4;
[loadingX1,scoresX1,~,~,explX1,muX1]=pca(X1,'NumComponent',N);
figure(100);
subplot(131);
plot(cumsum(explX1),'xr'), hold all;
xlabel('PCA #'); ylabel('Variation explained');
subplot(132);

plot(loadingX1,'xr'), hold all;
xlabel('Offset'); ylabel('Amplitude');
subplot(133);

plot(muX1,'xr'),  hold all;
xlabel('Offset'); ylabel('CEST');
%% Do PCA on signal + B0 
[loadingX2,scoresX2,~,~,explX2,muX2]=pca([X1,B0_inhom],'NumComponent',N);
figure(100);
subplot(131);
plot(cumsum(explX2),'ob'), 
subplot(132);
plot(loadingX2,'ob'), 
subplot(133);
plot(muX2,'ob'),

%% PLS Z spectra only
N=10;
[~,~,~,~,~,PCTVAR1] = plsregress(zscore(X1),[Y],N);
figure(1);
subplot(3,3,8);
plot(1:N,cumsum(100*PCTVAR1(1,:)),'-o'); xlim([0 10]); ylim([0 101]);
xlabel('Number of component in Partial Least Square');
ylabel('% of Variance explained');
title('Dimensionality reduction of CEST spectra');
%%
figure(101);
plot(1:N,cumsum(100*PCTVAR1(1,:)),'-o'); hold all;
plot(1:N,cumsum(100*PCTVAR1(2,:)),'-x'); 
xlabel('Number of PLS components');
ylabel('% Explained in y');
 title('20 Components with and without B0');
legend('pH Method 1','Conc Method 1');
%%
Xpredictors=[X1];

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(zscore(Xpredictors),[Y,C],N);
figure(102);
plot(1:N,cumsum(100*PCTVAR(1,:)),'-o'); hold all;
plot(1:N,cumsum(100*PCTVAR(2,:)),'-x'); 

legend('pH Method 1','Conc Method 1','pH Method 2','Conc Method 2');
xlabel('Number of PLS components');
ylabel('% Explained in y');
legend('pH Method 2','Conc Method 2');
 title('Including Cocentration as a Predictor');

%%
figure(1);
subplot(3,3,9);
scatter(Y,stats.Yresiduals(:,1) + Y,100,C,'filled') ;
colorbar; colormap('jet'); xlabel('Measured'); ylabel('Predicted');


