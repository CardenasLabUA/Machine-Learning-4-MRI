clear all; close all;
% pH6.25		pH6.77		pH7.19		pH8.00
% pH6.02		pH6.60		pH7.03		pH7.40

pHL={'pH6.02','pH6.60','pH6.77','pH7.03','pH7.19','pH7.40','pH8.00'};
 ppmexpected=linspace(-5,20,101);
 
[signal{1,1},ppm{1,1},B0{1,1}]=fit_CEST( '../data/pH6.02/CEST/' ); 
[signal{2,1},ppm{2,1}]=fit_CEST( '../data/pH6.60/CEST/' ); 
[signal{3,1},ppm{3,1}]=fit_CEST( '../data/pH6.77/CEST/' ); 
[signal{4,1},ppm{4,1}]=fit_CEST( '../data/pH7.03/CEST/' ); 
[signal{5,1},ppm{5,1}]=fit_CEST( '../data/pH7.19/CEST/' ); 
[signal{6,1},ppm{6,1}]=fit_CEST( '../data/pH7.40/CEST/' ); 
[signal{7,1},ppm{7,1}]=fit_CEST( '../data/pH8.00/CEST/' ); 

%%
for q=1:7
subplot(4,2,q);
plot(ppm{q,1}(:,5:end)',signal{q,1}(:,5:end)','.-'); xlim([ -5 20]);
title(pHL{q});
end
%% Preapre for PCA
X1=(cell2mat(signal));  X1=X1(:,5:end);
X2=(cell2mat(ppm));     X2=X2(:,5:end);

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
%% Do PCA on signal only
[loadingX1,scoresX1,~,~,explX1,muX1]=pca(X1,'NumComponent',10);
figure(100);
subplot(131);
plot(cumsum(explX1),'-o'), title('X1');
xlabel('PCA #'); ylabel('Variation explained');
subplot(132);
plot(loadingX1,'-o'), title('Components X1 only');
xlabel('Offset'); ylabel('Amplitude');
subplot(133);
plot(muX1,'-o'), title('Components X1 only');
xlabel('Offset'); ylabel('CEST');
%% Do PCA on offset only
[loadingX2,scoresX2,~,~,explX2,muX2]=pca(X2,'NumComponent',2);
figure(100);
subplot(131);
plot(cumsum(explX2),'-x'), title('X2');
xlabel('PCA #'); ylabel('Variation explained');
subplot(132);
plot(loadingX2,'-x'), title('Components X2 only'); legend({'PCA-1','PCA-2'})
xlabel('Offset'); ylabel('Amplitude');
subplot(133);
plot(muX2,'-x'), title('Components X2 only');
xlabel('Offset'); ylabel('CEST');
%%
%% Do PCA on offset AND signal
[loadingX3,scoresX3,~,~,explBoth,muBoth]=pca([X1,X2],'NumComponent',2);
figure(100);
subplot(131);
plot(cumsum(explBoth),'-x'), title('Both');
xlabel('PCA #'); ylabel('Variation explained');
subplot(132);
plot(loadingX3,'-x'), title('Components Both'); legend({'PCA-1','PCA-2'})
xlabel('Offset'); ylabel('Amplitude');
subplot(133);
plot(muBoth,'-x'), title('Components Both');
xlabel('Offset'); ylabel('CEST');
%% PLS
for q=1:20
N=q;
[XL,yl,XS,YS,beta,PCTVAR] = plsregress(zscore([X1,X2]),[Y,C],N);
subplot(4,5,q);
plot(1:N,cumsum(100*PCTVAR(1,:)),'-o'); hold all;
plot(1:N,cumsum(100*PCTVAR(2,:)),'-x'); 
legend('pH','Conc');
%xlabel('Number of PLS components');
%ylabel('%e Explained in y');
 title(num2str(q));
end
%%



