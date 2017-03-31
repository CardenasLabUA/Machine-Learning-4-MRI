clear all; clc; close all
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
%% Prepare data
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
xdata=ppm{1,1}(1,5:end)';
%% Normalize data
Xnormalized=X1;
%    for q=1:size(X1,1)
%    Xnormalized(q,:)=X1(q,:)  ./  trapz(xdata,Xnormalized(q,:)');
%    end
%% Components for PLS
Components=20;
%% PLS on pH
[~,~,XS_pH,~,beta,pctvarpH] = plsregress(Xnormalized,Y,Components);
yfit = [ones(size(Xnormalized,1),1) Xnormalized]*beta;
figure(100);
subplot(321); plot(cumsum(pctvarpH'),'o-'); legend('X','Y')
ylabel('% of variation');   xlabel('Component'); title('PLS on pH only'); ylim([0 1]);
subplot(322);
scatter(Y,yfit(:,1),100,C,'filled'); colorbar; colormap('jet'); lsline; 
xlabel('Measured pH');   ylabel('Predicted pH');  title('PLS on pH only'); 
%% PLS on Conc
[~,~,XS_Conc,~,beta,pctvarConc] = plsregress(Xnormalized,C,Components);
yfit = [ones(size(Xnormalized,1),1) Xnormalized]*beta;
figure(100);
subplot(323); plot(cumsum(pctvarConc'),'o-'); legend('X','Y')
ylabel('% of variation');   xlabel('Component'); title('PLS on Concentration Only'); ylim([0 1]);
subplot(324);
scatter(Y,yfit(:,1),100,C,'filled'); colorbar; colormap('jet'); lsline; 
xlabel('Measured pH');   ylabel('Predicted pH');title('PLS on Concentration Only');
%% PLS on Both
[~,~,XS_all,~,beta,pctvarall] = plsregress(Xnormalized,[Y,C],Components);
yfit = [ones(size(Xnormalized,1),1) Xnormalized]*beta;
figure(100);
subplot(325); plot(cumsum(pctvarall'),'o-'); legend('X','Y')
ylabel('% of variation');   xlabel('Component'); title('PLS on BOTH'); ylim([0 1]);
subplot(326);
scatter(Y,yfit(:,1),100,C,'filled'); colorbar; colormap('jet'); lsline; 
xlabel('Measured pH');   ylabel('Predicted pH');title('PLS on BOTH'); 
