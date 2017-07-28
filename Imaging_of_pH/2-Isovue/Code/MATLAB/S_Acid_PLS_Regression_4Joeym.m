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
C=Conc * ones(1,7); C=C(:);
% pH
pH=[6.02,6.60,6.77,7.03,7.19,7.40,8.00]';
Y=pH * ones(1,9);   Y=Y(:);
%% Align data
[X1, intervals, indexes] = icoshift (1-X1(1,:), 1-X1,'whole');
xdata=ppm{1,1}(1,5:end)';

%% Plot
cool=reshape(X1',length(xdata),9,7);
for q=1:7
    subplot(2,4,q); plot(xdata,1-cool(:,:,q));
    title(pHL{q}); ylim([0 1]);
    xlabel('Offset (ppm)')
    ylabel('Mz/Mo')
end
subplot(2,4,8); plot(xdata,X1');title('All');
%% Normalize data
Xnormalized=X1;
%      for q=1:size(X1,1)
%      Xnormalized(q,:)=X1(q,:)  ./  trapz(xdata,Xnormalized(q,:)');
%      end

%% PLS on pH
Components=2:2:16;

for q=1:length(Components)
    [~,~,XS_pH,~,beta,pctvarpH] = plsregress(Xnormalized,Y,Components(q),'CV',63);
    yfit = [ones(size(Xnormalized,1),1) Xnormalized]*beta;
    rmse_=(sum((Y-yfit(:,1)).^2)./ length(Y))^(1/2);
    figure(2);
    subplot(4,2,q);
    scatter(Y,yfit(:,1),100,C,'filled');
    h = colorbar;
    set(get(h,'title'),'string','mM'); colormap('jet'); lsline;
    xlabel('Measured pH');   ylabel('Predicted pH');
    title(['Components = ', num2str(Components(q)),'   rmse = ',num2str(rmse_)]);
end
%%