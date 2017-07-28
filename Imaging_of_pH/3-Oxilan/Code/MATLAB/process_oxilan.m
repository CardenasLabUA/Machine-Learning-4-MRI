clear; close all;
data_folder = 'D:\Users\Joey D\Desktop\Joey\GitHub\Machine-Learning-4-MRI\Imaging_of_pH\3-Oxilan\Data\';
expected_cest_range = [-8,15];
[c_centered,c_xpred,pH,conc,ppm,c_rsq,c_signal,~,offset]=fit_CEST_all_pH(data_folder,expected_cest_range);
% [T1,T1_error,T1_rsq,T1_signal,TR]=fit_T1_all_pH(data_folder);
% [T2,T2_error,T2_rsq,T2_signal,TE]=fit_T2_all_pH(data_folder);

%% Plotting
figure('Name','Centered Data');

pHs = pH(1:9:end);
cool=reshape(c_centered',length(c_xpred),9,8);
for q=1:8
    subplot(3,3,q); plot(c_xpred,1-cool(:,:,q));
    title(['(Centered Data) pH = ',num2str(pHs(q))]); ylim([0 1]);
    xlabel('Offset (ppm)')
    ylabel('Mz/Mo')
    xlim(expected_cest_range);
end
subplot(3,3,9); plot(c_xpred,1-c_centered');title('All');

figure('Name','Raw Data');
cool=reshape(c_signal',size(ppm,2),9,8);
cool_c=reshape(c_centered',length(c_xpred),9,8);
cool_ppm=reshape(offset',size(ppm,2),9,8);
for q=1:8
    subplot(3,3,q); plot(cool_ppm(:,:,q),cool(:,:,q),'o'); hold on; plot(c_xpred,1-cool_c(:,:,q),'-'); hold off;
    title(['(Raw Data) pH = ',num2str(pHs(q))]); ylim([0,1]);
    xlabel('Offset (ppm)')
    ylabel('Mz/Mo')
    xlim(expected_cest_range);
end
subplot(3,3,9); plot(ppm',c_signal');title('All');
    
% figure();
% hold on;
% 
% for j=1:9
%     plot(pHs,1./T2(j:9:end));
%     plot(pHs,1./T1(j:9:end));
%     conc_labels{j}=num2str(conc(j));
%     conc_labels{j+9}=num2str(conc(j));
% end
% legend(conc_labels); xlabel('pH'); ylabel('R1 or R2');

% %%
% figure();
% hold on;
% concs = conc(1:9);
% 
% for j=1:8
%     plot(concs,1./T2(9*j-8:9*j));
%     ob=fitlm(concs,1./T2(9*j-8:9*j));
%     vars=ob.Coefficients.Variables;
%     r2(j)=vars(2,1);
%     plot(concs,1./T1(9*j-8:9*j));
%     ob=fitlm(concs,1./T1(9*j-8:9*j));
%     vars=ob.Coefficients.Variables;
%     r1(j)=vars(2,1);
%     pH_labels{j}=num2str(pH(9*j));
%     pH_labels{j+8}=num2str(pH(9*j));
% end
% legend(pH_labels); xlabel('Concentration (mM)'); ylabel('R1 or R2');
% 
% figure();
% plot(pHs,r2); xlabel('pH'); ylabel('r2');
% 
% figure();
% plot(pHs,r1); xlabel('pH'); ylabel('r1');
% 
% 
% %% PLS Regression
% 
% Components=[2,52,60];
% figure();
% 
% for q=1:length(Components)
%     [~,~,XS_pH,~,beta,pctvarpH] = plsregress(c_centered,pH,Components(q),'cv',size(c_centered,1));
%     yfit = [ones(size(c_centered,1),1) c_centered]*beta;
%     rmse_=(sum((pH-yfit(:,1)).^2)./ length(pH))^(1/2) ./ mean(pH) *100;
%     subplot(2,2,q,'FontSize',18);
%     scatter(pH,yfit(:,1),100,conc,'filled');
%     h = colorbar;
%     set(get(h,'title'),'string','mM'); colormap('jet'); lsline;
%     xlabel('Measured pH','FontSize',18);   
%     
%     ylabel('Predicted pH','FontSize',18);
%     title(['Components = ', num2str(Components(q)),'   rmse = ',num2str(rmse_)],'FontSize',18);
% end
% %%
% Components=2:8:64;
% figure();
% 
% for q=1:length(Components)
%     [~,~,XS_conc,~,beta,pctvarconc] = plsregress(c_centered,conc,Components(q),'cv',size(c_centered,1));
%     yfit = [ones(size(c_centered,1),1) c_centered]*beta;
%     rmse_=(sum((conc-yfit(:,1)).^2)./ length(conc))^(1/2);
%     subplot(4,2,q);
%     scatter(conc,yfit(:,1),100,pH,'filled');
%     h = colorbar;
%     set(get(h,'title'),'string','pH'); colormap('jet'); lsline;
%     xlabel('True Concentration (mM)');   ylabel('Predicted Concentration (mM)');
%     title(['Components = ', num2str(Components(q)),'   rmse = ',num2str(rmse_)]);
% end
% 
% 
% 
% %% rmse vs components
% Components=2:1:71;
% figure();
% 
% for q=1:length(Components)
%     [~,~,XS_pH,~,beta,pctvarpH] = plsregress(c_centered,pH,Components(q),'cv',size(c_centered,1));
%     yfit = [ones(size(c_centered,1),1) c_centered]*beta;
%     rmse_(q)=(sum((pH-yfit(:,1)).^2)./ length(pH))^(1/2);
% end
% 
% plot(Components,rmse_); xlabel('Components'); ylabel('RMSE for pH Prediction'); ylim([0,0.35]); xlim([1,71]);
% 
% %% rmse vs components no CV
% Components=2:1:71;
% figure();
% 
% for q=1:length(Components)
%     [~,~,XS_pH,~,beta,pctvarpH] = plsregress(c_centered,pH,Components(q));
%     yfit = [ones(size(c_centered,1),1) c_centered]*beta;
%     rmse_(q)=(sum((pH-yfit(:,1)).^2)./ length(pH))^(1/2);
% end
% 
% plot(Components,rmse_); xlabel('Components'); ylabel('RMSE for pH Prediction'); ylim([0,0.35]); xlim([1,71]);
% 
% %% lower res
% figure(98);
% hold on;
% figure(99);
% hold on;
% for k=1:4
% zspec_low = c_centered(:,1:k:end);
% figure(99);
% plot(c_xpred(1:k:end),zspec_low(72,:),'o-');
% if size(zspec_low,2)<=71
%     Components=2:1:size(zspec_low,2);
% else
%     Components=2:1:71;
% end
% rmse_=[];
% 
% for q=1:length(Components)
%     [~,~,XS_pH,~,beta,pctvarpH] = plsregress(zspec_low,pH,Components(q),'cv',size(zspec_low,1));
%     yfit = [ones(size(zspec_low,1),1) zspec_low]*beta;
%     rmse_(q)=(sum((pH-yfit(:,1)).^2)./ length(pH))^(1/2);
% end
% figure(98);
% plot(Components,rmse_); xlabel('Components'); ylabel('rmse');
% end
% 
% xlim([2,25]); ylim([0,1]);