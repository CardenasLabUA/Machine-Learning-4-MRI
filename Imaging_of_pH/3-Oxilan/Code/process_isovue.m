clear; close all;
data_folder = 'F:\GitHub\Machine-Learning-4-MRI\Imaging_of_pH\2-Isovue\Data\';
expected_cest_range = [-5,10];
[c_centered,c_xpred,pH,conc,ppm,c_rsq,c_signal]=fit_CEST_all_pH(data_folder,expected_cest_range);
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
end
subplot(3,3,9); plot(c_xpred,1-c_centered');title('All');

figure('Name','Raw Data');
cool=reshape(c_signal',size(ppm,2),9,8);
cool_ppm=reshape(ppm',size(ppm,2),9,8);
for q=1:8
    subplot(3,3,q); plot(cool_ppm(:,:,q),cool(:,:,q));
    title(['(Raw Data) pH = ',num2str(pHs(q))]);
    xlabel('Offset (ppm)')
    ylabel('Mz/Mo')
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
% 
% figure();
% hold on;
% concs = conc(1:9);
% 
% for j=1:8
%     plot(concs,1./T2(9*j-8:9*j));
%     plot(concs,1./T1(9*j-8:9*j));
%     pH_labels{j}=num2str(pH(9*j));
%     pH_labels{j+8}=num2str(pH(9*j));
% end
% legend(pH_labels); xlabel('Concentration (mM)'); ylabel('R1 or R2');


%% PLS Regression

Components=2:2:16;
figure();

for q=1:length(Components)
    [~,~,XS_pH,~,beta,pctvarpH] = plsregress(c_centered,pH,Components(q),'CV',size(c_centered,1));
    yfit = [ones(size(c_centered,1),1) c_centered]*beta;
    rmse_=(sum((pH-yfit(:,1)).^2)./ length(pH))^(1/2);
    subplot(4,2,q);
    scatter(pH,yfit(:,1),100,conc,'filled');
    h = colorbar;
    set(get(h,'title'),'string','mM'); colormap('jet'); lsline;
    xlabel('Measured pH');   ylabel('Predicted pH');
    title(['Components = ', num2str(Components(q)),'   rmse = ',num2str(rmse_)]);
end