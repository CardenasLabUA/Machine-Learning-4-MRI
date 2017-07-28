clear; close all; clc;
folder_path = 'D:\Users\Joey D\Desktop\Joey\GitHub\Machine-Learning-4-MRI\Imaging_of_pH\3-Oxilan\Code\ML_python\SA\Images\';
% Set up different pHs
here=pwd;
cd(folder_path)
pH_folders=dir('pH*');
n_pH=length(pH_folders);
conc=[13.42;16.77;20.97;26.21;32.76;40.96;51.20;64.00;80.00];
load('mask_SA.mat');

xpred=linspace(-5,20,101);
te=linspace(10,500,50)/1000;
% Extract CEST data
c_centered=[];
c_uncentered=[];
c_ppm=[];
c_rsq=[];
c_conc=[];
c_pH=[];
c_pars=[];

%% Method for Lorentzian fitting
% Prescribe number of pools for lorentzian fit
method.Npools=2;

% Prescribe initial guesses and lower and upper bounds for amplitudes,
        % widths, and offsets of lorentzian peaks from the 3 pools
        method.x0=[0.4, .7, 1.0;... Water amplitude - 1
        0, 0.3, 1;... SA pool 1 amplitude - 2
        0.5, 4, 10;... Water width - 3
        .5, 4, 10;... SA pool 1 width - 4
        -1, 0, 1;... Water offset - 5
        8, 10, 11.00];... SA pool 1 offset - 6

for p=1:n_pH
    cest_img = reshape(Bruker_reader([folder_path,pH_folders(p).name,'\1']),128,128,105);
    pHs(p)=str2num(pH_folders(p).name(3:end));
    mask=squeeze(mask_SA(:,:,p));
    
    %% Determine concentration of each voxel
    scan_h_mask=zeros(128,128);
    scan_v_mask=scan_h_mask;
    concs=zeros(3,3);
    conc_sort=sort(conc(1:9),'descend');
    concs(1,:)=conc_sort(1:3);
    concs(2,:)=conc_sort(4:6);
    concs(3,:)=conc_sort(7:end);
    conc_dims=zeros(1,6);

    for j=1:128
        h_num=0;
        h_count=0;
        for k=1:128
            if mask(j,k)==1
                scan_h_mask(j,k)=h_num;
                h_count=0;
            else
                h_count=h_count+1;
            end
            if h_count==6
                h_num=h_num+1;
            end
        end
    end
    row=find(sum(scan_h_mask~=0,2)==max(sum(scan_h_mask~=0,2)),1);
    conc_dims(1,2)=find(scan_h_mask(row,:)==1,1,'last')+5;
    conc_dims(1,4)=find(scan_h_mask(row,:)==2,1,'last')+5;
    conc_dims(1,6)=find(scan_h_mask(row,:)==3,1,'last')+5;

    for j=1:128
        v_num=0;
        v_count=0;
        for k=1:128
            if mask(k,j)==1
                scan_v_mask(k,j)=v_num;
                v_count=0;
            else
                v_count=v_count+1;
            end
            if v_count==6
                v_num=v_num+1;
            end
        end
    end
    col=find(sum(scan_v_mask~=0,1)==max(sum(scan_v_mask~=0,1)),1);
    conc_dims(1,1)=find(scan_v_mask(:,col)==1,1,'last')+5;
    conc_dims(1,3)=find(scan_v_mask(:,col)==2,1,'last')+5;
    conc_dims(1,5)=find(scan_v_mask(:,col)==3,1,'last')+5;

    conc_mask=zeros(128,128);
    conc_mask(1:conc_dims(1),1:conc_dims(2))=concs(1,1);
    conc_mask(conc_dims(1):conc_dims(3),1:conc_dims(2))=concs(2,1);
    conc_mask(conc_dims(3):conc_dims(5),1:conc_dims(2))=concs(3,1);
    conc_mask(1:conc_dims(1),conc_dims(2):conc_dims(4))=concs(1,2);
    conc_mask(conc_dims(1):conc_dims(3),conc_dims(2):conc_dims(4))=concs(2,2);
    conc_mask(conc_dims(3):conc_dims(5),conc_dims(2):conc_dims(4))=concs(3,2);
    conc_mask(1:conc_dims(1),conc_dims(4):conc_dims(6))=concs(1,3);
    conc_mask(conc_dims(1):conc_dims(3),conc_dims(4):conc_dims(6))=concs(2,3);
    conc_mask(conc_dims(3):conc_dims(5),conc_dims(4):conc_dims(6))=concs(3,3);
    
    mask(conc_mask==0)=0;
    
    %% Organize all CEST and T2 spectra by pH and concentration
    indices=find(reshape(mask,[],1));
    cest_img=reshape(cest_img,[],105);
    centered=nan(length(indices),101);
    uncentered=centered;
    ppms=centered;
    d_pars=nan(length(indices),6);
    rsq=nan(length(indices),1);
    pH=rsq;
    d_conc=rsq;
    
    for q=1:length(indices)
        s=cest_img(indices(q),:);
        s=s(5:end);
        [~,i]=max(1-s);
        ppm=xpred-xpred(i);

        control=1;
        method.range=[1,length(s)];
        fit=cf_Lorentzian(s,ppm,method,s(control));
        ppm_cent=ppm-fit.pars(5);
        fit=cf_Lorentzian(s,ppm_cent,method,s(control));
        fit.pars(5:6)=fit.pars(5:6)-fit.pars(5);
        uncentered(q,:)=s./s(control);
        ppms(q,:)=ppm_cent;
        centered(q,:)=lorentzian(fit.pars,xpred)';
        rsq(q,1)=fit.rsq;
        pH(q,1)=pHs(p);
        d_conc(q,1)=conc_mask(indices(q));
        d_pars(q,1)=fit.pars(2);
        d_pars(q,2)=fit.pars(4);
        d_pars(q,3)=fit.pars(6);
        d_pars(q,4)=fit.pars(1);
        d_pars(q,5)=fit.pars(3);
        d_pars(q,6)=fit.pars(5);
        
    
   
    end
    
    %%
    c_uncentered=[c_uncentered;uncentered];
    c_ppm=[c_ppm;ppms];
    c_centered=[c_centered;centered];
    c_rsq=[c_rsq;rsq];
    c_pH=[c_pH;pH];
    c_conc=[c_conc;d_conc];
    c_pars=[c_pars;d_pars];      
    cd(here);
    save('all_voxel_centered_spectra_no_T2_SA.mat','xpred','c_centered','c_rsq','c_pH','c_conc','c_pars',...
        'c_ppm','c_uncentered','pHs','concs','-v7');
    cd(folder_path);
    p
end

cd(here);


%%
folder_path = 'D:\Users\Joey D\Desktop\Joey\GitHub\Machine-Learning-4-MRI\Imaging_of_pH\1-Salicylic Acid\Data\PBS_only\CEST\';
% Set up different pHs
here=pwd;
center_range = [-5,20];
cd(folder_path);
files=dir('pH*.txt');
load('all_voxel_centered_spectra_no_T2_SA.mat');
c_dif = c_centered;

for p=1:length(pHs)
    offset=importdata('ppm.txt');
    offset=offset(5:end)';
    
    s= importdata(files(p).name)'; 
    s=s(5:end);
    [~,i]=max(1-s);
    ppm_pbs=offset-offset(i);
    
    %% Center data via Lorentzian fitting
    % Prescribe number of pools for lorentzian fit
    method.Npools=1;

    % Prescribe initial guesses and lower and upper bounds for amplitudes,
    % widths, and offsets of lorentzian peaks from the 3 pools
    method.x0=[0.4, 1, 1.5;... Water amplitude
    0, 4, 10;... Water width
    -2, 0, 2];... Water offset

    control=1;
    
    method.range=[1,length(s)];
    fit=cf_Lorentzian(s,ppm_pbs,method,s(control));
    fit.pars(3)=0;
    c_pbs_fit(p,:)=lorentzian(fit.pars,xpred)';
    c_pbs_unfit(p,:)=s./s(control);
    c_pbs_rsq(p)=fit.rsq;
    c_pbs_ppm(p,:)=ppm_pbs;
    
    c_dif(c_pH==pHs(p),:)=c_dif(c_pH==pHs(p),:)-repmat(c_pbs_fit(p,:),size(c_dif(c_pH==pHs(p),:),1),1);
end

cd(here);
save('all_voxel_centered_spectra_no_T2_SA.mat','xpred','c_centered','c_rsq','c_pH','c_conc','c_pars',...
        'c_ppm','c_uncentered','pHs','concs','c_dif','-v7');


%%
figure();
for p=1:8
    subplot(4,2,p); hold on;
    c_plot=c_dif(c_pH==pHs(p),:);
    for j=1:size(c_plot,1)
        plot(xpred,1-c_plot(j,:)); ylim([0,1]); xlim([-8,15]);
    end
    hold off;
end

%%
figure();
for p=1:8
    subplot(4,2,p); hold on;
    c_plot=c_pbs_fit(p,:);
    plot(xpred,1-c_plot); ylim([0,1]); xlim([-8,15]);
    c_plot=c_pbs_unfit(p,:);
    plot(c_pbs_ppm(p,:),c_plot,'o'); ylim([0,1]); xlim([-8,15]);
    hold off;
end

%% plot uncentered cest
figure();
rsq_cut=-1;
for p=1:8
    subplot(4,2,p); hold on;
    c_plot=c_uncentered(c_rsq>rsq_cut,:);
    c_plot=c_plot(c_pH(c_rsq>rsq_cut)==pHs(p),:);
    c_plot_ppm=c_ppm(c_rsq>rsq_cut,:);
    c_plot_ppm=c_plot_ppm(c_pH(c_rsq>rsq_cut)==pHs(p),:);
    for j=1:size(c_plot,1)
        plot(c_plot_ppm(j,:),c_plot(j,:)); ylim([0,1]);
    end
    hold off;
end

%% plot centered cest
rsq_cut=0.90;
figure();
for p=1:8
    subplot(4,2,p); hold on;
    c_plot=c_centered(c_rsq>rsq_cut,:);
    c_plot=c_plot(c_pH(c_rsq>rsq_cut,:)==pHs(p),:);
    for j=1:size(c_plot,1)
        plot(xpred,1-c_plot(j,:)); ylim([0,1]); xlim([-8,15]);
    end
    hold off;
end


%% Trainer ########################################################
load('all_voxel_centered_spectra_no_T2_SA.mat');
rsq_cut = -1;
training_set=c_centered(c_rsq>rsq_cut,:);
Components=2:1:size(training_set,2);
figure();
rmse_=[];
response=c_pH(c_rsq>rsq_cut);
other=c_conc(c_rsq>rsq_cut);

for q=1:length(Components)
    [~,~,XS_pH,~,beta,pctvarpH] = plsregress(training_set,response,Components(q),'cv','resubstitution');
    yfit = [ones(size(training_set,1),1) training_set]*beta;
    rmse_(q)=(sum((response-yfit(:,1)).^2)./ length(response))^(1/2) ./ mean(response)*100;
end
[~,i]=min(rmse_);

plot(Components,rmse_); xlabel('Components'); ylabel('NRMSE for Prediction (%)');

rmse_=[];
Components=Components(i);
[~,~,XS_pH,~,beta,pctvarpH] = plsregress(training_set,response,Components,'CV',round(size(training_set,1)/10));
yfit = [ones(size(training_set,1),1) training_set]*beta;
rmse_=(sum((response-yfit(:,1)).^2)./ length(response))^(1/2) ./ mean(response)*100;

figure();
scatter(response,yfit(:,1),100,other,'filled');
h = colorbar; colormap('jet'); lsline;
xlabel('Measured','FontSize',18);   

ylabel('Predicted','FontSize',18);
title(['Components = ',num2str(Components),'   NRMSE = ',num2str(rmse_),'%'],'FontSize',18);

for p=1:length(pHs)
    rmse_pH(p)=(sum((response(response==pHs(p))-yfit(response==pHs(p),1)).^2)./ length(response(response==pHs(p))))^(1/2) ./ mean(response(response==pHs(p)))*100;
end

figure();

plot(pHs,rmse_pH);

%% Compare to avg
data_folder = 'D:\Users\Joey D\Desktop\Joey\GitHub\Machine-Learning-4-MRI\Imaging_of_pH\1-Salicylic Acid\Data\';
expected_cest_range = [-5,20];
[a_centered,xpred,a_pH,a_conc,a_ppm,a_rsq,a_signal]=fit_CEST_all_pH_SA(data_folder,expected_cest_range);
% training_set=a_centered;
% response=a_pH;
% other=a_conc;

training_set_a=a_centered;
Components=2:1:71;
rmse_=[];
response_a=a_pH;

figure();
for q=1:length(Components)
    [~,~,XS_pH,~,beta,pctvarpH] = plsregress(training_set_a,response_a,Components(q),'cv','resubstitution');
    yfit = [ones(size(training_set,1),1) training_set]*beta;
    rmse_(q)=(sum((response-yfit(:,1)).^2)./ length(response))^(1/2) ./ mean(response)*100;
%     rsq_(q)=rsq(a_pH,yfit);
end
[~,i]=min(rmse_);

plot(Components,rmse_); xlabel('Components'); ylabel('NRMSE for Prediction (%)');

rmse_=[];
Components=Components(i);
[~,~,XS_pH,~,beta,pctvarpH] = plsregress(training_set_a,response_a,Components,'cv','resubstitution');
yfit = [ones(size(training_set,1),1) training_set]*beta;
rmse_=(sum((response-yfit(:,1)).^2)./ length(response))^(1/2) ./ mean(response)*100;

figure();
scatter(response,yfit(:,1),100,other,'filled');
h = colorbar; colormap('jet'); lsline;
xlabel('Measured','FontSize',18);   

ylabel('Predicted','FontSize',18);
title(['Components = ',num2str(Components),'   NRMSE = ',num2str(rmse_),'%'],'FontSize',18);

for p=1:length(pHs)
    rmse_pH(p)=(sum((response(response==pHs(p))-yfit(response==pHs(p),1)).^2)./ length(response(response==pHs(p))))^(1/2) ./ mean(response(response==pHs(p)))*100;
end

figure();

plot(pHs,rmse_pH);