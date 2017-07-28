clear; close all; clc;
folder_path = 'C:\Users\Luis\Desktop\Oxilan\Averages\PBS\CEST\';
% Set up different pHs
here=pwd;
center_range = [-8,15];
cd(folder_path);
files=dir('pH*.txt');
load('all_voxel_centered_spectra_with_T2.mat');
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
    method.Npools=3;

    % Prescribe initial guesses and lower and upper bounds for amplitudes,
    % widths, and offsets of lorentzian peaks from the 3 pools
    method.x0=[0.4, 1, 1.5;... Water amplitude
    0, 0, 0.7;... Oxilan pool 1 amplitude
    0, 0, 1;... Oxilan hydroxyls amplitude
    0, 4, 10;... Water width
    0, 1, 15;... Oxilan pool 1 width
    0, 1, 5;... Oxilan hydroxyls width
    -0.5, 0, 0.5;... Water offset
    3.8, 4.1, 5.2;... Oxilan pool 1 offset
    0.25, 0.7, 3];... Oxilan hydroxyls offset

    control=find(round(ppm_pbs)==center_range(1),1);
    if isempty(control)
        control=1;
    end
    
    method.range=[1,length(s)];
    fit=cf_Lorentzian(s,ppm_pbs,method,s(control));
    fit.pars(7:9)=fit.pars(7:9)-fit.pars(7);
    c_pbs_fit(p,:)=lorentzian(fit.pars,xpred)';
    c_pbs_unfit(p,:)=s./s(control);
    c_pbs_rsq(p)=fit.rsq;
    c_pbs_ppm(p,:)=ppm_pbs;
    
    c_dif(c_pH==pHs(p),:)=c_dif(c_pH==pHs(p),:)-repmat(c_pbs_fit(p,:),size(c_dif(c_pH==pHs(p),:),1),1);
end

cd(here);


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
