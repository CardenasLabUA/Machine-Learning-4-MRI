function [centered,xpred,pH,conc,ppm,c_rsq,signal,B0_error,ppm_cent]=fit_CEST_all_pH( folder_path, center_range, method )
% Set up different pHs
here=pwd;
cd(folder_path)
pH_folders=dir('pH*');
n_pH=length(pH_folders);

% Extract CEST data
last_n=0;
for p=1:n_pH
    cd([folder_path,pH_folders(p).name,'\CEST']);
    offset=importdata('ppm.txt');
    offset=offset(5:end);
    conc_files=dir('*mM*.txt');
    n_conc=length(conc_files);

     for q=1:n_conc
        pH(q+last_n,1)=str2num(pH_folders(p).name(3:end));
        if isempty(str2num(conc_files(q).name(1:6)))
            conc(q+last_n,1)=str2num(conc_files(q).name(1:5));
        else
            conc(q+last_n,1)=str2num(conc_files(q).name(1:6));
        end
        s= importdata(conc_files(q).name); 
        s=s(5:end);
        [~,i]=max(1-s);
        B0_error(q,1)=offset(i);
        ppm(q+last_n,:)=offset-offset(i);
        control=1;
        s=s./s(control);
        signal(q+last_n,:)=s;
     end
     last_n=last_n+n_conc;
end

%% Center data via Lorentzian fitting
% Prescribe number of pools for lorentzian fit
method.Npools=3;

% Prescribe initial guesses and lower and upper bounds for amplitudes,
% widths, and offsets of lorentzian peaks from the 3 pools
method.x0=[0.4, 1, 1.0;... Water amplitude
    0, 0.2, 1;... Oxilan pool 1 amplitude
    0, 0.1, 1;... Oxilan hydroxyls amplitude
    .2, 4, 15;... Water width
    0, 1, 15;... Oxilan pool 1 width
    0, 1, 15;... Oxilan hydroxyls width
    -0.5, 0, 0.5;... Water offset
    3.8, 4.1, 5.2;... Oxilan pool 1 offset
    0.25, 0.7, 2.5];... Oxilan hydroxyls offset

xpred=linspace(center_range(1),center_range(2),length(offset));
for j=1:size(signal,1)
    
    method.range=[1,length(signal(j,:))];
    fit=cf_Lorentzian(signal(j,:),ppm(j,:),method,signal(j,control));
    ppm_cent(j,:)=ppm(j,:)-fit.pars(7);
    fit=cf_Lorentzian(signal(j,:),ppm_cent(j,:),method,signal(j,control));
    fit.pars(7:9)=fit.pars(7:9)-fit.pars(7);
    centered(j,:)=lorentzian(fit.pars,xpred)';
    c_rsq(1,j)=fit.rsq;
end
    
 cd(here)  
end
 



 
 


