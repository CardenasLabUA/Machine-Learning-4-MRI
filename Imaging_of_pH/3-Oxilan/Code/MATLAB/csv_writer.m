clear;
load('all_voxel_centered_spectra_with_T2.mat');
load('mask_oxi.mat');
here=pwd;
cd('D:\Users\Joey D\Desktop\Joey\GitHub\Machine-Learning-4-MRI\Imaging_of_pH\3-Oxilan\Code\ML_python\Oxilan');
csvwrite('fitted_cest.csv',c_centered);
csvwrite('unfitted_cest.csv',c_uncentered);
csvwrite('ppm.csv',c_ppm);
csvwrite('pars.csv',c_pars);
csvwrite('rsq.csv',c_rsq);
csvwrite('pH.csv',c_pH);
csvwrite('pHs.csv',pHs);
csvwrite('conc.csv',c_conc);
csvwrite('concs.csv',reshape(concs,9,1));
csvwrite('diff.csv',c_dif);
csvwrite('mask.csv',mask_oxi);
cd(here);

%% Isovue Voxel
clear;
load('all_voxel_centered_spectra_no_T2_isovue.mat');
load('mask_iso.mat');
here=pwd;
cd('D:\Users\Joey D\Desktop\Joey\GitHub\Machine-Learning-4-MRI\Imaging_of_pH\3-Oxilan\Code\ML_python\Isovue');
csvwrite('fitted_cest.csv',c_centered);
csvwrite('unfitted_cest.csv',c_uncentered);
csvwrite('ppm.csv',c_ppm);
csvwrite('pars.csv',c_pars);
csvwrite('rsq.csv',c_rsq);
csvwrite('pH.csv',c_pH);
csvwrite('pHs.csv',pHs);
csvwrite('conc.csv',c_conc);
concs=[13.42;16.77;20.97;26.21;32.76;40.96;51.20;64.00;80.00];
csvwrite('concs.csv',concs);
csvwrite('diff.csv',c_dif);
csvwrite('mask.csv',mask_iso);
cd(here);

%% Isovue Avg
clear; close all;
data_folder = 'D:\Users\Joey D\Desktop\Joey\GitHub\Machine-Learning-4-MRI\Imaging_of_pH\2-Isovue\Data\';
expected_cest_range = [-5,10];
[avg_centered,avg_xpred,avg_pH,avg_conc,~,avg_rsq,avg_signal,avg_pars,~,avg_ppm]=fit_CEST_all_pH_iso(data_folder,expected_cest_range);
here=pwd;
avg_pars = avg_pars(linspace(1,72,72)~=40,:);
avg_pH = avg_pH(linspace(1,72,72)~=40);
avg_conc = avg_conc(linspace(1,72,72)~=40);
avg_ratio_log = log10(avg_pars(:,2)./avg_pars(:,4));
cd('D:\Users\Joey D\Desktop\Joey\GitHub\Machine-Learning-4-MRI\Imaging_of_pH\3-Oxilan\Code\ML_python\Isovue');
csvwrite('avg_pars.csv',avg_pars);
csvwrite('avg_pH.csv',avg_pH);
csvwrite('avg_conc.csv',avg_conc);
csvwrite('avg_ratio_log.csv',avg_ratio_log);
cd(here);


%% SA Voxel
clear;
load('all_voxel_centered_spectra_no_T2_SA.mat');
load('mask_SA.mat');
here=pwd;
cd('D:\Users\Joey D\Desktop\Joey\GitHub\Machine-Learning-4-MRI\Imaging_of_pH\3-Oxilan\Code\ML_python\SA');
csvwrite('fitted_cest.csv',c_centered);
csvwrite('unfitted_cest.csv',c_uncentered);
csvwrite('ppm.csv',c_ppm);
csvwrite('pars.csv',c_pars);
csvwrite('rsq.csv',c_rsq);
csvwrite('pH.csv',c_pH);
csvwrite('pHs.csv',pHs);
csvwrite('conc.csv',c_conc);
concs=[13.42;16.77;20.97;26.21;32.76;40.96;51.20;64.00;80.00];
csvwrite('concs.csv',concs);
csvwrite('diff.csv',c_dif);
csvwrite('mask.csv',mask_SA);
cd(here);