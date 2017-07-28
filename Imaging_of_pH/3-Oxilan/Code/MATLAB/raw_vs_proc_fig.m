folder_path='F:\GitHub\Machine-Learning-4-MRI\Imaging_of_pH\3-Oxilan\Data\';
center_range=[-8,15];

here=pwd;
cd(folder_path)
pH_folders=dir('pH*');
n_pH=length(pH_folders);

% Extract CEST data
last_n=0;
for p=1:n_pH
    cd([folder_path,pH_folders(p).name,'\CEST']);
    conc_files=dir('*mM*.txt');
    n_conc=length(conc_files);

     for q=1:n_conc
        s= importdata(conc_files(q).name); 
        s=s(5:end);
        signal(q+last_n,:)=s;
     end
     last_n=last_n+n_conc;
end

cd(here);

figure();
pHs = pH(1:9:end);
q=5;
cool=reshape(signal',length(c_xpred),9,8);
subplot(2,1,1,'FontSize',18); plot(c_xpred,cool(:,:,q),'--');
title('Raw Z Spectra (Ioxilan, pH 7.02)'); xlim([-8,15]);
xlabel('Offset (ppm)');
ylabel('M_z/M_o');

cool=reshape(c_centered',length(c_xpred),9,8);
subplot(2,1,2,'FontSize',18); plot(c_xpred,1-cool(:,:,q),'-','LineWidth',1.5);
title('Processed Z Spectra (Ioxilan, pH 7.02)'); ylim([0 1]); xlim([-8,15]);
xlabel('Offset (ppm)');
ylabel('M_z/M_o');

%%
figure();
pHs = pH(1:9:end);
q=5;
cool=reshape(c_signal',length(c_xpred),9,8);
plot(cool_ppm(:,:,q),cool(:,:,q),'-','LineWidth',2);
title('Raw Z Spectra (Ioxilan, pH 7.02)'); xlim([-8,15]); ylim([0,1]);
xlabel('Offset (ppm)');
ylabel('M_z/M_o');

%%
figure();
pHs = pH(1:9:end);
q=5;
cool=reshape(signal'./10^4,length(c_xpred),9,8);
subplot(1,2,1,'FontSize',18);
plot(c_xpred,cool(:,:,q),'o','MarkerSize',4);
title('Raw CEST Spectra (pH 7.02)'); xlim([-8,15]);
xlabel('Offset (ppm)'); ylim([0,2.62]);
ylabel('M_z (x 10^4)');

pHs = pH(1:9:end);
q=5;
cool=reshape(c_signal',length(c_xpred),9,8);
subplot(1,2,2,'FontSize',18);
plot(cool_ppm(:,:,q),cool(:,:,q),'o','MarkerSize',4);
title('Normalized/Corrected CEST Spectra (pH 7.02)'); xlim([-8,15]);
xlabel('Offset (ppm)'); ylim([0,1.05]);
ylabel('M_z/M_o');