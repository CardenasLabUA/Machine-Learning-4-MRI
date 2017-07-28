clear; clc; close all;
chosen_pH=7.02;
cest_dir='C:\Users\Luis\Desktop\Oxilan\Images\pH7.02\CEST';
t2_dir='C:\Users\Luis\Desktop\Oxilan\Images\pH7.02\T2';
cest_img = reshape(Bruker_reader(cest_dir),128,128,105);
anatomical=cest_img(:,:,1);
t2_img = reshape(Bruker_reader(t2_dir),128,128,50);
xpred=linspace(-8,15,101);
expected_cest_range=[-8,15];

mask=t2_img(:,:,20)./max(max(t2_img(:,:,20)))>=.19;

% %%
% mask=zeros(128,128);
% imshow(anatomical,[]);
% for j=1:9
%     s_roi=roipoly;
%     mask=mask+s_roi;
% end
% %%

load('all_voxel_centered_spectra_with_T2.mat');

%% PLS Regression

% Components=2:71;
% rmse_pH=nan(size(Components));
% rmse_conc=rmse_pH;
% beta_pH_var=zeros(size(c_centered,2)+1,length(Components));
% beta_conc_var=zeros(size(c_centered,2)+1,length(Components));

Components=5;

% for j=1:length(Components)
% [~,~,~,~,beta_pH,~] = plsregress(c_centered,c_pH,Components(j),'cv',size(c_centered,1));
% beta_pH_var(:,j)=beta_pH;
% [~,~,~,~,beta_conc,~] = plsregress(c_centered,c_conc(c_rsq>0.98,:),Components(j),'cv',size(c_centered,1));
% beta_conc_var(:,j)=beta_conc;
% end

[~,~,~,~,beta_pH,~] = plsregress(c_centered,c_pH,Components,'CV',size(c_centered,1));
% yfit = [ones(size(c_centered,1),1) c_centered]*beta_pH;
[~,~,~,~,beta_conc,~] = plsregress(c_centered,c_conc,Components,'CV',size(c_centered,1));


indices=find(reshape(mask,[],1));
cest_img=reshape(cest_img,[],105);
rsq_map=nan(128,128);
% pH_map=nan(128,128,length(Components));
% conc_map=pH_map;
pH_map=rsq_map;
conc_map=pH_map;

for q=1:length(indices)
    s=cest_img(indices(q),:);
    s=s(5:end);
    [~,i]=max(1-s);
    ppm=xpred-xpred(i);
    first=find(round(ppm)==expected_cest_range(1),1);
    if ~isempty(first)
        s=s./s(first);
        control=first;
    else
        s=s./s(1);
        control=1;
    end
    
    %% Center data via Lorentzian fitting
    % Prescribe number of pools for lorentzian fit
    method.Npools=3;

    % Prescribe initial guesses and lower and upper bounds for amplitudes,
    % widths, and offsets of lorentzian peaks from the 3 pools
    method.x0=[0.4, 1, 1.5;... Water amplitude
    0, 0.3, 1;... Oxilan pool 1 amplitude
    0, 0.2, 1;... Oxilan hydroxyls amplitude
    0, 4, 10;... Water width
    0, 4, 15;... Oxilan pool 1 width
    0, 2, 15;... Oxilan hydroxyls width
    -0.5, 0, 0.5;... Water offset
    3.5, 4.1, 5.5;... Oxilan pool 1 offset
    0.5, 2.5, 3.5];... Oxilan hydroxyls offset
    
    method.range=[1,length(s)];
    fit=cf_Lorentzian(s,ppm,method,s(control));
    fit.pars(7:9)=fit.pars(7:9)-fit.pars(7);
    centered=lorentzian(fit.pars,xpred)';
    rsq_map(indices(q))=fit.rsq;
    
    %% Apply PLS model
%     for k=1:length(Components)
%     pH_map_dummy=pH_map(:,:,k);
%     pH_map_dummy(indices(q))=[1 centered]*beta_pH_var(:,k);
%     pH_map(:,:,k)=pH_map_dummy;
%     conc_map_dummy=conc_map(:,:,k);
%     conc_map_dummy(indices(q))=[1 centered]*beta_conc_var(:,k);
%     conc_map(:,:,k)=conc_map_dummy;
%     end
    pH_map(indices(q))=[1 centered]*beta_pH;
    conc_map(indices(q))=[1 centered]*beta_conc;
end

%% Plot
figure(); imagesc(pH_map,[0 10]); colormap('jet'); colorbar;
figure(); imagesc(conc_map,[0 170]); colormap('jet'); colorbar;
figure(); imagesc(rsq_map,[0 1]); colormap('jet'); colorbar;

%% Determine concentration of each voxel
scan_h_mask=zeros(128,128);
scan_v_mask=scan_h_mask;
load('conc.mat');
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

%% Find RMSE
mask=logical(mask);
rmse_pH=(sum((chosen_pH-pH_map(mask)).^2)./ length(pH_map(mask)))^(1/2) ./ mean(pH_map(mask))*100;
rmse_conc=(sum((conc_mask(mask)-conc_map(mask)).^2)./ length(conc_map(mask)))^(1/2) ./ mean(conc_map(mask))*100;

% for p=1:length(Components)
%     pH_map_dummy=pH_map(:,:,p);
%     rmse_pH(p)=(sum((7.02-pH_map_dummy(mask)).^2)./ length(pH_map_dummy(mask)))^(1/2);
%     conc_map_dummy=conc_map(:,:,p);
%     rmse_conc(p)=(sum((conc_mask(mask)-conc_map_dummy(mask)).^2)./ length(conc_map_dummy(mask)))^(1/2);
% end

%%
% figure();
% plot(Components,rmse_pH); ylim([0,1]);
% figure();
% plot(Components,rmse_conc); ylim([0,50]);

%% Parametric overlay
overlayparametricmap(anatomical,pH_map,mask,[5,9]); title(['pH = ',num2str(chosen_pH),'  RMSE = ',num2str(rmse_pH)]); set(gcf,'Position',[1 41 1920 963]);
c_bar=get(gcf,'Children');
c_bar=c_bar(1);
set(c_bar,'FontSize',24);
%print(['By_voxel_pH_map_',num2str(chosen_pH)], '-dtiff', '-r300');
overlayparametricmap(anatomical,conc_map,mask,[15,170]); title(['pH = ',num2str(chosen_pH),'  RMSE = ',num2str(rmse_conc)]); set(gcf,'Position',[1 41 1920 963]);
c_bar=get(gcf,'Children');
c_bar=c_bar(1);
set(c_bar,'FontSize',36);
%print(['By_voxel_conc_map_',num2str(chosen_pH)], '-dtiff', '-r300');