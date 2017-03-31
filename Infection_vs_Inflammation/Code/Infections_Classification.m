clear all;
clc;
Data_Structure=importdata('CEST_infections.csv');
%%
X=Data_Structure.data;
X_infection= X(X(:,end) == 0,  :);
X_ster_inf= X(X(:,end) == 2,  :);
Xfinal=[X_infection(:,1:end-1);X_ster_inf(:,1:end-1)];
%%
Table_Infec_vs_Inf=array2table(Xfinal);
Table_Infec_vs_Inf.Tissue=categorical([X_infection(:,end);X_ster_inf(:,end)]);

%%
[coeff,score,latent,tsquared,explained,mu] = pca(Xfinal,'NumComponents',3);
%%
G=Table_Infec_vs_Inf.Tissue;
subplot(1,2,1); plot(cumsum(explained),'o-');
subplot(1,2,2); gscatter(score(:,1),score(:,2),G,[],[],50)
%%
[coeff,score,latent,tsquared,explained,mu] = pca(X(:,1:end-1),'NumComponents',3);
figure(2); scatter3(score(:,1),score(:,2),score(:,3),100,X(:,end),'filled'); colorbar;
colormap(jet(4));
%%
G=Table_Infec_vs_Inf.Tissue;
[coeff,score,latent,tsquared,explained,mu] = pca(Xfinal,'NumComponents',3);
figure(); scatter3(score(:,1),score(:,2),score(:,3),100,G,'filled'); colorbar;
colormap(jet(2));
%%
gscatter(ratings(:,1),ratings(:,2),group,'br','xo')
gscatter