% close all; clc;
% %mask_SA=zeros(128,128,8);
% load('all_voxel_centered_spectra_with_T2_SA.mat');
% for p=5:length(pHs)
%     cest_img=reshape(Bruker_reader(['C:\Users\Luis\Desktop\Oxilan\SA\pH',num2str(pHs(p),'%.2f'),'\CEST']),128,128,105);
%     figure();
%     imshow(cest_img(:,:,1),[]);
%     for c=1:9
%         mask_dum=roipoly;
%         mask_SA(:,:,p)=mask_SA(:,:,p)+mask_dum;
%     end
%     close;
% end
% mask_SA=logical(mask_SA);
% save('mask_SA.mat','mask_SA');

%mask_iso=zeros(128,128,8);
% load('all_voxel_centered_spectra_with_T2_isovue.mat');
% for p=8:length(pHs)
%     cest_img=reshape(Bruker_reader(['C:\Users\Luis\Desktop\Oxilan\Isovue\pH',num2str(pHs(p),'%.2f'),'\CEST']),128,128,105);
%     figure();
%     imshow(cest_img(:,:,1),[]);
%     for c=1:9
%         mask_dum=roipoly;
%         mask_iso(:,:,p)=mask_iso(:,:,p)+mask_dum;
%     end
%     close;
% end
% mask_iso=logical(mask_iso);
% save('mask_iso.mat','mask_iso');

% %%
% figure();
% for p=1:8
%     subplot(3,3,p); imshow(mask_iso(:,:,p),[]);
% end

mask_oxi=zeros(128,128,8);
load('all_voxel_centered_spectra_with_T2.mat');
load('conc.mat');
for p=1:length(pHs)
    cest_img=reshape(Bruker_reader(['D:\Users\Joey D\Desktop\Joey\GitHub\Machine-Learning-4-MRI\Imaging_of_pH\3-Oxilan\Code\ML_python\Oxilan\Images\pH',num2str(pHs(p),'%.2f'),'\1']),128,128,105);
    t2_img=reshape(Bruker_reader(['D:\Users\Joey D\Desktop\Joey\GitHub\Machine-Learning-4-MRI\Imaging_of_pH\3-Oxilan\Code\ML_python\Oxilan\Images\pH',num2str(pHs(p),'%.2f'),'\2']),128,128,50);
    mask=t2_img(:,:,20)./max(max(t2_img(:,:,20)))>=.19;
    
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
    mask_oxi(:,:,p)=mask;
end

mask_oxi=logical(mask_oxi);
save('mask_oxi.mat','mask_oxi');

%%
figure();
for p=1:8
    subplot(3,3,p); imshow(mask_oxi(:,:,p),[]);
end
