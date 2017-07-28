uiopen();

DCE = squeeze(DCE(:,:,1));
imshow(DCE,[]);


print('anatomical', '-dtiff', '-r300');
close all;

tumor = mask.TOI;
cbar = false(192,192);

for j = 1:length(TempRes)
    par = imgaussfilt(LRRMmaps(j).RKtrans.*tumor);
    h=imshow(par,[0,10]); colormap(jet(256));
    set(h,'AlphaData',tumor);
    
%     print(['LRRM_RKtrans_',num2str(TempRes(j))], '-dtiff', '-r300');
%     close all;
    set(h,'AlphaData',cbar);
    colorbar;
    print('RKtrans_colorbar','-dtiff','-r300');
    close all;
    
    par = imgaussfilt(LRRMmaps(j).kepTOI.*tumor);
    h=imshow(par,[0,2]); colormap(jet(256));
    set(h,'AlphaData',tumor);
    
%     print(['LRRM_kepTOI_',num2str(TempRes(j))], '-dtiff', '-r300');
%     close all;
    set(h,'AlphaData',cbar);
    colorbar;
    print('kepTOI_colorbar','-dtiff','-r300');
    close all;
    
%     par = imgaussfilt(NRRMmaps(j).RKtrans.*tumor);
%     h=imshow(par,[0,10]); colormap(jet(256));
%     set(h,'AlphaData',tumor);
%     
%     print(['NRRM_RKtrans_',num2str(TempRes(j))], '-dtiff', '-r300');
%     close all;
%     
%     par = imgaussfilt(NRRMmaps(j).kepTOI.*tumor);
%     h=imshow(par,[0,2]); colormap(jet(256));
%     set(h,'AlphaData',tumor);
%     
%     print(['NRRM_kepTOI_',num2str(TempRes(j))], '-dtiff', '-r300');
%     close all;
    
    par = imgaussfilt(LTMmaps(j).Ktrans.*tumor,0.2);
    h=imshow(par,[0,7]); colormap(jet(256));
    set(h,'AlphaData',tumor);
    
%     print(['LTM_Ktrans_',num2str(TempRes(j))], '-dtiff', '-r300');
%     close all;

    set(h,'AlphaData',cbar);
    colorbar;
    print('Ktrans_colorbar','-dtiff','-r300');
    close all;
    
%     par = imgaussfilt(LTMmaps(j).kepTOI.*tumor,0.2);
%     h=imshow(par,[0,2]); colormap(jet(256));
%     set(h,'AlphaData',tumor);
%     
%     print(['LTM_kepTOI_',num2str(TempRes(j))], '-dtiff', '-r300');
%     close all;
%     
%     par = imgaussfilt(NLTMmaps(j).Ktrans.*tumor,0.2);
%     h=imshow(par,[0,7]); colormap(jet(256));
%     set(h,'AlphaData',tumor);
%     
%     print(['NLTM_Ktrans_',num2str(TempRes(j))], '-dtiff', '-r300');
%     close all;
%     
%     par = imgaussfilt(NLTMmaps(j).kepTOI.*tumor,0.2);
%     h=imshow(par,[0,2]); colormap(jet(256));
%     set(h,'AlphaData',tumor);
%     
%     print(['NLTM_kepTOI_',num2str(TempRes(j))], '-dtiff', '-r300');
%     close all;
    
end