function [] = overlayparametricmap( AnatomicalImage,MAP,mask,range )
%UNTITLED Summary of this function goes here
%   [] = overlayparametricmap( AnatomicalImage,MAP,mask )
% AnatomicalImage = grey anatonical reference
% MAP= parametric map
% mask = binary mask
 
 
baseImage=AnatomicalImage;
parameterROIImage=MAP;
 
baseImage = baseImage/(max(baseImage(:))); % normalize base (anatomical) image
rgbSlice  = baseImage(:,:,[1 1 1]);        % converting to RGB (ignore colormaps)
figure();
imshow(rgbSlice);                      % superimpose anatomical image
hold on;
h = imagesc(parameterROIImage, range);             % show parametric image
colormap('jet');                           % apply colormap
set(h, 'AlphaData', mask);      % make pixels outside the ROI transparent
colorbar; 
 
end