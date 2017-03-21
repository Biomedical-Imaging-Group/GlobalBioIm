function imdisp(im,tilt,newfig)
%% function imdisp(im,tilt)
%
% Display a 2D image with gray colormap, removing the grid 
% and keeping the proportions of the image when displaying
% newfig is a boolean true to create a new figure
% -- Example
% im=imread('cameraman.tif')
% imdisp(im,'Input image',1);
%
% Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch

if newfig, figure; end
imagesc(im); colormap gray; axis off; axis image;
title(tilt);
end
