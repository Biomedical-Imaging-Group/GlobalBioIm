function imdisp(im,titl,newfig)
%% function imdisp(im,titl,newfig)
%
% Display a 2D image with gray colormap, removing the grid 
% and keeping the proportions of the image when displaying
% - im is the image to display
% - titl is the title of the figure
% - newfig is a boolean true to create a new figure (default true)
%
% -- Example
% im=imread('cameraman.tif');
% imdisp(im,'Cameraman',1);
%
% Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch

if nargin ==2, newfig=1; end;
if newfig, figure; end
imagesc(im); colormap gray; axis off; axis image;
title(titl);set(gca,'FontSize',14);
end
