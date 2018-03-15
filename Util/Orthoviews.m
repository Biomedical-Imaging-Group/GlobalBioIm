function Orthoviews(im,pos,tt,newfig)
%% function Orthoviews(im,pos,tt,newfig)
%
% Display a orthoviews of a 3D image with gray colormap, removing the grid 
% and keeping the proportions of the image when displaying. 
% - im is the image to display
% - pos is a vector containing the position to extract planes in (i,j,k)
%   coordinates (i.e. (row,colunm,frame))
% - tt is for a global tilte
% - newfig is a boolean true to create a new figure (default true)
%
% -- Example
% load mri.mat; D=squeeze(D);
% Orthoviews(D);
% Orthoviews(D,[30,70,10]);
%
% Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch

sz=size(im);im=double(im)/double(max(im(:)));
margin=10;
if nargin <2 || isempty(pos), pos=floor(sz/2); end;
if nargin <3, tt=' '; end;
if nargin <4, newfig=1; end;
if newfig, figure; end
assert(length(sz)==3,'Wrong length of pos vector');
assert(all(pos>0) && all(pos<=sz),'Given position out of image range');
imm=ones_(sz(1)+sz(3)+margin,sz(2)+sz(3)+margin)*0.93;
imm(1:sz(1),1:sz(2))=im(:,:,pos(3));
imm(sz(1)+margin+1:end,1:sz(2))=transpose(squeeze(im(pos(1),:,:)));
imm(1:sz(1),sz(2)+margin+1:end)=squeeze(im(:,pos(2),:));
imagesc(imm); axis image; axis off;colormap gray;caxis([0,1])
line([pos(2) pos(2)],[1 sz(1)],'Color','b','LineWidth',2);
line([1 sz(2)],[pos(1) pos(1)],'Color','b','LineWidth',2);
line([sz(2)+margin+pos(3) sz(2)+margin+pos(3)],[1 sz(1)],'Color','b','LineWidth',2);
line([sz(2)+margin+1 sz(2)+margin+sz(3)+1],[pos(1) pos(1)],'Color','b','LineWidth',2);
line([1 sz(2)],[sz(1)+margin+pos(3) sz(1)+margin+pos(3)],'Color','b','LineWidth',2);
line([pos(2) pos(2)],[sz(1)+margin+1 sz(1)+margin+sz(3)],'Color','b','LineWidth',2);
text(-5,pos(1),['I=',num2str(pos(1))],'FontSize',14,'HorizontalAlignment','right')
text(pos(2),-7,['J=',num2str(pos(2))],'FontSize',14,'HorizontalAlignment','center')
text(sz(2)+margin+pos(3),-7,['K=',num2str(pos(3))],'FontSize',14,'HorizontalAlignment','center')
text(-5,sz(1)+margin+pos(3),['K=',num2str(pos(3))],'FontSize',14,'HorizontalAlignment','right')

ax=axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
set(get(ax,'Title'),'Visible','on','FontSize',14)
title(tt);

end
