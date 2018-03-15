function im=StarLikeSample(dim,N,w,s,pow,x0)
%------------------------------------------------
% function im=StarLikeSample(dim,N,w,s,pow,x0)
% 
% Generates a Star-Like object using
%
% Inputs : dim -> dimension of the ouput (2 or 3)
%          N   -> size of the square (or cube), must be an even number
%          w   -> number of branches of the object will be 4*w
%          s   -> slope of the sigmoid function attenuating the boundaries
%          pow -> power to the final image (to have smoother edges)
%          x0  -> shift of the sigmoid function 1/(1+exp(w*(x-x0)))
%
% Output : im -> generated image 
%
% Example : 
%     im=StarLikeSample(2,256,8,20,3,0.7);
%
% Copyright (C) 2018, E. Soubies emmanuel.soubies@epfl.ch
%------------------------------------------------

% -- Test inputs
assert(mod(N,2)==0,'In generateStarLikeSample : dimension N must be an even number'); 
assert(mod(w,2)==0,'In generateStarLikeSample : parameter w must be an even'); 

ll=linspace(-1,1,N);
if dim==2
	% Grid generation
	[X,Y]=meshgrid(ll,ll);
    [th,~] = cart2pol(X,Y);
    im=1+cos(4*w*th);
	im=im./(1+exp(s*(sqrt(X.^2+Y.^2)-x0)));
elseif dim==3
    % Grid generation
    [X,Y,Z]=meshgrid(ll,ll,ll);
    [~,ph1,~]=cart2sph(Y,Z,X);
    [~,ph2,~]=cart2sph(X,Y,Z);
    [~,ph3,~]=cart2sph(X,Z,Y);
    im=(1+cos(4*w*ph1)) + (1+cos(4*w*ph2)) + (1+cos(4*w*ph3));
    im=im./(1+exp(s*(sqrt(X.^2+Y.^2+Z.^2)-x0)));
end
im=max(im,0);
im=(im/2).^pow;
im=im/max(im(:));
end
