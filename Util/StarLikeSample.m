function im=StarLikeSample(dim,N,w,s,pow,x0)
%------------------------------------------------
% function im=StarLikeSample(dim,N,w,s,pow,x0)
% 
% Generates a Star-Like object using
%
% Inputs : dim -> dimension of the ouput (2 or 3)
%          N   -> size of the square (or cube), must be an even number
%          w   -> number of branches of the object
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

% -- Test if N is even
if mod(N,2)
	error('In generateStarLikeSample : dimension N must be an even number'); 
end

ll=linspace(-1,1,N);
if dim==2
	% Grid generation
	[X,Y]=meshgrid(ll,ll);
    [th,r] = cart2pol(X,Y);
    im=1+cos(w*th);
	im=im./(1+exp(s*(sqrt(X.^2+Y.^2)-x0)));
elseif dim==3
	% Grid generation
	[X,Y,Z]=meshgrid(ll,ll,ll);
    [th,ph,r]=cart2sph(X,Y,Z);
    im=(1+cos(w*th))+(1+cos(w*ph));
    im=im./(1+exp(s*(sqrt(X.^2+Y.^2+Z.^2)-x0)));
end
im=max(im,0);
im=(im/2).^pow;
end
