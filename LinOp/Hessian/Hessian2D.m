function Hf=Hessian2D(F,bc)
%F denotes the vector field component. This should be of size MxN. The
%indices indicate the spatial coordinates.

[r,c]=size(F);


Hf=zeros(r,c,3); %Hessian of the vector field component. For every spatial
%coordinate we store the Hessian which is a 3x1 matrix.
%Since Hfk is symmetric we can just save the upper triangular
%part of every Hfk.
%The convention is that the last dimension holds the elements of the
%Hessian in a column-wise fashion, i.e.,
%Hf(m,n,:)=[d^2F/dxx;d^2F/dxy;d^2F/dyy];

%d^2F_k/dxx
Hf(:,:,1)=(F-2*shift(F,[-1,0],bc)+shift(F,[-2,0],bc));
%d^2F_k/dxy
Hf(:,:,2)=(F-shift(F,[0,-1],bc)-shift(F,[-1,0],bc)+shift(F,[-1,-1],bc));
%d^2F_k/dyy
Hf(:,:,3)=(F-2*shift(F,[0,-1],bc)+shift(F,[0,-2],bc));





