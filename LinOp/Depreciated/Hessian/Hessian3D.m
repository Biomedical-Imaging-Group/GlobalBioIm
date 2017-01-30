function Hf=Hessian3D(F,bc)
%F denotes one component of the 3D vector field. This should be of size MxNxK. The 
%three indices indicate the spatial coordinates. 

[r,c,s]=size(F);

Hf=zeros(r,c,s,6); %Hessian of the vector field component. For every spatial 
%coordinate we store the Hessian which is a 6x1 matrix.
%Since Hf is symmetric we can just save the upper triangular
%part of every vector field component. 
%The convention is that the last dimension holds the elements of the 
%Hessian in a column-wise fashion, i.e., 
%Hf(m,n,l,:)=[d^2F/dxx;d^2F/dxy;d^2F/dxz;d^2F/dyy;d^2F/dyz;d^2F/dzz];


%d^2F_k/dxx
Hf(:,:,:,1)=(F-2*shift(F,[-1,0,0],bc)+shift(F,[-2,0,0],bc));
%d^2F_k/dxy
Hf(:,:,:,2)=(F-shift(F,[0,-1,0],bc)-shift(F,[-1,0,0],bc)+shift(F,[-1,-1,0],bc));
%d^2F_k/dxz
Hf(:,:,:,3)=(F-shift(F,[0,0,-1],bc)-shift(F,[-1,0,0],bc)+shift(F,[-1,0,-1],bc));
%d^2F_k/dyy
Hf(:,:,:,4)=(F-2*shift(F,[0,-1,0],bc)+shift(F,[0,-2,0],bc));
%d^2F_k/dyz
Hf(:,:,:,5)=(F-shift(F,[0,0,-1],bc)-shift(F,[0,-1,0],bc)+shift(F,[0,-1,-1],bc));
%d^2F_k/dzz
Hf(:,:,:,6)=(F-2*shift(F,[0,0,-1],bc)+shift(F,[0,0,-2],bc));





