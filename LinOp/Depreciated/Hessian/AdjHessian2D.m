function F=AdjHessian2D(Hf,bc)
%Hf denotes the Hessian of a vector field component. This should be of size MxNx3.
%The first two indices indicate the spatial coordinates. The last dimension
%holds the Hessian components for a specific coordinate.
%The convention is that the last dimension holds the elements of the
%Hessian in a column-wise fashion, i.e.,
%Hf(m,n,:)=[d^2F/dxx;d^2F/dxy;d^2F/dyy;d^2F/dxx;d^2F/dxy;d^2F/dyy];

[r,c,v]=size(Hf);

if ~isequal(v,3)
  error('AdjVHessian2D: The first input argument is not a valid Hessian of a 2D vector field');
end

%F=zeros(r,c); %A 2D vector field. For every spatial coordinate we store
%the two vector components of the vector field F.


Pxx=Hf(:,:,1);
Pxx=(Pxx-2*shiftAdj(Pxx,[-1,0],bc)+shiftAdj(Pxx,[-2,0],bc));
Pxy=Hf(:,:,2);
Pxy=(Pxy-shiftAdj(Pxy,[0,-1],bc)-shiftAdj(Pxy,[-1,0],bc)+...
  shiftAdj(Pxy,[-1,-1],bc));
Pyy=Hf(:,:,3);
Pyy=(Pyy-2*shiftAdj(Pyy,[0,-1],bc)+shiftAdj(Pyy,[0,-2],bc));
F=Pxx+2*Pxy+Pyy;
