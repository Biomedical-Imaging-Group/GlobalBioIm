function F=AdjHessian3D(Hf,bc)
%Hf denotes the Hessian of a component of the 3D vector field. This should
%be of size MxNxKx6.
%The first three indices indicate the spatial coordinates. The last dimension
%holds the Hessian components for a specific coordinate.
%The convention is that the last dimension holds the elements of the
%Hessian in a column-wise fashion, i.e.,
%Hf(m,n,l,:)=[d^2F/dxx;d^2F/dxy;d^2F/dxz;d^2F/dyy;d^2F/dyz;d^2F/dzz];

[r,c,s,v]=size(Hf);

if ~isequal(v,6)
  error('AdjHessian2D: The first input argument is not a valid Hessian of a 3D vector field component');
end

%F=zeros(r,c,s); %A 3D vector field component.

Pxx=Hf(:,:,:,1);
Pxx=(Pxx-2*shiftAdj(Pxx,[-1,0,0],bc)+shiftAdj(Pxx,[-2,0,0],bc));
Pxy=Hf(:,:,:,2);
Pxy=(Pxy-shiftAdj(Pxy,[0,-1,0],bc)-shiftAdj(Pxy,[-1,0,0],bc)+shiftAdj(Pxy,[-1,-1,0],bc));
Pxz=Hf(:,:,:,3);
Pxz=(Pxz-shiftAdj(Pxz,[0,0,-1],bc)-shiftAdj(Pxz,[-1,0,0],bc)+shiftAdj(Pxz,[-1,0,-1],bc));
Pyy=Hf(:,:,:,4);
Pyy=(Pyy-2*shiftAdj(Pyy,[0,-1,0],bc)+shiftAdj(Pyy,[0,-2,0],bc));
Pyz=Hf(:,:,:,5);
Pyz=(Pyz-shiftAdj(Pyz,[0,0,-1],bc)-shiftAdj(Pyz,[0,-1,0],bc)+shiftAdj(Pyz,[0,-1,-1],bc));
Pzz=Hf(:,:,:,6);
Pzz=(Pzz-2*shiftAdj(Pzz,[0,0,-1],bc)+shiftAdj(Pzz,[0,0,-2],bc));
F=Pxx+Pyy+Pzz+2*(Pxy+Pxz+Pyz);
