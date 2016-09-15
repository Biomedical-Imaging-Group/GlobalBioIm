 % x = im2double(imread('rice.png'));
% H = Hessian(size(x),'reflexive');
% Hx = H*x; 
% A = Schatten(1); 
% AHx = A.Apply(Hx,10^-5);
% figure, imagesc(AHx);

x = im2double(imread('rice.png'));
x = repmat(x, [1 1 10]); 
H = Hessian(size(x),'reflexive');
Hx = H*x; 
A = Schatten(1); 
AHx = A.Apply(Hx,10^-5);
figure, imagesc(Hx(:,:,5,3)), colorbar, colormap gray;