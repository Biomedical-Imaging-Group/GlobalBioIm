 x = im2double(imread('rice.png'));
H = Hessian(size(x),'reflexive');
Hx = H*x; 
A = Schatten(1); 
AHx = A.Apply(Hx,10^-5);
figure, imagesc(AHx);

%%
x = im2double(imread('rice.png'));
x = repmat(x, [1 1 10]); 
H = Hessian(size(x),'reflexive');
Hx = H*x; 
A = Schatten(1); 
AHx = A.Apply(Hx,10^-5);
figure, imagesc(Hx(:,:,5,3)), colorbar, colormap gray;

%%
x = im2double(imread('rice.png'));
%[dx, dy] = gradient(x);
%x = abs(dx+dy);

filter = ifftshift(fspecial('gaussian', size(x), 3));

%A = Convolution(filter);
A = Identity(size(x));

g = A*x + randn(size(x))*.1;


% Finite difference operator
%D = Grad(size(g));
%zProx = JointL1(3); % JointL1( D*x)  == Total variation 

D = Hessian(size(g), 'reflexive'); %
zProx = Schatten(2); % HS


mu = 1e-4; % hyperparameter
x0 = zeros(size(g));
cgmaxiter = 5;
maxiter =100;

rho =  mu / .01;


xHat=ADMM_Restore(A, D ,Identity(),g, zProx,mu, rho,x0,maxiter,cgmaxiter);
imagesc([g xHat x])