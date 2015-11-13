% example of the X-ray reconstruction using invpblib
clear all
restoredefaultpath
addpath(genpath('.'));

%% parameters

% space domain
D = 2;
N = 20;
xStart = -ones(1,D) * N;
xStep = ones(1,D);
xSize = 2*N*ones(1,D) + 1;

% projection domain
M = 40;
yStart = -ones(1,D-1) * M;
yStep = ones(1,D-1);
ySize =  2*M*ones(1,D-1) + 1;

% projections
numThetas = 50;
t = linspace(0, pi, numThetas);
P = zeros(D-1, D, numThetas);
[r1, r2] = Euler3D(t-pi/2, -pi/2 * ones(size(t)), zeros(size(t)));
R = cat(3, r2, r1); % swap is intentional
P = reshape(R(1:D, :, 1:D-1), size(P));

% discretization
phi = KaiserBesselWindow(2, 10.6, 2, D);

%% create LinOp
H = XRay(xStart, xStep, xSize, yStart, yStep, ySize, P, phi);

%% create a phantom, take measurements
ells = [1, -10, 2, 11, 5, -pi/4
.4, 5, 5, 0, 0, 0
2, 10, 7, -1, -10, 2*pi/3
2 .5 .5 -5 7 0];
phantom = EllipsoidPhantom(ells);

fTruth = phantom.eval(H.xStart, H.xStep/4, H.xSize*4 );

g = phantom.xray(H.yStart, H.yStep, H.ySize, H.Ps);

sigma = 1e0;
gNoisy = g + sigma*randn(size(g));

figure
prettyPlot(H.xStart, H.xStep, H.xSize, fTruth);
title('truth');

figure
prettyPlot([H.yStart 1], [H.yStep 1], [H.ySize H.numProj], gNoisy, {'y_0', 'y_1'});
title('sinogram')

%% CG optimization
% DECONVOLUTION
% x =argmin_x 1/2|| H.x - y||_2^2 + mu  || D.x||_2^2  
% normal equation  (H'H + mu D'D) x = H'.y
%                               A x = b
mu = 10; % hyperparameter

% Finite difference operator
D = Grad( H.xSize );

b = H' * gNoisy;
A = H'*H + mu * (D'*D);
% Solving Ax = b using conjugate gradient 
maxiter = 10;
x = ConjGrad(A, b, zeros(H.xSize), maxiter);

f = phi.reconstruct(H.xStep/4, upsample(upsample(x,4)',4)' );
figure
prettyPlot(H.xStart, H.xStep, H.xSize, f);
title('CG')

%% ADMM optimization 
maxiter = 50;
cgmaxiter = 5;

W = Identity(H.xSize); % data fit weight matrix

D = Grad(H.xSize); % LinOp for first prox
zProxes = {L2()
	L1()
	JointL1(3)};

B = Identity(H.xSize); % LinOp for second prox
tProx = NonNegativity();


rho2 = 1e4;
for mu = [1e6 1e1 1e3]
	
	rho1 = mu;
	
	for proxInd = 1:length(zProxes)
		zProx = zProxes{proxInd};
		
		
		x0 = zeros(H.xSize);
		
		
		x = ADMM_Restore(H' * H, D, B, W, H'*gNoisy , zProx, tProx, mu, rho1, rho2, x0, maxiter, cgmaxiter);
		f = phi.reconstruct(H.xStep/4, upsample(upsample(x,4)',4)' );
		SNR = snr(f, f-fTruth);
		
		
		figure
		prettyPlot(H.xStart, H.xStep, H.xSize, f);
		title(sprintf('%s of %s and %s, mu=%g. SNR = %g', zProx.name, D.name, tProx.name, rho1, SNR))
		
	end
end
