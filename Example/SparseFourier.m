%% Script using the inv. Pb. Library for sparse fourier measurement.
%
%

% Create a sparse image
sz = [256, 256];
lvl = 0.001;  % level of sparsity of the image
im = zeros(sz);
nz = rand(sz)< lvl;
nnz = sum(nz(:));
im(nz) = (randn([nnz 1])).^2;
figure; imagesc(im);

% Create the forward model
lvl = 0.5;  % level of sparsity of the data
sel = rand(sz)< lvl;
nel = sum(sel(:));
im(sel) = randn([nel 1]);
S = Selector(sel);
F = DFT();
F.Apply(im);
H = S*F;


% Simulate the data
data = H.Apply(im);
figure; imagesc(fftshift( abs(S.Adjoint(data))));

% Build the prox 
ProxIm = L1('NonNegativity');
ProxData = L2(data);

% Optimize...
maxiter = 500;
gamma = [1e-5 , 1];
lambda = 1;
x0 = ones(sz);
x = DouglasRachford(ProxIm, ProxData,H,x0,maxiter,gamma, lambda);