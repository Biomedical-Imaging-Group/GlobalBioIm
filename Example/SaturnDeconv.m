saturn = fitsread('saturn.fits');
psf = fitsread('saturn_psf.fits');
psf = padarray(psf,[230,230]);
psf = psf(1:end-1, 1:end-1);
psf = psf/ sum(sum(psf));
cnvl = Convolution(fftshift(psf));
D = Grad(size(saturn));
mu =1e-2;

A = cnvl'*cnvl + mu * D'*D;
b = cnvl'* saturn;
x = ConjGrad(A,b,  zeros(size(saturn)),100);

figure;imshow(saturn,[0,max(max(saturn))]);
figure;imshow(x,[0,max(max(x))])