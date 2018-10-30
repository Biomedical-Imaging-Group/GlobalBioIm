%%
m = rand(15, 15) + 1i * rand(15, 15);
H = LinOpMatrix(m);

normEst = estimateNorm(H, 999, 300); % very accurate, slow

assert(abs(H.norm - normEst) < 1e-14);

%%
[im,psf,y]=GenerateData('Gaussian',20);
H=LinOpConv(fft2(psf));

normEst = estimateNorm(H, 999, 300); % very accurate, slow

assert(abs(H.norm - normEst) < 1e-14);

%%
hHat = rand(21, 21);
H=LinOpConv(hHat);

[normEst, v] = estimateNorm(H, 10000, 300, rand(H.sizein) + 1i * rand(H.sizein)); % very accurate, slow


assert(abs(H.norm - normEst) < 1e-6); % this one is hard to make accurate