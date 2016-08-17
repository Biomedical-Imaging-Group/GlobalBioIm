function x=ADMM_Restore(H, proxOps, W, y, proxes, mu ,rhos, x0,maxiter,cgmaxiter)
% function x=ADMM_Restore(H,D,B,W,y, zprox, tProx, mu ,rho1, rho2, x0,maxiter,cgmaxiter)

x = x0;

if ~iscell(proxOps)
	proxOps = {proxOps};
end
if ~iscell(proxes)
	proxes = {proxes};
end

fprintf('******************************************\n');
fprintf('**  ADMM Restore   **\n');
fprintf('******************************************\n');
fprintf('#iter  Likelihood   \t primal norm z    dual norm z \t primal norm t    dual norm t    \n')
fprintf('====================================================================\n');

A = OneToMany([{H}, proxOps], [1, rhos]);
Wy = W*y;

zs = cell(1, length(proxes)); % z t
for i = 1:length(proxes)
	zs{i} = zeros(proxOps{i}.sizeout);
end
us = zs; % u1 u2
z_prevs = zs;
resids = zs;

zus = cell(1, length(proxes)); % zu1, tu2


for k=1:maxiter
% Sub problem 1
% || Hx - y||_w^2 + rho1/2 || Dx - z + u1/rho1 ||_2^2 + rho2/2 || Bx - t + u2/rho2 ||_2^2

for i = 1:length(proxes)
	zus{i} = zs{i} - us{i}/rhos(i);
end
%b = H'* wy + rho1*D'*zu1 + rho2 *B'*tu2;
b = A.Adjoint([Wy, zus]);
x = ConjGrad(A,b,x,cgmaxiter, [{W}, num2cell(ones(1, length(proxes)))] ); 

lkl = norm( reshape(W*(H*x - y),numel(y),1)); % likelihood

% Sub problems
for i = 1:length(proxes)
	z_prevs{i} = zs{i};
	Dx = proxOps{i} * x;
	xu = Dx + us{i} / rhos(1);
	zs{i} = proxes{i}.Apply( xu, mu/rhos(i) );
	
	resids{i} = Dx - zs{i};
	us{i} = us{i} + rhos(i) * resids{i};
end


for i = 1:length(proxes)
	rnorm(i) = norm( reshape(resids{i},numel(zs{i}),1));% Primal residual norm
	snorm(i) = norm( reshape(rhos(i) * proxOps{i}'*(zs{i} - z_prevs{i}),numel(x),1));% Dual residual norm
end

fprintf('%3d \t%12.6g \t', k,lkl);
for i = 1:length(proxes)
	fprintf('%12.6g \t%12.6g \t', rnorm(i),snorm(i));
end
fprintf('\n');
end
end