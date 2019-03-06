x.size = [125 125];
x.step = [1 1];

numThetas = 100;
thetas = linspace(0, pi, numThetas+1);
thetas(end) = []; % we don't want to have both 0 and pi

H = LinOpXRay(x, thetas);

% make some arbitrary phantom
[X, Y] = ndgrid(1:x.size(1), 1:x.size(2));
c = (X-50).^2 + (Y-50).^2 < 40^2;
c = c & (((X-40)/1.5).^2 + (Y-24).^2 > 10^2); 
c = c & (((X-40)/1.5).^2 + (Y-60).^2 > 20^2); 
c = c & ((X-60).^2 + (Y-35).^2 > 5^2); 
c = c & ((X-77).^2 + ((Y-35)./1.5).^2 > 8^2); 
c = c .* (sqrt( (X - x.size(1)).^2 + (Y- x.size(2)).^2 ) );
c = reshape( double(c), x.size);

%%

g = H*c;
HTg = H'*g;

%%
d = zeros(x.size);
d(ceil(end/2), ceil(end/2)) = 1;

HTHd = H'*H*d;


%%
checkLinOp(H)


