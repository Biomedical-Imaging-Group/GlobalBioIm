function x=ADMM_Restore(H,D,B,W,y, zProx, tProx,mu, rho1, rho2,x0,maxiter,cgmaxiter)
x = x0;
z = zeros(D.sizeout);
t = zeros(B.sizeout);
u1 =  zeros(D.sizeout);
u2 = zeros(B.sizeout);

    fprintf('******************************************\n');
    fprintf('**  ADMM Restore   **\n');
    fprintf('******************************************\n');
    fprintf('#iter  Likelihood   \t primal norm z    dual norm z \t primal norm t    dual norm t    \n')
    fprintf('====================================================================\n');

A = OneToMany({H,D,B},[1, rho1, rho2]);
Wy = W*y;
for k=1:maxiter
% Sub problem 1
% || Hx - y||_w^2 + rho1/2 || Dx - z + u1/rho1 ||_2^2 + rho2/2 || Bx - t + u2/rho2 ||_2^2
zu1 = z - u1/rho1;
tu2 = t - u2/rho2;
%b = H'* wy + rho1*D'*zu1 + rho2 *B'*tu2;
b = A.Adjoint({Wy, zu1,tu2});
x = ConjGrad(A,b,  x,cgmaxiter,{W,1,1}); 

lkl = norm( reshape(W*(H*x - y),numel(y),1)); % likelihood

% Sub problem 2
z_prev = z;
Dx = D*x;
xu1 = Dx + u1/rho1;
z = zProx.Apply(xu1,mu/rho1);

% Sub problem 3
t_prev = t;
Bx = B*x;
xu2 = Bx + u2/rho2;
t = tProx.Apply(xu2, 1./rho2);

% Residuals of the constraints
res1 = Dx - z;
res2 = Bx - t;

% Lagrange parameters update
u1 = u1 + rho1 * res1;
u2 = u2 + rho2 * res2;

% First constraints
rnorm1 = norm( reshape(res1,numel(z),1));% Primal residual norm
snorm1 = norm( reshape(rho1 * D'*(z - z_prev),numel(x),1));% Dual residual norm 

% Second constraints
rnorm2 = norm( reshape(Bx - t,numel(t),1));
snorm2 = norm( reshape(rho2 * B'*(t - t_prev),numel(x),1));

  fprintf('%3d \t%12.6g \t%12.6g \t%12.6g \t%12.6g \t%12.6g\n',k,lkl,rnorm1,snorm1,rnorm2,snorm2);
end
end