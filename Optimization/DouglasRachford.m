function [lkl,RMS,cost,bestx] = DouglasRachford(Prox1, Prox2,L,y0,maxiter,gamma, lambda,Ref1,Ref2)
% Implement the Douglas Rachford splitting algorithm that solves:
% x^+ = argmin_x( f1(x) + f2(L x))
% with Prox1 and Prox2 the proximal operator of f1 and f2 respectively
% L a linear operator  such  L.L^T = nu I
% y0 the initialization
% maxiter the number of iteration
% gamma \in [0,+\inf[
% lambda\in ]0,2[ the relaxation parmeter.
%
nu =1;
useL = 0;
bestx = y0;
Lx = L*y0;
bestrms = sqrt(sum(abs(Ref1(:) - y0(:)).^2) + sum(abs(Ref2(:) - Lx(:)).^2));

lkl = zeros([maxiter,1]);
RMS = zeros([maxiter,1]);
cost = zeros([maxiter,1]);
if isa(L,'LinOp')
    useL = 1;
    r = randn(L.sizeout);
    nu = r ./ L.HHt(r);
    assert(std(nu(:)) <1e-6, 'LLt != nu I');
    nu = mean(nu(:));
    if nu==1
        useL = 2;
    end
end
y = y0;
x= y0;
for n=1:maxiter
    x_prev = x;
    y_prev = y;
    if useL
        Ly = L*y;
       if useL==2
        x = L.Adjoint( Prox2.Apply(Ly, gamma));
       else
        x = y + 1./nu.* L.Adjoint( Prox2.Apply(Ly, nu.*gamma) - Ly);
       end
    else
        x = Prox2.Apply(y, gamma);
    end
    y = y + lambda .* ( Prox1.Apply(2.*x- y,gamma) - x);
    tmp1 =  Prox1.Apply(x,gamma) - x;
    Lx = L*x;
    tmp2 =  Prox2.Apply(Lx,gamma) - Lx;
    lkl(n) = sum(abs(tmp1(:)).^2 + abs(tmp2(:)).^2);
    RMS(n) = sqrt(sum(abs(Ref1(:) - x(:)).^2) + sum(abs(Ref2(:) - Lx(:)).^2));
    cost(n) = Prox1.FCost(x) + Prox2.FCost(Lx);
    if RMS(n) < bestrms
        bestrms =RMS(n);
        bestx = x;
    end
end
end