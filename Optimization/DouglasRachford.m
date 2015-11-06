function [lkl,RMS,x] = DouglasRachford(Prox1, Prox2,L,y0,maxiter,gamma, lambda,Ref)
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
lkl = zeros([maxiter,1]);
RMS = zeros([maxiter,1]);
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
        x = L.Adjoint( Prox1.Apply(Ly, gamma));
       else
        x = y + 1./nu.* L.Adjoint( Prox1.Apply(Ly, nu.*gamma) - Ly);
       end
    else
        x = Prox1.Apply(y, gamma);
    end
    y = y + lambda .* ( Prox2.Apply(2.*x- y,gamma) - x);
    tmp =  Prox2.Apply(x,gamma) - x;
    lkl(n) = sum(abs(x(:) - x_prev(:)).^2 +abs(tmp(:)).^2);
    RMS(n) = sqrt(sum(abs(Ref(:) - x(:)).^2));
end
end