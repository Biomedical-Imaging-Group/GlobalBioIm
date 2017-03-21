%% Conjugate gradient
%
% Example
% x = ConjGrad(A,b,  x0,maxiter,W)
% Estimate the solution of the eqution A*x = b with XO the initilization
% and maxiter the maximum number of iteration.
% if parameters W is given (  x = ConjGrad(A,b,  x0,maxiter, W)
% it solves  A.HtWH(x) = b

%     Copyright (C) 2015 F. Soulez ferreol.soulez@epfl.ch
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or (at
%     your option) any later version.
%
%     This program is distributed in the hope that it will be useful, but
%     WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

function x = ConjGrad(A,b,  x0,maxiter,W)
HtH = false;
if (nargin>4) && (~isempty(W))
    HtH = true;
end
x = x0;% starting point
if HtH
    r = b - A.HtWH(x,W);
else
    r = b - A*x;
end
k = 0;
while k< maxiter
    rho = dot(r(:), r(:));
    if k == 0
        tol = eps*rho;
        p = r;
    elseif rho <= tol
        break
    else
        beta = rho/rho_prec;
        p = r + beta*p;
    end
    if HtH
        q = A.HtWH(p,W);
    else
        q = A*p;
    end
    
    alpha = rho/dot(p(:), q(:));
    x = x + alpha*p;
    r = r - alpha*q;
    rho_prec = rho;
    k = k + 1;
end
