function x = steepest(A, b, x, iMax ,W)
% solves Ax = b with the steepest descent method 
%
% from https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
% appendix B1


HtH = false;
if (nargin>4) && (~isempty(W))
    HtH = true;
end
if HtH
    r = b - A.HtWH(x,W);
else
    r = b - A*x;
end

delta = r(:)'*r(:);
delta0 = delta;
for i = 1:iMax
	if delta <= eps*delta0
		break;
	end
	if HtH
		q = A.HtWH(r,W);
	else
		q = A*r;
	end
	alpha = delta/(r(:)'*q(:));
	x = x + alpha * r;
	r = r - alpha * q;
	delta = r(:)'*r(:);
	
end