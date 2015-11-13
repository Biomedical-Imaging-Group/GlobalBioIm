function CheckLinOp(LinOp)
%% CHECKLINOP function
% Matlab Linear Operator Library
% Function checking the correctness of Adjoint and Inverse implementation
% of the linear operator LinOp
%
% Example:
% I = Identity()
% CheckLinop(I) % will check if adjoint gram and inverse gives coherent results 
%
% See also LinOp


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


if isnumeric(LinOp.sizein)
    sizein = LinOp.sizein;
else
    sizein = input('What is the input size of the LinOp?');
end
sizeout = LinOp.sizeout;


In = rand(sizein);
Out = rand(sizeout);

if LinOp.iscomplex
	disp('Complex LinOp');
	In = In + 1i * rand(sizein);
	Out = Out + 1i * rand(sizeout);	
end

HIn = LinOp.Apply(In);
HtOut = LinOp.Adjoint(Out);


tol = 1e-3*max( max(max(abs(In(:))),max(abs(HIn(:)))),max(max(abs(Out(:))),max(abs(HtOut(:))))); % Tolerance for numerical equality

diff = dot(In(:) ,HtOut(:)) - dot(HIn(:) , Out(:));
normDiff = diff(:)' * diff(:);
if normDiff < tol
    disp('Adjoint OK');
else
    error('Adjoint error: <Hx.y> ~= <x.H^*y> : diff = %d', normDiff);
end

diff = LinOp.Adjoint(HIn) - LinOp.HtH(In);
normDiff = diff(:)' * diff(:);
if  normDiff < tol
    disp('HtH matrix OK');
else
    error('Error in HtH matrix computation: diff = %d', normDiff);
end

if LinOp.issquare
	diff = LinOp.Apply(HtOut) - LinOp.HHt(HtOut);
	if  diff(:)' * diff(:) < tol
		disp('HHt matrix OK');
	else
		error('Error in HtH matrix computation');
	end
else
	disp('LinOp non square');
end

if LinOp.isinvertible
    if abs(LinOp.Inverse(HIn)-In)<tol
        disp('Inverse OK');
    else
        error('Error in inverse computation: H^-1 H ~= I');
    end
    HIIn = LinOp.AdjointInverse(In);
    HItOut = LinOp.Inverse(Out);
    if abs(dot(In(:) ,HItOut(:) ) - dot(HIIn(:) , Out(:)))<10*tol
        disp('Adjoint Inverse  OK');
    else
        error('Adjoint Inverse error: <H^-1 x.y> ~= <x.H^-* y>');
    end
else
    disp('LinOp non invertible');
end


end
