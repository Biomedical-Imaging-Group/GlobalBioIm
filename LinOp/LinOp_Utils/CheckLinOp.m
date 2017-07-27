function CheckLinOp(LinOp)
%% CHECKLINOP function
% Matlab Linear Operator Library
% Function checking the correctness of adjoint and Inverse implementation
% of the linear operator LinOp
%
% Example:
% I = Identity()
% CheckLinop(I) % will check if adjoint gram and Inverse gives coherent results 
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

if LinOp.isComplex
	disp('Complex LinOp');
	In = In + 1i * rand(sizein);
	Out = Out + 1i * rand(sizeout);	
end

HIn = LinOp.apply(In);
HtOut = LinOp.adjoint(Out);


tol = 1e-3*max( max(max(abs(In(:))),max(abs(HIn(:)))),max(max(abs(Out(:))),max(abs(HtOut(:))))); % Tolerance for numerical equality

diff = dot(In(:) ,HtOut(:)) - dot(HIn(:) , Out(:));
normDiff = diff(:)' * diff(:);
if normDiff < tol
    disp('adjoint OK');
else
    error('adjoint error: <Hx.y> ~= <x.H^*y> : diff = %d', normDiff);
end

diff = LinOp.adjoint(HIn) - LinOp.HtH(In);
normDiff = diff(:)' * diff(:);
if  normDiff < tol
    disp('HtH matrix OK');
else
    error('Error in HtH matrix computation: diff = %d', normDiff);
end

if  (dot((LinOp.apply(HtOut) - LinOp.HHt(Out)),conj(LinOp.apply(HtOut) - LinOp.HHt(Out)))<tol)
    disp('HHt matrix OK');
else
    error('Error in HHt matrix computation');
end


if LinOp.isInvertible
    if abs(LinOp.inverse(HIn)-In)<tol
        disp('inverse OK');
    else
        error('Error in Inverse computation: H^-1 H ~= I');
    end
    HIIn = LinOp.adjointInverse(In);
    HItOut = LinOp.inverse(Out);
    if abs(dot(In(:) ,HItOut(:) ) - dot(HIIn(:) , Out(:)))<10*tol
        disp('adjoint Inverse  OK');
    else
        error('adjoint Inverse error: <H^-1 x.y> ~= <x.H^-* y>');
    end
else
    disp('LinOp non invertible');
end


end
