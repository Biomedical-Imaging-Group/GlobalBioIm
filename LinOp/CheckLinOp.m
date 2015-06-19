function CheckLinOp(LinOp)
%% CHECKLINOP function
% Matlab Linear Operator Library
% Function checking the correctness of Adjoint and Inverse implementation
% of the linear operator LinOp
%
% see also LinOp


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

if LinOp.iscomplex
    disp('Complex LinOp');
    In = complex(rand(sizein),rand(sizein));
    HIn = LinOp.Apply(In);
    sizeout = LinOp.sizeout;
    Out = complex(rand(sizeout),rand(sizeout));
    HtOut = LinOp.Adjoint(Out);
else
    In = rand(sizein);
    HIn = LinOp.Apply(In);
    sizeout = LinOp.sizeout;
    Out = rand(sizeout);
    HtOut = LinOp.Adjoint(Out);
end
tol = 1e-8*max( max(max(In(:)),max(HIn(:))),max(max(Out(:)),max(HtOut(:)))); % Tolerance for numerical equality

if abs(dot(In(:) ,HtOut(:)) - dot(HIn(:) , Out(:)))<tol
    disp('Adjoint OK');
else
    error('Adjoint error: <Hx.y> ~= <x.H^*y>');
end

if abs( (LinOp.Adjoint(HIn) - LinOp.Gram(In)))<tol
    disp('Gram matrix OK');
else
    error('Error in Gram matrix computation');
end

if LinOp.isinvertible
    if abs(LinOp.Inverse(HIn)-In)<tol
        disp('Inverse OK');
    else
        error('Error in inverse computation: H^-1 H ~= I');
    end
    HItOut = LinOp.AdjointInverse(Out);
    HIIn = LinOp.Inverse(In);
    if abs(dot(In(:) ,HItOut(:) ) - dot(HIIn(:) , Out(:)))<tol
        disp('Adjoint Inverse  OK');
    else
        error('Adjoint Inverse error: <H^-1 x.y> ~= <x.H^-* y>');
    end
else
    disp('LinOp non invertible');
end


end
