
classdef Identity <  LinOp
    %% Identity : identity operator
    %  Matlab Linear Operator Library 
    % Build the identity operator such that x = Id.apply(x)
    %
    % Example
    % Id = Identity()
    %
    % Please refer to the LinOp superclass for documentation
    % See also LinOp
     
%     Copyright (C) 2015 F. Soulez  ferreol.soulez@epfl.ch
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

    methods
        function this = Identity()
            this.name='Identity';
            this.iscomplex=false;
        end
        function y = Apply(this,x)
            this.sizeout=size(x);
            this.sizein=size(x);
            y =x;
        end
        function y = Adjoint(this,x)
            this.sizeout=size(x);
            this.sizein=size(x);
            y =x;
        end
        function y = Gram(~,x)
            y =x;
        end
        function y = Inverse(~,x)
            y =x;
        end
        function y = AdjointInverse(~,x)
            y =x;
        end
    end
end

