
classdef LinOpIdentity <  LinOp
    %% LinOpIdentity : identity operator
    %  Matlab Linear Operator Library
    % Build the identity operator such that x = Id.apply(x)
    %
    % Example
    % Id = LinOpIdentity()
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
        function this = LinOpIdentity(sz)
            this.name='LinOp Identity';
            this.iscomplex=false;
            this.isinvertible=true;
            this.issquare = true;
            this.norm=1;
            if nargin>0
                this.sizeout=sz;
                this.sizein=sz;
            end
        end
        function y = apply(this,x)
            this.sizeout=size(x);
            this.sizein=size(x);
            y =x;
        end
        function y = adjoint(this,x)
            this.sizeout=size(x);
            this.sizein=size(x);
            y =x;
        end
        function y = HtH(~,x)
            y =x;
        end
        function y = inverse(~,x)
            y =x;
        end
        function y = adjointInverse(~,x)
            y =x;
        end
    end
end

