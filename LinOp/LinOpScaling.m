classdef LinOpScaling <  LinOp
    %% LinOpScaling : LinOp scaling operator
    %  Matlab Linear Operator Library
    %
    % Example
    % Obj = LinOpScaling(scale)
    %
    % Build the scaling operator that scale the input by the scalar factor
    % scale
    %
    % Please refer to the LinOp superclass for documentation
    % See also LinOp
    
    %     Copyright (C) 2015 F. Soulez ferreol.soulez@epfl.ch
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
    
    properties (SetAccess = protected,GetAccess = public)
        scale % scale factor
    end
    methods
        function this = LinOpScaling(scale,sz)
            this.name ='LinOp Scaling';
            if isscalar(scale)
                this.scale = scale;
                if isreal(scale)
                    this.isComplex=false;
                else
                    this.isComplex=true;
                end
                if scale == 0
                    this.isInvertible= false;
                end
            else
                error('LinOpScale value must be a scalar');
            end
            assert(issize(sz),'The input size sz should be a conformable  to a size ');
            
            this.sizeout=sz;
            this.sizein=sz;
            
		end
	end
	
	methods (Access = protected)
        function y = apply_(this,x)
            this.sizeout=size(x); % todo: remove this delayed size setting
            this.sizein=size(x);
            y =this.scale .* x;
        end
        function y = adjoint_(this,x)
            this.sizeout=size(x);
            this.sizein=size(x);
            if this.isComplex
                y =this.scale .*x;
            else
                y =conj(this.scale) .*x;
            end
        end
        function y = inverse_(this,x)

            y =(1./this.scale) .*x;
        end
        function y = adjointInverse_(this,x)
           
            if this.isComplex
                y = (1./this.scale) .*x;
            else
                y =conj(1./this.scale) .*x;
            end
        end
    end
end
