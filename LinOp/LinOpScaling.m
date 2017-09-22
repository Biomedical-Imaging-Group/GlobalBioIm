classdef LinOpScaling <  LinOp
    % LinOpScaling : LinOp scaling operator
    %
    % Scale the input by a scalar factor
    %
    % :param sz: input size of the operator
    % :param scale: scaling parameter
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Example** S= LinOpScaling(scale)
    %
    % See also :class:`LinOp`, :class:`Map`
    
    %%    Copyright (C) 2015 
    %     F. Soulez ferreol.soulez@epfl.ch
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
    
    %% Constructor
    methods
        function this = LinOpScaling(scale,sz)
            assert(issize(sz),'The input size sz should be a conformable  to a size ');   
            this.name ='LinOpScaling';
            if isscalar(scale)
                this.scale = scale;
                if scale == 0
                    this.isInvertible= false;
                end
            else
                error('LinOpScale value must be a scalar');
            end
            this.isDifferentiable=true;                   
            this.sizeout=sz;
            this.sizein=sz;          
		end
    end
	
    %% Core Methods containing implementations (Protected)
	methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y =this.scale .* x;
        end
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y =conj(this.scale) .*x;
        end
        function y = applyInverse_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if this.isInvertible
                y =(1./this.scale) .*x;
            else
                y = applyInverse_@LinOp(this,x);
            end
        end
        function y = applyAdjointInverse_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            if this.isInvertible
                y =conj(1./this.scale) .*x;
            else
                y = applyAdjointInverse_@LinOp(this,x);
            end           
        end
    end
end
