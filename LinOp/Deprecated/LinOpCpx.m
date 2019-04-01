classdef LinOpCpx <  LinOp
    % LinOpCpx : Complex representation  operator
    %
    % Build the operator transforming a complex vector as a 2D real vector with 
    % [Real, Im] the imaginary part of the input vector
    %
    % :param sz: sizein of the operator.
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Example** C=LinOpCpx(sz)
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
    
    properties (SetAccess = protected,GetAccess = protected)
        nbDim
    end
    
    %% Constructor
    methods
        function this = LinOpCpx(sz)
            warning('LinOpCpx Deprecated');
            this.name ='LinOpCpx';
            assert(issize(sz),'The input size sz should be a conformable  to a size ');
            this.isInvertible = true;
            this.isDifferentiable=true;           
            this.sizeout = [sz 2];
            this.sizein = sz;
            this.nbDim = numel(sz);
		end
    end

    %% Core Methods containing implementations (Protected)
	methods (Access = protected)
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y = cat(this.nbDim+1, real(x),imag(x));
        end
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOp`.
            y = reshape(x,[],2);
            y = complex(y(:,1),y(:,2));
            y = reshape(y,this.sizein);
        end
        function y = applyHtH_(~,x)
            % Reimplemented from parent class :class:`LinOp`.
            y =x;
        end
        function y = applyHHt_(~,x)
            % Reimplemented from parent class :class:`LinOp`.
            y =x;
        end      
        function y = applyInverse_(~,x)
            % Reimplemented from parent class :class:`LinOp`.
            y = complex(x(:,:,1),x(:,:,2));
        end
        function M = makeHHt_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M=LinOpDiag(this.sizein);
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M=LinOpDiag(this.sizein);
        end
    end
end
