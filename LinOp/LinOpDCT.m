classdef LinOpDCT <  LinOp
    % LinOpDCT : Discrete cosine transform
    %
    % :param sz: sizein of the operator.
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Note** Currently only implemented in 2D
    %
    % **Example** DCT=LinOpDCT(sz)
    %
    % See also :class:`LinOp`, :class:`Map`, dct2, idct2, 

    %%    Copyright (C) 2021 
    %     E. Soubies emmanuel.soubies@irit.fr
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
        index % index along wich dimension are computed the FFT
        Notindex% ~index
        N % Number of element
        ndms % number of dimensions
        unitary
    end
    
    %% Constructor
    methods
        function this = LinOpDCT(sz)
            assert(issize(sz),'The input size sz should be a conformable  to a size ');
            this.sizein=sz;
            this.sizeout=sz;
            this.name ='LinOpDCT';
            this.isDifferentiable=true;            
            this.ndms = length(this.sizein);
            this.norm=1;              
        end
    end
    
    %% Core Methods containing implementations (Protected)
	methods (Access = protected)
        function y = apply_(~,x)
            % Reimplemented from parent class :class:`LinOp`.
            y=dct2(x);
        end
        function y = applyAdjoint_(~,x)
            % Reimplemented from parent class :class:`LinOp`.
            y=idct2(x);
        end
        function y = applyInverse_(~,x)
            % Reimplemented from parent class :class:`LinOp`.
            y=idct2(x);
        end
        function y = applyHtH_(~,x)
            % Reimplemented from parent class :class:`LinOp`.
            y = x;
        end
        function y = applyHHt_(~,x)
            % Reimplemented from parent class :class:`LinOp`.
            y=x;
        end
        function y = applyAdjointInverse_(~,x)
            % Reimplemented from parent class :class:`LinOp`. 
            y=x;
        end
        function M = makeHHt_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M=LinOpDiag(this.sizein);
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M=LinOpDiag(this.sizein);
        end
        function M = makeInversion_(this)
            % Reimplemented from parent class :class:`LinOp`.
            M = this.makeAdjoint;
        end
    end
end

