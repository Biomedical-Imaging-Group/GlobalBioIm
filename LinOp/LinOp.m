classdef LinOp < handle
    %% LinOP : Linear Operator generic class
    %  Matlab Linear Operator Library
    % The LinOp meta class implement generic methods for all linear
    % operators $\mathbf{H}$: $\mathcal{X}\rightarrow\mathcal{Y}$
    %
    % $$ \mathbf{y} = \mathbf{H} \mathbf{x} $$
    %
    %% Properties
    % * |name|          - name of the linear operator $\mathbf{H}$
    % * |sizein|        - dimension of the right hand side vector space $\mathcal{X}$
    % * |sizeout|       - dimension of the left hand side vector space $\mathcal{Y}$
    % * |isinvertible|  - true if the operator is invertible
    % * |issquare|      - true if the operator is square
    % * |iscomplex|     - true is the operator is complex
    % * |norm|          - norm of the operator (if known, otherwise -1)
    %
    %% Methods
    % * |Apply|    - Apply the operator $\mathbf{H}$
    % * |Adjoint|  - Apply the adjoint  $\mathbf{H}^{\mathrm{*}}$ defined
    % such  $<\mathbf{H} \mathbf{x} . \mathbf{y} > = <\mathbf{x} .
    % \mathbf{H}^{\mathrm{*}} \mathbf{y} >$
    % * |HtH|      - Apply the HtH matrix: the operator followed by its adjoint
    % $\mathbf{H}^{\mathrm{*}} \cdot \mathbf{H}$    
    % * |HHt|      - Apply the HHt matrix: the adjoint operator followed by its 
    % $ \mathbf{H} \cdot \mathbf{H}^{\mathrm{*}} $
    % * |Inverse|  - Apply the inverse  $\mathbf{H}^{\mathrm{-1}}$ of the
    % operator if it exist
    % * |AdjointInverse|  - Apply the adjoint of the inverse  $\mathbf{H}^{\mathrm{-*}}$ of the
    % operator if it exist
    %%
     
    
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
        name = 'none'           % name of the linear operator
        sizein;                 % dimension of the right hand side vector space
        sizeout;                % dimension of the left hand side vector space
        isinvertible = false;   % true if the operator is invertible
        issquare = false;       % true if the operator is square
        iscomplex = false;      % true is the operator is complex
        norm=-1;                % norm of the operator
    end
    
    methods
        function Apply(~,~) % Apply the operator
            error('Operator not implemented');
        end
        function Adjoint(~,~) % Apply the adjoint
            error('Adjoint not implemented');
        end
        function y = HtH(this,x) %  Apply the HtH matrix
            y = this.Adjoint(this.Apply(x));
        end
        function y = HHt(this,x) %  Apply the HHt matrix
            if this.issquare   % HtH =HHt
                y = this.HtH(x);
            else
            y = this.Apply(this.Adjoint(x));
            end
        end
        function y = HtWH(this,x,W) %  Apply the HtH matrix
            if (isscalar(W) && isreal(W))
                y = W.*this.HtH(x);
            else
                assert(isa(W,'LinOp'),'W must be a LinOp');
                y = this.Adjoint(W.Apply(this.Apply(x)));
            end
        end
        function Inverse(this,~) % Apply the inverse
            if this.isinvertible
                error('Inverse not implemented');
            else
                error('Operator not invertible');
            end
        end
        function AdjointInverse(this,~) % Apply the inverse
            if this.isinvertible
                error('AdjointInverse not implemented');
            else
                error('Operator not invertible');
            end
        end
        function this = transpose(this) % Overloading for transpose 
            if this.iscomplex
            warning('Warning: Do you mean adjoint? For LinOp object transpose() is an alias of Adjoint method');
            end
           this = adjoint(this);
        end
        function this = ctranspose(this) % Overloading for ctranspose 
           this = adjoint(this);
        end
        function y =	mtimes(this,x)% Overloading for *
            if isa(x,'LinOp')
            y = MulLinOp(this, x);
            else
                y = this.Apply(x);
            end
        end
        function y =	plus(this,x)% Overloading for +
            assert(isa(x,'LinOp'),'addition of LinOp is only define with other LinOp');
            y = SumLinOp({this,x});
		end
	end
	
	
	methods (Static, Access=protected)
		function checkSize(x, sz)
			assert( isequal(size(x),sz),  'x does not have the correct size, [%s]', num2str(sz));
			
		end
	end
end

