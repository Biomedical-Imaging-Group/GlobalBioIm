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
    % * |iscomplex|     - true is the operator is complex
    %
    %% Methods
    % * |Apply|    - Apply the operator $\mathbf{H}$
    % * |Adjoint|  - Apply the adjoint  $\mathbf{H}^{\mathrm{*}}$ defined
    % such  $<\mathbf{H} \mathbf{x} . \mathbf{y} > = <\mathbf{x} .
    % \mathbf{H}^{\mathrm{*}} \mathbf{y} >$
    % * |Gram|      - Apply the Gram matrix: the operator followed by its adjoint
    % $\mathbf{H}^{\mathrm{*}} \cdot \mathbf{H}$
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
        name = 'none'   % name of the linear operator
        sizein;         % dimension of the right hand side vector space
        sizeout;        % dimension of the left hand side vector space
        isinvertible = false; % true if the operator is invertible
        iscomplex;      % true is the operator is complex
    end
    
    methods
        function Apply(~,~) % Apply the operator
            error('Operator not implemented');
        end
        function Adjoint(~,~) % Apply the adjoint
            error('Adjoint not implemented');
        end
        function y = Gram(this,x) %  Apply the Gram matrix
            y = this.Adjoint(this.Apply(x));
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
    end
end

