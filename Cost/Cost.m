classdef Cost < handle
    %% Cost : Cost function generic class
    %  Matlab Inverse Problems Library
    %  The Cost meta class implements generic methods for all cost function
    %  Cost functions  $\mathbb{C}^N\rightarrow\mathbb{R}$
    %                                 $cost \rightarrow f( H.x , y )$
    % return a real scalar  f( H.x , y )
    % -- Properties
    % * |name|       - name of the function
    % * |H|          - LinOp composed with the functional if H is a size it
    % implicitly make H=LinOpIdentity(H);
    % * |y|          - a vector of size H.sizeout
    % * |lip|        - Lipschitz constant of the gradient (if known, otherwise -1)
    % * |isconvex|   - boolean true is the function is convex
    %
    % -- Methods
    % * |eval|       - evaluates the functional
    % * |grad|       - evaluates the gradient of the functional
    % * |o|          - compose with a LinOp
    % * |prox|       - computes the proximity operator
    % * |prox_fench| - computes the proximity operator of the fenchel transform
    %                  (default for convex Cost: uses the Moreau's identity
    %                       prox_{alpha F*}(y) = y - alpha prox_{F/alpha}(y/alpha)
    %
    %     Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch & F. Soulez ferreol.soulez@epfl.ch
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
    
    % Protected Set and public Read properties
    properties (SetAccess = protected,GetAccess = public)
        name = 'none'       % name of the functional
        lip=-1;             % Lipschitz constant of the gradient
        % LinOp Infos
        y=0;
        H=[];  % linear operator
        isconvex=false;
    end
    
    methods
        %% Evaluation of the Functional
        function eval(~,~)
            error('Eval not implemented');
        end
        %% Gradient of the Functional
        function grad(~,~)
            error('Grad not implemented');
        end
        %% Compute the functional and its gradient
        function [cost , gradient] = eval_grad(this,x)
            cost = this.eval(x) ;
            gradient = this.grad(x) ;
        end
        %% Proximity operator of the functional
        function prox(~,~,~)
            error('Prox not implemented');
        end
        %% Proximity operator of the Fenchel transform of the functional prox_{alpha F*}
        function y=prox_fench(this,x,alpha)
            assert(isscalar(alpha),'alpha must be a scalar');
            if this.isconvex
                y= x - alpha*(this.prox((x-this.y)/alpha,1/alpha)+this.y);
            else
                error('Prox Fenchel not implemented');
            end
        end
        %% Operator compose with a LinOp
        function v=	o(this,x)
            assert(isa(x,'LinOp'),' Composition is only with a LinOp');
            v=ComposeLinOpCost(this,x);
        end
        %% Overload the operator +
        function y = plus(this,x)
            assert(isa(x,'Cost'),'Addition is only between two Costs');
            y = SumCost({this,x});
        end
        %% Overload the operator -
        function y = minus(this,x)
            assert(isa(x,'Cost'),'Subtraction is only between two Costs');
            y = SumCost({this,x},[1,-1]);
        end
        %% Overload the operator *
        function y = mtimes(this,x)
            if isa(x,'Cost')
                y = MulCost(this,x);
            elseif isscalar(x)
                y = MulCost(x,this);
            else
                error('x must be a Cost object or a scalar');
            end
        end
        %% Set data y (must be conformable with H sizeout)
        function set_y(this,y)
            assert(isnumeric(y),' y must be a numeric');
            assert(isempty(this.H) || ((isscalar(y)) || (isequal(this.H.sizeout,size(y)))),'size y must be equal to this.H.sizeout');
            this.y=y;
        end
        %% Set operator H (sizein of H must be conformable with size of data y)
        % H can be a LinOp, a size (to define identity) or empty (identity
        % will be defined with the size of this.y)
        function set_H(this,H)
            if isa(H, 'LinOp')
                assert( (isscalar(this.y)) || (isequal(H.sizeout,size(this.y))),'H.sizeout must be equal to size of this.y');
            elseif issize(H)
                assert( (isscalar(this.y)) || (isequal(H,size(this.y))),'the given size must be equal to the size of this.y');
                H=LinOpIdentity(H);
            elseif isempty(H)
                H =LinOpIdentity(size(this.y));
            end
            this.H=H;
        end
    end
end
