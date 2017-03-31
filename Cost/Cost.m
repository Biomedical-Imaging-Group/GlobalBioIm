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
            error('Prox not implemented');
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
            assert(isa(x,'LinOp'),' Composition of Cost(.o) is only define with a LinOp');
            v=ComposeLinOpCost(this,x);
        end
        %% Overload the operator +
        function y = plus(this,x)
            assert(isa(x,'Cost'),'Addition of Cost is only define with other Cost');
            y = SumCost({this,x});
        end
        %% Overload the operator -
        function y = minus(this,x)
            assert(isa(x,'Cost'),'Subtraction of Cost is only define with other Cost');
            y = SumCost({this,x},[1,-1]);
        end
        %% Function that set properly the operator H (has to be modified if new properties is???H are added)
        function set_H(this,H,y)
  
            if (nargin <3 || isempty(y))
                y =0;
            end
            
            assert(isnumeric(y),' y must be a numeric');
            
            if isa(H, 'LinOp')  
                assert( (isscalar(y)) || (isequal(H.sizeout,size(y))),'y must be equal to H.sizeout');
            else if issize(H)
                H=LinOpIdentity(H);
                assert( (isscalar(y)) || (isequal(H.sizeout,size(y))),'y must be equal to H.sizeout');
                else if isempty(H)
                        H =LinOpIdentity(size(y));
                    end
                end
            end
            
            this.H=H;
            this.y = y;
            
            
        end
    end
end
