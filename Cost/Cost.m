classdef Cost < handle
    % Abstract class for Cost functions
    % $$ C : \\mathrm{X} \\longrightarrow \\mathbb{R}$$
    % with the following special structure
    % $$ C(\\mathrm{x}) := F( \\mathrm{Hx} , \\mathrm{y} ) $$
	% where \\(F\\) is a function takin two variables.
	%
    % :param H: a :class:`LinOp` object (default :class:`LinOpIdentity`)
    % :param y: data vector of size H.sizeout (default 0)
    % :param name: name of the cost function
    % :param lip: Lipschitz constant of the gradient (when applicable and known, otherwise -1)
	% :param isconvex: boolean true is the function is convex
    %
	% See also :class:`LinOp`.

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
        y=0;				% data y
        H=[];  				% linear operator
        isconvex=false;     % boolean
    end
    
    methods
        function f=eval(this,x)
        	% **(Abstract Method)** Evaluates the cost function
       	 	%
        	% :param x: \\(\\in \\mathrm{X}\\)
        	% :returns y: \\(= C(\\mathrm{x})\\)
           
            error('Eval not implemented');         
        end

        function g=grad(this,x)
        	% **(Abstract Method)** Evaluates the gradient of the cost function (when applicable)
        	%
        	% :param x: \\(\\in \\mathrm{X}\\)
        	% :returns y: \\(= \\nabla C(\\mathrm{x})\\)
        	
            error('Grad not implemented');
        end

        function [cost , gradient] = eval_grad(this,x)
            % Evaluates both the cost function and its gradient (when applicable)
            %
        	% **Note**: For some derived classes this function is reimplemented in a faster way than running both :meth:`eval` and :meth:`grad` successively (default).
        	%
        	% :param x: \\(\\in \\mathrm{X}\\)
        	% :returns y: \\(= \\left[C(\\mathrm{x}), \\nabla C(\\mathrm{x})\\right]\\)
        	
            cost = this.eval(x) ;
            gradient = this.grad(x) ;
        end
        
        function y=prox(this,x,alpha)
        	% **(Abstract Method)** Evaluates the proximity operator of the cost (when applicable)
            % $$ \\mathrm{prox}_{\\alpha C}(\\mathrm{x}) =  \\mathrm{arg} \\, \\mathrm{min}_{\\mathrm{u} \\in \\mathrm{X}} \\; \\frac{1}{2\\alpha} \\| \\mathrm{u} - \\mathrm{x} \\|_2^2 + C(\\mathrm{u}). $$
        	%
        	% :param x: \\(\\in \\mathrm{X}\\)
        	% :param alpha: \\(\\in \\mathbb{R}\\)
        	% :returns y: \\(= \\mathrm{prox}_{\\alpha C}(\\mathrm{x})\\)
        	
            error('Prox not implemented');
        end

        function y=prox_fench(this,x,alpha)
        	% Evaluates the proximity operator of the Fenchel Transform \\(C^*\\) (when applicable) which is computed using Moreau's identity:
        	% $$\\mathrm{prox}_{\\alpha C^*}(\\mathrm{x}) = \\mathrm{x} - \\alpha \\,\\mathrm{prox}_{\\frac{1}{\\alpha}C}\\left(\\frac{\\mathrm{x}}{\\alpha}\\right).$$
        	% **Note-1**: Only defined if :attr:`isconvex` =True.
        	%
        	% **Note-2** When defining a new class :class:`Cost`, one only requires to implements the prox.
        	%
        	% :param x: \\(\\in \\mathrm{X}\\)
        	% :param alpha: \\(\\in \\mathbb{R}\\)
        	% :returns y: \\(= \\mathrm{prox}_{\\alpha C^*}(\\mathrm{x})\\)
        	
            assert(isscalar(alpha),'alpha must be a scalar');
            if this.isconvex
                y= x - alpha*(this.prox((x-this.y)/alpha,1/alpha)+this.y);
            else
                error('Prox Fenchel not implemented');
            end
        end

        function Cnew=	o(this,L)
        	% Compose the cost \\(C\\) with a :class:`LinOp` \\(\\mathrm{L}\\)
        	% $$ C_{new}(\\mathrm{x}) := C(\\mathrm{Lx}) $$
        	%
        	% :param L: :class:`LinOp` object.
        	% :returns Cnew: the new :class:`Cost`.
        	%
        	% See also: :class:`ComposeLinOpCost`
        	
            assert(isa(L,'LinOp'),' Composition is only with a LinOp');
            Cnew=ComposeLinOpCost(this,L);
        end
        
        function Cnew = plus(this,C2)
        	% Overload operator (+) for :class:`Cost` objects
        	% $$ C_{new}(\\mathrm{x}) := C(\\mathrm{x}) + C_2(\\mathrm{x})$$
        	%
        	% :param C2: :class:`Cost` object.
        	% :returns Cnew: :class:`Cost`.
        	%
        	% See also: :class:`SumCost`
        	
            assert(isa(C2,'Cost'),'Addition is only between two Costs');
            Cnew = SumCost({this,C2});
        end

        function C = minus(this,C2)
            % Overload operator (-) for :class:`Cost` objects
        	% $$ C_{new}(\\mathrm{x}) := C(\\mathrm{x}) - C_2(\\mathrm{x})$$
        	%
        	% :param C2: :class:`Cost` object.
        	% :returns Cnew: :class:`Cost`.
        	%
        	% See also: :class:`SumCost`
        	
            assert(isa(C2,'Cost'),'Subtraction is only between two Costs');
            V = SumCost({this,C2},[1,-1]);
        end

        function y = mtimes(this,C2)
            % Overload operator (-) for :class:`Cost` objects
        	% $$ C_{new}(\\mathrm{x}) := C(\\mathrm{x}) \\times C_2(\\mathrm{x})$$
        	%
        	% :param C2: :class:`Cost` object or a scalar in \\(\\mathbb{R}\\).
        	% :returns Cnew: :class:`Cost`.
        	%
        	% See also: :class:`MultCost`
        	
            if isa(C2,'Cost')
                y = MulCost(this,C2);
            elseif isscalar(C2)
                y = MulCost(C2,this);
            else
                error('C2 must be a Cost object or a scalar');
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
