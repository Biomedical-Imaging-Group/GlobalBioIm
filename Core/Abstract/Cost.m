classdef Cost < handle
    % Abstract class for Cost functions
    % $$ C : \\mathrm{X} \\longrightarrow \\mathbb{R}$$
    % with the following special structure
    % $$ C(\\mathrm{x}) := F( \\mathrm{Hx} , \\mathrm{y} ) $$
    % where \\(F\\) is a function takin two variables.
    %
    % :param H: a :class:`Map` object (default :class:`LinOpIdentity`)
    % :param y: data vector of size H.sizeout (default 0)
    % :param name: name of the cost function
    % :param lip: Lipschitz constant of the gradient (when applicable and known, otherwise -1)
    % :param isConvex: boolean true is the function is convex
    %
    % See also :class:`Map`, :class:`LinOp`.
    
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
        isConvex=false;     % boolean
    end
    
    %% Interface Methods (cannot be overloaded in derived classes: Sealed)
    % - apply(this,x)
    % - applyGrad(this,x)
    % - applyProx(this,z,alpha)
    % - makeMapComposition(this,M)
    methods (Sealed)
        function y = apply(this, x)
            % Computes \\(\\mathrm{y}=C(\\mathrm{x})\\) for the given
            % \\(\\mathrm{x} \\in \\mathrm{X}\\).
            if ~checkSize(x, this.H.sizein) % check input size
                error('Input to apply was size [%s], didn''t match stated sizein of this.H: [%s].',...
                    num2str(size(x)), num2str(this.H.sizein));
            end
            y = this.apply_(x);
        end
        function g = applyGrad(this,x)
            % Computes the gradient of the cost function at  \\(\\mathrm{x} \\in \\mathrm{X}\\) (when applicable)
            % $$ \\mathrm{g} = \\nabla C(\\mathrm{x}) $$
            if ~checkSize(x, this.H.sizein) % check input size
                error('Input to applyGrad was size [%s], didn''t match stated sizein of this.H: [%s].',...
                    num2str(size(x)), num2str(this.H.sizein));
            end
            g = this.applyGrad_(x);
        end
        function x=applyProx(this,z,alpha)
            % Computes the proximity operator of the cost at \\(\\mathrm{z} \\in \\mathrm{X} \\) (when applicable)
            % $$ \\mathrm{prox}_{\\alpha C}(\\mathrm{z}) =  \\mathrm{arg} \\, \\mathrm{min}_{\\mathrm{u} \\in \\mathrm{X}} \\; \\frac{1}{2\\alpha} \\| \\mathrm{u} - \\mathrm{z} \\|_2^2 + C(\\mathrm{u}). $$
            assert(isscalar(alpha),'alpha must be a scalar');
            if ~checkSize(z, this.H.sizein) % check input size
                error('Input to applyProx was size [%s], didn''t match stated sizein of this.H: [%s].',...
                    num2str(size(z)), num2str(this.H.sizein));
            end
            x = this.applyProx_(z,alpha);
        end
        function y=prox_fench(this,z,alpha)
            % Computes the proximity operator of the Fenchel Transform \\(C^*\\) at \\(\\mathrm{z} \\in \\mathrm{Y} \\) (when applicable)
            assert(isscalar(alpha),'alpha must be a scalar');
            if ~checkSize(z, this.H.sizeout) % check input size
                error('Input to applyProxFench was size [%s], didn''t match stated sizeout of this.H: [%s].',...
                    num2str(size(z)), num2str(this.H.sizeout));
            end
            y = this.applyProxFench_(z,alpha);
        end
        function Cnew=	makeMapComposition(this,M)
            % Compose the cost \\(C\\) with a :class:`Map` \\(\\mathrm{M}\\)
            % $$ C_{new}(\\mathrm{x}) := C(\\mathrm{Mx}) $$
            % See also: :class:`ComposeMapCost`
            assert(isa(M,'Map'),'Composition is only with a Map');
            assert(isequal(M.sizeout,this.H.sizein),'The Map M to compose has sizeout inconsistent with the sizein of this.H');
            Cnew=this.makeMapComposition_(M);
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyGrad_(this,x)
    % - applyProx_(this,z,alpha)
    % - applyProxFench_(this,z,alpha)
    % - makeMapComposition_(this,M)
    methods (Access = protected)
        function y = apply_(this, x)
            error('apply_ method not implemented');
        end
        function g=applyGrad_(this,x)
            error('applyGrad_ method not implemented');
        end
        function x=applyProx_(this,z,alpha)
            % By default, if the cost \\(C\\) :attr:`isConvex`, computes the proximity operator of  \\(C^*\\)
            % at \\(\\mathrm{z} \\in \\mathrm{X} \\) using the
            % Moreau's identity which uses the :meth:`applyProxFench` method
            % $$\\mathrm{prox}_{\\alpha C}(\\mathrm{z}) = \\mathrm{z} - \\alpha \\,\\mathrm{prox}_{\\frac{1}{\\alpha}C^*}\\left(\\frac{\\mathrm{z}}{\\alpha}\\right).$$
            if this.isConvex
                x= z - alpha*(this.applyProxFench(z/alpha,1/alpha));
            else
                error('applyProx_ method not implemented');
            end
        end
        function y=applyProxFench_(this,z,alpha)
            % By default, if the cost \\(C\\) :attr:`isConvex`, computes the proximity operator of the Fenchel Transform
            % \\(C^*\\) at \\(\\mathrm{z} \\in \\mathrm{Y} \\) using the
            % Moreau's identity which uses the :meth:`applyProx` method
            % $$\\mathrm{prox}_{\\alpha C^*}(\\mathrm{z}) = \\mathrm{z} - \\alpha \\,\\mathrm{prox}_{\\frac{1}{\\alpha}C}\\left(\\frac{\\mathrm{z}}{\\alpha}\\right).$$
            if this.isConvex
                y= z - alpha*(this.applyProx(z/alpha,1/alpha));
            else
                error('applyProxFench not implemented');
            end
        end
        function Cnew=	makeMapComposition_(this,M)
            Cnew=this;
            Cnew.set_H(this.H.makecomposition(M));
        end
    end
    
    %% Special Methods (mainly overload of Matlab operators)
    methods
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
            C = SumCost({this,C2},[1,-1]);
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
    end
    
    
    %% Utility methods
    % - set_y(this,y)
    % - set_H(this,H)
    methods
        function set_y(this,y)
            % Set the attribute \\(\\mathrm{y}\\) by checking the size
            % of \\(\\mathrm{y}\\) which
            %   - has to be conformable with the :attr:`sizeout` of the
            %     :class:`Map`\\(\\mathrm{H}\\),
            %   - can be anything if \\(\\mathrm{H}\\) is not yet set (empty),
            %   - can be a scalar.
            assert(isnumeric(y),'y must be a numeric');
            assert(isempty(this.H) || ((isscalar(y)) || checkSize(y,this.H.sizeout)),'size y must be equal to this.H.sizeout');
            this.y=y;
        end
        function set_H(this,H)
            % Set the :class:`Map` \\(\\mathrm{H}\\) by checking its type and size which can be
            %
            %  - a :class:`Map` with :attr:`sizeout` conformable to the size
            %    of the attribute \\(\\mathrm{y}\\) (when \\(\\mathrm{y}\\) is not a scalar),
            %  - a size conformable to the size of the attribute \\(\\mathrm{y}\\)
            %    (when \\(\\mathrm{y}\\) is not a scalar) to create a
            %    :class:`LinOpIdentity`,
            %  - empty to create a :class:`LinOpIdentity` with the size of \\(\\mathrm{y}\\)
            %
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
