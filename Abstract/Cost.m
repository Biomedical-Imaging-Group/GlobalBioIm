classdef Cost < Map
    % Abstract class for Cost functions
    % $$ C : \\mathrm{X} \\longrightarrow \\mathbb{R}$$
    % with the following special structure
    % $$ C(\\mathrm{x}) := F( \\mathrm{x} , \\mathrm{y} ) $$
    % where \\(F\\) is a function takin two variables, both in \\(\\mathrm{X}\\).
    %
    % :param y: data vector of size H.sizeout (default 0)
    % :param name: name of the cost function
    % :param lip: Lipschitz constant of the gradient (when applicable and known, otherwise -1)
    % :param isConvex: true if the cost is convex 
    %
    % All attributes of parent class :class:`Map` are inherited and
    % :attr:`norm` is fixed to -1, :attr:`sizeout`is fixed to 1, and
    % isComplexOut is fixed to false for all :class:`Cost`
    %
    % See also :class:`Map`, :class:`LinOp`.
    
    %%    Copyright (C) 2017 
    %     E. Soubies emmanuel.soubies@epfl.ch 
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
    
    %% Public properties
    properties 
        name = 'none'             % name of the functional
        isConvex=false;           % true if the Cost is convex
        lip=-1;                   % Lipschitz constant of the gradient
        y=0;				      % data y
        memoizeOpts = struct('apply', false, ...
                             'applyJacobianT', false, ...
                             'applyInverse', false,...
                             'applyGrad',false,...
                             'applyProx',false,...
                             'applyProxFench',false);
    end
    
    %% Properties inherited from Map that are fixed in Cost 
    properties (SetAccess = private,GetAccess = private)
        isComplexOut = false;     % true if the space Y is complex valued
        sizeout=1;                % dimension of the left hand side vector space   
        norm=-1;                  % norm of the operator
        memoCache = struct('apply', struct('in', [], 'out', []),...
                           'applyJacobianT', struct('in', [], 'out', []), ...
                           'applyInverse', struct('in', [], 'out', []), ...
                           'applyGrad', struct('in', [], 'out', []),...
                           'applyProx', struct('in', [], 'out', []),...
                           'applyProxFench', struct('in', [], 'out', []));
    end
     
    %% Interface Methods (cannot be overloaded in derived classes: Sealed)
    % In addition to inherited methods from Map
    % - applyGrad(this,x)
    % - applyProx(this,z,alpha)
    % - applyProxFench(this,z,alpha)
    methods (Sealed)
        function g =applyGrad(this,x)
            % Computes the gradient of the cost function at  \\(\\mathrm{x} \\in \\mathrm{X}\\) (when applicable)
            % $$ \\mathrm{g} = \\nabla C(\\mathrm{x}) $$
            %
            % Calls the method :meth:`applyGrad_`
            if ~checkSize(x, this.sizein) % check input size
                error('Input to applyGrad was size [%s], didn''t match stated sizein [%s].',...
                    num2str(size(x)), num2str(this.sizein));
            end
            % memoize
            g = this.memoize('applyGrad', @this.applyGrad_, x); 
            if ~checkSize(g, this.sizein) % check output size
                error('Output of applyGrad was size [%s], didn''t match stated sizein [%s].',...
                    num2str(size(g)), num2str(this.sizein));
            end
        end
        function x=applyProx(this,z,alpha)
            % Computes the proximity operator of the cost at \\(\\mathrm{z} \\in \\mathrm{X} \\) (when applicable)
            % $$ \\mathrm{prox}_{\\alpha C}(\\mathrm{z}) =  \\mathrm{arg} \\, \\mathrm{min}_{\\mathrm{u} \\in \\mathrm{X}} \\; \\frac{1}{2\\alpha} \\| \\mathrm{u} - \\mathrm{z} \\|_2^2 + C(\\mathrm{u}). $$
            %
            % Calls the method :meth:`applyProx_`
            assert(isscalar(alpha),'alpha must be a scalar');
            if ~checkSize(z, this.sizein) % check input size
                error('Input to applyProx was size [%s], didn''t match stated sizein [%s].',...
                    num2str(size(z)), num2str(this.sizein));
            end
            % memoize
            x = this.memoize('applyProx', @this.applyProx_,{z,alpha}); 
            if ~checkSize(x, this.sizein) % check output size
                error('Output of applyProx was size [%s], didn''t match stated sizein [%s].',...
                    num2str(size(x)), num2str(this.sizein));
            end
        end
        function y=applyProxFench(this,z,alpha)
            % Computes the proximity operator of the Fenchel Transform \\(C^*\\) at \\(\\mathrm{z} \\in \\mathrm{Y} \\) (when applicable)
            %
            % Calls the method :meth:`applyProxFench_`
            assert(isscalar(alpha),'alpha must be a scalar');
            if ~checkSize(z, this.sizein) % check input size
                error('Input to applyProxFench was size [%s], didn''t match stated sizein [%s].',...
                    num2str(size(z)), num2str(this.sizein));
            end
            % memoize
            y = this.memoize('applyProxFench', @this.applyProxFench_,{z,alpha});
            if ~checkSize(y, this.sizein) % check output size
                error('Output of applyProxFench was size [%s], didn''t match stated sizein [%s].',...
                    num2str(size(y)), num2str(this.sizein));
            end
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - applyGrad_(this,x)
    % - applyProx_(this,z,alpha)
    % - applyProxFench_(this,z,alpha)
    % - makeComposition_(this,G)
    methods (Access = protected)
        function g=applyGrad_(this,x)
            % Not implemented in this Abstract class
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
                error('applyProxFench_ method not implemented');
            end
        end
        function M=makeComposition_(this,G)
            % Constructs a :class:`CostComposition` object to compose the
            % current Cost (this) with the given :class:`Map`\\(\\mathrm{G}\\). 
            M = CostComposition(this,G);
        end
    end  
    
    %% Methods of superclass Map that do not need to be reimplemented in derived Costs
    % - applyJacobianT_(this, y, v)
    methods (Access = protected, Sealed)
        function x = applyJacobianT_(this, y, v)
            % Uses the method applyGrad (hence do not need to be
            % reimplemented in derived classes)
            x=y*this.applyGrad(v);
        end
    end
    
    %% Utility methods
    % - set_y(this,y)
    methods
        function set_y(this,y)
            % Set the attribute \\(\\mathrm{y}\\)
            %   - has to be conformable with the :attr:`sizeout` of the
            %     :class:`Map`\\(\\mathrm{H}\\),
            %   - can be anything if \\(\\mathrm{H}\\) is not yet set (empty),
            %   - can be a scalar.
            assert(isnumeric(y),'y must be a numeric');
            assert(isscalar(y) || checkSize(y,this.sizein),'size y must be a scalar or equal to this.sizein');
            this.y=y;
        end
    end
end
