classdef Cost < Map
    % Abstract class for Cost functions
    % $$ C : \\mathrm{X} \\longrightarrow \\mathbb{R}$$
    % with the following special structure
    % $$ C(\\mathrm{x}) := F( \\mathrm{x} , \\mathrm{y} ) $$
    % where \\(F\\) is a function takin two variables, both in \\(\\mathrm{X}\\).
    %
    % :param y: data vector  (default 0)
    % :param name: name of the cost function
    % :param lip: Lipschitz constant of the gradient (when applicable and known, otherwise -1)
    % :param isConvex: true if the cost is convex 
    % :param isSeparable: true if the cost is separable (R^n basis) 
    %
    % All attributes of parent class :class:`Map` are inherited and
    % :attr:`norm` is fixed to -1, :attr:`sizeout` is fixed to   for all :class:`Cost`
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
    
    %% Write protected properties
    properties  (SetAccess = protected,GetAccess = public)
        isConvex=false;           % true if the Cost is convex
        isSeparable=false;        % true is the Cost is separable (R^n basis)
        lip=-1;                   % Lipschitz constant of the gradient
        y=0;				      % data y
    end    
    
    properties (SetAccess = protected,GetAccess = protected)
        recProxProxFench=0;      % boolean to control infinite recursion between Prox and ProxFench when not implemented in subclasses
    end
    
    %% Constructor
    methods
        function this=Cost(sz,y)
            % set data y
            if nargin <2, y=0; end
            if nargin<1 || isempty(sz), sz=size(y); end;
            assert(issize(sz),'The input size is not conformable to a size');
            this.sizein=sz;
            this.set_y(y);
            % Add new fields to memoizeOpts and memoCache
            this.memoizeOpts.applyGrad=false;
            this.memoizeOpts.applyProx=false;
            this.memoizeOpts.applyProxFench=false;
            this.memoCache.applyGrad=struct('in', [], 'out', []);
            this.memoCache.applyProx=struct('in', [], 'out', []);
            this.memoCache.applyProxFench=struct('in', [], 'out', []);
            % Properties fixed for costs
            this.sizeout=1;                % dimension of the left hand side vector space
            this.norm=-1;                  % norm of the operator
        end
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
            if this.memoizeOpts.applyGrad
                g = this.memoize('applyGrad', @this.applyGrad_, x);
            else
                g =this.applyGrad_(x);
            end
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
            assert(isnumeric(alpha),'alpha must be a numeric');
            if ~checkSize(z, this.sizein) % check input size
                error('Input to applyProx was size [%s], didn''t match stated sizein [%s].',...
                    num2str(size(z)), num2str(this.sizein));
            end
            % memoize
            if this.memoizeOpts.applyProx
                x = this.memoize('applyProx', @this.applyProx_,{z,alpha});
            else
                x =this.applyProx_(z,alpha);
            end
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
            if this.memoizeOpts.applyProxFench
                y = this.memoize('applyProxFench', @this.applyProxFench_,{z,alpha});
            else
                y =this.applyProxFench_(z,alpha);
            end
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
    % - plus_(this,G)
    % - minus_(this,G)
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
            if this.isConvex && ~this.recProxProxFench
                this.recProxProxFench=1;
                x= z - alpha*(this.applyProxFench(z/alpha,1/alpha));
                this.recProxProxFench=0;
            else
                error('applyProx_ and applyProxFench_ methods not implemented');
            end
        end
        function y=applyProxFench_(this,z,alpha)
            % By default, if the cost \\(C\\) :attr:`isConvex`, computes the proximity operator of the Fenchel Transform
            % \\(C^*\\) at \\(\\mathrm{z} \\in \\mathrm{Y} \\) using the
            % Moreau's identity which uses the :meth:`applyProx` method
            % $$\\mathrm{prox}_{\\alpha C^*}(\\mathrm{z}) = \\mathrm{z} - \\alpha \\,\\mathrm{prox}_{\\frac{1}{\\alpha}C}\\left(\\frac{\\mathrm{z}}{\\alpha}\\right).$$
            if this.isConvex && ~this.recProxProxFench
                this.recProxProxFench=1;
                y= z - alpha*(this.applyProx(z/alpha,1/alpha));
                this.recProxProxFench=0;
            else
                error('applyProx_ and applyProxFench_  methods not implemented');
            end
        end
        function M = plus_(this,G)
            % If \\(\\mathrm{G}\\) is a :class:`Cost`, constructs a :class:`CostSummation` object to sum the
            % current :class:`Cost` \\(C\\) with the given \\(G\\).
            % Otherwise the summation will be a :class:`MapSummation`.
            if isa(G,'Cost')
                M = CostSummation({this,G},[1,1]);
            else
                M = plus_@MapSummation({this,G},[1,1]);
            end
        end
        function M = minus_(this,G)
            % If \\(\\mathrm{G}\\) is a :class:`Cost`, constructs a :class:`CostSummation` object to subtract to the
            % current :class:`Cost` \\(C\\), the given \\(G\\).
            % Otherwise the summation will be a :class:`MapSummation`.
            if isa(G,'Cost')
                M = CostSummation({this,G},[1,-1]);
            else
                M = minus_@MapSummation({this,G},[1,-1]);
            end
        end
        function M=makeComposition_(this,G)
            % Reimplemented from parent class :class:`Map`.
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
            x=y.*this.applyGrad(v);
        end
    end
    
    %% Utility methods
    % - set_y(this,y)
    methods
        function set_y(this,y)
            % Set the attribute \\(\\mathrm{y}\\)
            %
            %  - has to be conformable with the :attr:`sizein` of the cost
            %  - can be a scalar.
            assert(isnumeric(y),'y must be a numeric');
            assert(isscalar(y) || checkSize(y,this.sizein),'The size of y must be equal to the input size of the Cost');
            this.y=y;
        end
    end
    
    methods (Access = protected)
        %% Copy
        function this = copyElement(obj)
            this = copyElement@Map(obj);
            this.memoCache.applyGrad=struct('in', [], 'out', []);
            this.memoCache.applyProx=struct('in', [], 'out', []);
            this.memoCache.applyProxFench=struct('in', [], 'out', []);
        end
    end  
end
