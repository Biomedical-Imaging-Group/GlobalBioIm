classdef CostKullLeib < Cost
    % CostKullLeib: KullbackLeibler divergence
    % $$ C(\\mathrm{x}) :=\\sum_n D_{KL}(\\mathrm{x}_n)$$
    % where
    % $$ D_{KL}(\\mathrm{z}_n) := \\left\\lbrace \\begin{array}[ll]
    % \\mathrm{z}_n - \\mathrm{y}_n \\log(\\mathrm{z}_n + \\beta) & \\text{ if } \\mathrm{z}_n + \\beta >0  \\newline
    % + \\infty &  \\text{otherwise}.
    % \\end{array} \\right.$$
    %
    % All attributes of parent class :class:`Cost` are inherited. 
    %
    % :param bet: smoothing parameter \\(\\beta\\) (default 0) 
    %
    % **Example** C=CostKullLeib(sz,y,bet)
    %
    % See also :class:`Map` :class:`Cost`, :class:`LinOp`

    %%    Copyright (C) 2017 
    %     E. Soubies emmanuel.soubies@epfl.ch
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
    
    %% Properties
    % - Public
    properties (SetObservable, AbortSet)
        bet;     % smoothing parameter, if bet=0 then the unsmoothed version is used
    end
    
    %% Constructor and handlePropEvents method
    methods       
        function this = CostKullLeib(sz,y,bet)
            % Default values
            if nargin<2 || isempty(y), y=0; end
            if nargin<3, bet=0; end
            % Call superclass constructor
            this@Cost(sz,y);
            % Listeners to PostSet events
            addlistener(this,'bet','PostSet',@this.handlePropEvents);
            % Set properties
            this.name='CostKullLeib';
            this.bet=bet;
            this.isConvex=true;
            this.isSeparable=true;
        end
        function handlePropEvents(this,src,~)
            % Reimplemented from superclass :class:`Cost`
            switch src.Name
                case {'bet','y'}
                    % Compute Lipschitz constant of the gradient
                    if (this.bet>0)
                        this.lip=max(this.y(:))./this.bet^2;
                        this.isDifferentiable=true;
                    else
                        this.isDifferentiable=false;
                    end
            end
            % Call superclass method (important to ensure the right execution order)
            handlePropEvents@Cost(this,src);
        end
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyGrad_(this,x)
    % - applyProx_(this,x,alpha)
	methods (Access = protected)
        function f=apply_(this,x)
        	% Reimplemented from parent class :class:`Cost`.
        	
            if any(x(:)<0)
                f=Inf;
            else
                if (this.bet~=0)
                    f=sum(-this.y(:).*log(x(:)+this.bet) + x(:));
                else
                    ft = zeros_(this.sz);
                    zidx = (x~=0);
                    ft(zidx)=-this.y(zidx).*log(x(zidx)) + x(zidx);
                    f=sum(ft(:));
                end
            end
        end
        function g=applyGrad_(this,x)
        	% Reimplemented from parent class :class:`Cost`.
            
            g= 1 - this.y./(x+this.bet);
        end  
        function z=applyProx_(this,x,alpha)
            % Reimplemented from parent class :class:`Cost`.
            
            if (this.bet~=0)
                delta=(x-alpha-this.bet).^2+4*(x*this.bet + alpha*(this.y-this.bet));
                z=zeros_(size(x));
                mask=delta>=0;
                z(mask)=0.5*(x(mask)-alpha-this.bet + sqrt(delta(mask)));
            else
                z=applyProx_@Cost(this,x,alpha);
            end
        end
    end
end
