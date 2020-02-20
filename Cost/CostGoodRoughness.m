classdef CostGoodRoughness < Cost
    % CostGoodRoughness: see [1]
    % $$C(\\mathrm{x}) := \\sum_k \\frac{\\vert (\\nabla \\mathrm{x})_{.,k} \\vert^2}{\\sqrt{\\vert \\mathrm{x} \\vert^2 + \\beta}} $$
    % with \\( \\vert (\\nabla \\mathrm{x})_{.,k} \\vert^2 \\)  the gradient
    % magnitude at pixel k.
    %
    % :param G: :class:`LinOpGrad` object
    % :param bet: smoothing parameter (default 1e-1)
    %
    % All attributes of parent class :class:`Cost` are inherited. 
    %
    % **Example** GR=CostGoodRoughness(G)
    %
    % **Reference**
    % [1] Verveer, P. J., Gemkow, M. J., & Jovin, T. M. (1999). A comparison of image restoration
    % approaches applied to three?dimensional confocal and wide?field fluorescence microscopy. 
    % Journal of microscopy, 193(1), 50-61.
    %
    % See also :class:`Map`, :class:`Cost`, :class:`LinOpGrad`

    %%    Copyright (C) 2018
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
    

    properties 
        bet; % smoothing parameter
    end
    % Protected Set properties
    properties (SetAccess = protected,GetAccess = public)
        G;        % gradient operators
        sumOp;
    end

    %% Constructor
    methods
        function this = CostGoodRoughness(G,bet)
            this@Cost(G.sizein);            
            this.name='CostGoodRoughness';
            if nargin<2|| isempty(bet)
                bet=1e-1;
            end
            this.bet=bet;
            this.isConvex=false;
            this.isDifferentiable=true;
            this.G=G;
            this.sumOp=LinOpSum(G.sizeout,this.G.ndms+1);
            
            this.memoizeOpts.computeNum=true;
            this.memoCache.computeNum= struct('in',[],'out', []);
            this.memoizeOpts.computeDem=true;
            this.memoCache.computeDem= struct('in',[],'out', []);
        end        
    end
    
    %% Core Methods containing implementations (Protected)
    % - apply_(this,x)
    % - applyGrad_(this,x)
	methods (Access = protected)
        function y=apply_(this,x)
        	% Reimplemented from parent class :class:`Cost`. 
            num = this.memoize('computeNum',@this.computeNum_,x);
            dem = this.memoize('computeDem',@this.computeDem_,x);
            y=sum(num(:)./sqrt(dem(:)));
        end
        function g=applyGrad_(this,x)
        	% Reimplemented from parent class :class:`Cost`.
            num = this.memoize('computeNum',@this.computeNum_,x);
            dem = this.memoize('computeDem',@this.computeDem_,x);           
            gnum=this.G.applyAdjoint(2*this.G*x); % derivative cost numerator
            g=gnum./sqrt(dem)-num./dem.^(3/2).*x;
        end
        function num = computeNum_(this,x)
            gx=this.G*x;                  % gradient of x
            num=this.sumOp*abs(gx).^2;    % cost numerator
        end
        function dem = computeDem_(this,x)
            dem=abs(x).^2 + this.bet;                % squared denominator gdem=dem.^(-3/2)
        end
    end
    
    methods (Access = protected)
         %% Copy
         function this = copyElement(obj)
             this = copyElement@Cost(obj);
             this.G = copy(obj.G);
            this.sumOp = copy(obj.sumOp);          
      end
    end  
end
