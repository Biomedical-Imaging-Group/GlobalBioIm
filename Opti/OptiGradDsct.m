classdef OptiGradDsct < Opti
    % Gradient Descent optimization algorithm to minimize a differentiable :class:`Cost` \\(C(\\mathrm{x})\\)
    %
    % :param C: a differentiable :class:`Cost` (i.e. with an implementation of :meth:`applyGrad`).
    % :param gam: descent step
    %
    % All attributes of parent class :class:`Opti` are inherited. 
    %
    % **Note** If the cost \\(C\\) is gradient Lipschitz, convergence is ensured by taking 
    % \\(\\gamma \\in (0,2/L] \\) where \\(L\\) is the Lipschitz constant of \\(\\nabla C\\) (see [1]).
	% The optimal choice is \\(\\gamma = 1/L \\) (see [1]). If \\(L\\) is known (i.e. F.lip different from -1), 
    % parameter \\(\\gamma\\) is automatically set to \\(1/L\\).
    %
    % **Reference**
    %
    % [1] Nesterov, Yurii. "Introductory lectures on convex programming." Lecture Notes (1998): 119-120.
    %
    % **Example** GD=OptiGradDsct(F)
    %
    % See also :class:`Opti` :class:`OutputOpti` :class:`Cost`
    
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
    	gam=[];      % descent step
    end
    
    %% Constructor
    methods   	
    	function this=OptiGradDsct(F)
            % Set properties
    		this.name='Opti Gradient Descent';
    		this.cost=F;
            % Initialize
            this.initObject('OptiGradDsct');
        end
    end
    %% updateProp method (Private)
    methods (Access = protected)
        function updateProp(this,prop)
            % Reimplemented superclass :class:`Opti`
            
            % Call superclass method
            updateProp@Opti(this,prop);
            % Update current-class specific properties
            if strcmp(prop,'cost') || strcmp(prop,'all')
                if this.cost.lip>0 && (isempty(this.gam) || this.gam >this.cost.lip)
                    this.gam=1/this.cost.lip;
                end
            end
        end
    end
    
    %% Methods for optimization
    methods
        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.
            
            initialize@Opti(this,x0);
            if isempty(this.gam), error('Parameter gam is not setted'); end
        end
        function flag=doIteration(this)
            % Reimplementation from :class:`Opti`.  Performs:
            % $$ \\mathrm{x}^{k+1} = \\mathrm{x}^k - \\gamma \\nabla C(\\mathrm{x}^k) $$
            
            this.xopt=this.xopt-this.gam*this.cost.applyGrad(this.xopt);
            flag=this.OPTI_NEXT_IT;
        end
    end
end
