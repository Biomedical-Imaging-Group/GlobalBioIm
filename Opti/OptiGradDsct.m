classdef OptiGradDsct < Opti
    % Gradient Descent optimization algorithm to minimize a differentiable :class:`Cost` \\(C(\\mathrm{x})\\)
    %
    % :param C: a differentiable :class:`Cost` (i.e. with an implementation of :meth:`applyGrad`).
    % :param gam: descent step
    % :param nagd: boolean (default false) to activate the Nesterov accelerated gradient descent
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

	% Full public properties
    properties
        nagd=false;      % Nesterov accelerated gradient descent
    	gam=[];      % descent step
        y;
    end
    
    methods
        %% Constructor
        function this=OptiGradDsct(F)
            this.name='Opti Gradient Descent';
            this.cost=F;
            if F.lip>0
                this.gam=1/F.lip;
            end
        end
        %% Run the algorithm
        function initialize(this,x0)
            % Reimplementation from :class:`Opti`.
            
            initialize@Opti(this,x0);
            if this.nagd
                this.needxold = true;
                this.xold = x0;
                this.y=x0;
            end
            if isempty(this.gam), error('Parameter gam is not setted'); end
        end
        function flag=doIteration(this)
            % Reimplementation from :class:`Opti`.  Performs:
            % $$ \\mathrm{x}^{k+1} = \\mathrm{x}^k - \\gamma \\nabla C(\\mathrm{x}^k) $$
            
            if this.nagd
                this.xopt=this.y - this.gam.*this.cost.applyGrad(this.y);
                this.y = this.xopt +  (this.niter -1)/(this.niter+2)*(this.xopt - this.xold);
            else
                this.xopt=this.xopt-this.gam*this.cost.applyGrad(this.xopt);
            end
            
           
            flag=this.OPTI_NEXT_IT;
        end
	end
end
